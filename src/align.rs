use anyhow::bail;

use crate::{
    ceil_div, cigar::EditOp, config::AlignConfig, equal::EqualityDefinition, mode::AlignMode,
    peq::build_peq_table, task::AlignTask,
};

/// Alias for single u64 bitvec word.
pub type Word = u64;

/// Size of word in bits.
pub const WORD_SIZE: u32 = Word::BITS;

/// Word for 000...1
pub const WORD_1: Word = 1;

/// Word bit mask. 100..00
pub const HIGH_BIT_MASK: Word = WORD_1 << (WORD_SIZE - 1);

/// Max chars we expect. Note is ASCII only.
pub const MAX_UCHAR: usize = 256;

/// Sequence alignment.
///
/// Omitted fields from original implementation.
/// * `status` removed as unnecessary in Rust.
/// * `numLocations` not necessary as Rust stores container metadata.
/// * `alignmentLength` is omitted for same reason as `numLocations`
#[derive(Debug, Default)]
pub struct Alignment {
    /// Distance between target and query.
    pub edit_distance: Option<usize>,
    /// Zero-based positions in target sequence where optimal alignment paths end.
    /// * If gap after query is penalized, gap counts as part of query ([`AlignMode::NW`](crate::mode::AlignMode::NW)), otherwise not.
    pub end_locations: Option<Vec<isize>>,
    /// Zero-based positions in target sequence where optimal alignment paths start.
    /// * Correspond to [`AlignResult::end_locations`].
    /// * If gap before query is penalized, gap counts as part of query ([`AlignMode::NW`](crate::mode::AlignMode::NW)), otherwise not.
    pub start_locations: Option<Vec<isize>>,
    /// Alignment is found for first pair of start and end locations.
    ///
    /// Is a sequence of [`EditOp`]s.
    ///
    /// Alignment aligns query to target from beginning of query till end of query.
    ///
    /// If gaps are not penalized, they are not in alignment.
    pub alignment: Option<Vec<EditOp>>,
    /// Number of different characters in query and target together.
    pub alphabet_length: usize,
}

/// Transform sequences to sequences of indices.
pub fn transform_sequences(query: &str, target: &str) -> (String, Vec<usize>, Vec<usize>) {
    let mut alphabet = String::new();

    // NOTE: Original implentation assumes ASCII charset.
    let mut letter_idx: [Option<usize>; MAX_UCHAR] = [None; MAX_UCHAR];
    let mut in_alphabet: [bool; MAX_UCHAR] = [false; MAX_UCHAR];

    fn modify_seq(
        seq: impl AsRef<str>,
        alphabet: &mut String,
        letter_idx: &mut [Option<usize>; MAX_UCHAR],
        in_alphabet: &mut [bool; MAX_UCHAR],
    ) -> Vec<usize> {
        let mut seq_transformed = Vec::with_capacity(seq.as_ref().len());

        for elem in seq.as_ref().chars() {
            let elem_idx = elem as usize;
            // Idx into in_alphabet and check if not in alphabet.
            if let Some(false) = in_alphabet.get(elem_idx) {
                // SAFETY: Above ensures that elem exists in in_alphabet.
                // Both containers init to same size with MAX_UCHAR.
                unsafe {
                    // Set char to be in alphabet and assign letter idx.
                    *in_alphabet.get_unchecked_mut(elem_idx) = true;
                    *letter_idx.get_unchecked_mut(elem_idx) = Some(alphabet.len())
                }
                alphabet.push(elem)
            };
            // Safe to index in and unwrap as above adds transformed letter to letter_idx
            seq_transformed.push(letter_idx[elem_idx].unwrap());
        }
        seq_transformed
    }

    let transformed_query = modify_seq(query, &mut alphabet, &mut letter_idx, &mut in_alphabet);
    let transformed_target = modify_seq(target, &mut alphabet, &mut letter_idx, &mut in_alphabet);

    (alphabet, transformed_query, transformed_target)
}

/// Alignment information.
pub struct AlignmentData {
    pub ps: Vec<Option<Word>>,
    pub ms: Vec<Option<Word>>,
    pub scores: Vec<Option<isize>>,
    pub first_blocks: Vec<Option<usize>>,
    pub last_blocks: Vec<Option<usize>>,
}

impl AlignmentData {
    /// Init data.
    ///
    /// We build a complete table and mark first and last block for each column
    /// (because algorithm is banded so only part of each columns is used).
    /// TODO: do not build a whole table, but just enough blocks for each column.
    pub fn new(max_num_blocks: usize, target_len: usize) -> Self {
        AlignmentData {
            ps: vec![None; max_num_blocks * target_len],
            ms: vec![None; max_num_blocks * target_len],
            scores: vec![None; max_num_blocks * target_len],
            first_blocks: vec![None; target_len],
            last_blocks: vec![None; target_len],
        }
    }
}

impl Alignment {
    /// Aligns two sequences (query and target) using edit distance (levenshtein distance) returning an [`Alignment`].
    ///
    /// * @param `config`: [`AlignConfig`] configuration.
    /// * @param `query`: First sequence.
    /// * @param `target`: Second sequence.
    ///
    /// ### Example
    /// ```
    /// use rs_edlib::{align::Alignment, config::AlignConfig};
    ///
    /// let query: &str = "ACT";
    /// let target: &str = "CGT";
    /// let align_res = Alignment::run(
    ///     AlignConfig::default(),
    ///     query,
    ///     target
    /// );
    /// ```
    pub fn run(
        config: AlignConfig,
        query: impl AsRef<str>,
        target: impl AsRef<str>,
    ) -> anyhow::Result<Self> {
        let word_size = usize::try_from(WORD_SIZE)?;
        let mut alignment = Alignment::default();
        let (query, target) = (query.as_ref(), target.as_ref());

        let (alphabet, transformed_query, transformed_target) = transform_sequences(query, target);
        alignment.alphabet_length = alphabet.len();

        // Special case where one of seq is empty.
        if transformed_query.is_empty() || transformed_target.is_empty() {
            match config.mode {
                // Global alignment.
                AlignMode::NW => {
                    // Completely different.
                    alignment.edit_distance = Some(std::cmp::max(
                        transformed_query.len(),
                        transformed_target.len(),
                    ));
                    // Couldn't this potentially overflow?
                    // https://github.com/Martinsos/edlib/blob/931be2b0909985551eb17d767694a6e64e31ebfa/edlib/src/edlib.cpp#L165
                    alignment.end_locations = Some(vec![transformed_target.len() as isize - 1])
                }
                AlignMode::SHW | AlignMode::HW => {
                    alignment.edit_distance = Some(transformed_query.len());
                    alignment.end_locations = Some(vec![-1])
                }
            }
            return Ok(alignment);
        }

        // Banded block algorithm.
        // See https://en.wikipedia.org/wiki/Band_matrix
        // See supplementary material for https://academic.oup.com/bioinformatics/article/33/9/1394/2964763?login=false.

        // Initialization
        // B_max
        let max_num_blocks = ceil_div!(transformed_query.len(), word_size);
        // Number of redundant cells in last level blocks.
        let w = max_num_blocks * word_size - transformed_query.len();
        let equality_def = EqualityDefinition::new(&alphabet, Some(&config.added_equalities));
        let peq = build_peq_table(alphabet.len(), &transformed_query, &equality_def)?;

        // Main Calculation
        let mut position_nw = None;
        // let mut align_data = AlignmentData::new(max_num_blocks, target.len());
        let mut dynamic_k = false;
        let mut k = config.k.unwrap_or_else(|| {
            dynamic_k = true;
            word_size
        });

        loop {
            match config.mode {
                AlignMode::NW => {
                    alignment.calc_edit_dst_nw(
                        &peq,
                        w,
                        max_num_blocks,
                        transformed_query.len(),
                        &transformed_target,
                        k,
                        &mut position_nw,
                        None,
                        None,
                    )?;
                }
                AlignMode::SHW | AlignMode::HW => alignment.calc_edit_dst_semi_global(
                    &peq,
                    w,
                    max_num_blocks,
                    transformed_query.len(),
                    &transformed_target,
                    k,
                    &config.mode,
                )?,
            };

            k *= 2;

            if !(dynamic_k && alignment.edit_distance.is_none()) {
                break;
            }
        }

        // If there is a solution.
        if alignment.edit_distance.is_some() {
            // If NW mode, set end location explicitly.
            if config.mode == AlignMode::NW {
                alignment.end_locations = Some(vec![(transformed_target.len() - 1).try_into()?])
            }

            if matches!(config.task, AlignTask::Loc | AlignTask::Path) {
                if config.mode == AlignMode::HW {
                    let rev_transformed_target: Vec<usize> =
                        transformed_target.iter().rev().copied().collect();
                    let rev_transformed_query: Vec<usize> =
                        transformed_query.iter().rev().copied().collect();

                    let rev_peq =
                        build_peq_table(alphabet.len(), &rev_transformed_query, &equality_def)?;
                    for (i, loc) in alignment
                        .end_locations
                        .as_ref()
                        .expect("No end locations.")
                        .iter()
                        .enumerate()
                    {
                        if *loc == -1 {
                            // NOTE: Sometimes one of optimal solutions is that query starts before target, like this:
                            //                       AAGG <- target
                            //                   CCTT     <- query
                            //   It will never be only optimal solution and it does not happen often, however it is
                            //   possible and in that case end location will be -1. What should we do with that?
                            //   Should we just skip reporting such end location, although it is a solution?
                            //   If we do report it, what is the start location? -4? -1? Nothing?
                            // TODO: Figure this out. This has to do in general with how we think about start
                            //   and end locations.
                            //   Also, we have alignment later relying on this locations to limit the space of it's
                            //   search -> how can it do it right if these locations are negative or incorrect?
                            // I put 0 for now, but it does not make much sense.
                            if let Some(start_locs) = alignment.start_locations.as_mut() {
                                start_locs[i] = 0; // I put 0 for now, but it does not make much sense.
                            }
                        } else {
                            // TODO: Check index.
                            let rev_target_idx =
                                transformed_target.len() - usize::try_from(*loc)? - 1;
                            let mut rev_alignment = Alignment::default();
                            rev_alignment.calc_edit_dst_semi_global(
                                &rev_peq,
                                w,
                                max_num_blocks,
                                rev_transformed_query.len(),
                                &[rev_transformed_target[rev_target_idx]],
                                alignment.edit_distance.unwrap(),
                                &AlignMode::SHW,
                            )?;

                            if let (Some(start_locs), Some(rev_last_loc)) = (
                                alignment.start_locations.as_mut(),
                                rev_alignment.end_locations.map(|locs| locs[locs.len() - 1]),
                            ) {
                                // Taking last location as start ensures that alignment will not start with insertions
                                // if it can start with mismatches instead.
                                start_locs[i] = loc - rev_last_loc;
                            }
                        }
                    }
                } else {
                    // If mode is SHW or NW
                    if let Some(start_locs) = alignment.start_locations.as_mut() {
                        for loc in start_locs.iter_mut() {
                            *loc = 0
                        }
                    }
                }
            }
        }
        // Find alignment -> all comes down to finding alignment for NW.
        // Currently we return alignment only for first pair of locations.
        if config.task == AlignTask::Path {
            let (Some(aln_start_loc), Some(aln_end_loc)) = (
                alignment
                    .start_locations
                    .as_ref()
                    .map(|start_locs| start_locs[0]),
                alignment.end_locations.as_ref().map(|end_locs| end_locs[0]),
            ) else {
                bail!("No start locations.")
            };
            let adj_range = usize::try_from(aln_start_loc)?..usize::try_from(aln_end_loc)?;
            let aln_target = &transformed_target[adj_range];
            alignment.obtain_optimal_path(&transformed_query, aln_target, &equality_def)?;
        }
        Ok(alignment)
    }

    fn obtain_optimal_path(
        &mut self,
        query: &[usize],
        target: &[usize],
        equality_def: &EqualityDefinition,
    ) -> anyhow::Result<()> {
        // Create reversed copies.
        let rev_query = query.iter().rev();
        let rev_target = target.iter().rev();

        // Get alpha len.
        let alphabet_len = equality_def.alphabet.len();

        // Special case.
        if let Some(empty_seq_op) = if query.is_empty() {
            Some(EditOp::Delete)
        } else if target.is_empty() {
            Some(EditOp::Insert)
        } else {
            None
        } {
            self.alignment = Some(vec![empty_seq_op; target.len() + query.len()]);
            return Ok(());
        }

        let word_size = usize::try_from(WORD_SIZE)?;
        let max_num_blocks = ceil_div!(query.len(), word_size);
        let w = max_num_blocks * word_size - query.len();
        Ok(())
    }
    // pub fn as_cigar(format: CigarFormat) {}
}

mod test {
    use super::*;

    #[test]
    fn test_transform_sequences() {
        /*
        Original sequences: "ACT" and "CGT".
        Alphabet would be recognized as "ACTG". Alphabet length = 4.
        Transformed sequences: [0, 1, 2] and [1, 3, 2].
        */
        const QUERY: &str = "ACT";
        const TARGET: &str = "CGT";

        const EXP_TRANSFORMED_QUERY: [usize; 3] = [0, 1, 2];
        const EXP_TRANSFORMED_TARGET: [usize; 3] = [1, 3, 2];
        const EXP_ALPHABET: &str = "ACTG";

        let (alphabet, transformed_query, transformed_target) = transform_sequences(QUERY, TARGET);

        assert_eq!(alphabet, EXP_ALPHABET);
        assert_eq!(transformed_query, EXP_TRANSFORMED_QUERY);
        assert_eq!(transformed_target, EXP_TRANSFORMED_TARGET);
    }
}
