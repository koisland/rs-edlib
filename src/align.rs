use crate::{
    cigar::{CigarFormat, EditOp},
    config::AlignConfig,
    equal::EqualityDefinition,
    mode::AlignMode,
};
use anyhow::Error;

pub type Word = u64;
/// Size of word in bits.
const WORD_SIZE: u32 = Word::BITS;
const WORD_1: Word = 1;
// 100..00? Large num bits?
const HIGH_BIT_MASK: Word = WORD_1 << (WORD_SIZE - 1);
/// Max chars we expect. Note is ASCII only.
const MAX_UCHAR: usize = 256;

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

pub struct AlignmentData {
    pub ps: Vec<Word>,
    pub ms: Vec<Word>,
    pub scores: Vec<usize>,
    pub first_blocks: Vec<usize>,
    pub last_blocks: Vec<usize>,
}

impl AlignmentData {
    pub fn new(max_num_blocks: usize, target_len: usize) -> Self {
        AlignmentData {
            ps: Vec::with_capacity(max_num_blocks * target_len),
            ms: Vec::with_capacity(max_num_blocks * target_len),
            scores: Vec::with_capacity(max_num_blocks * target_len),
            first_blocks: Vec::with_capacity(target_len),
            last_blocks: Vec::with_capacity(target_len),
        }
    }
}

pub struct Block {
    /// Pvin
    pub p: Word,
    /// Mvin
    pub m: Word,
    /// Score of last cell in block.
    pub score: usize,
}

impl Alignment {
    /// Aligns two sequences (query and target) using edit distance (levenshtein distance) returning an [`Alignment`].
    ///
    /// * @param `config`: [`AlignConfig`] configuration.
    /// * @param `query`: First sequence.
    /// * @param `target`: Second sequence.
    ///
    /// ### Exampleq
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
        let mut alignment = Alignment::default();
        let (query, target) = (query.as_ref(), target.as_ref());

        let (alphabet, transformed_query, transformed_target) = transform_sequences(query, target);
        alignment.alphabet_length = alphabet.len();

        // Special case where one of seq is empty.
        if query.is_empty() || target.is_empty() {
            match config.mode {
                // Global alignment.
                AlignMode::NW => {
                    // Completely different.
                    alignment.edit_distance = Some(std::cmp::max(query.len(), target.len()));
                    // Couldn't this potentially overflow?
                    // https://github.com/Martinsos/edlib/blob/931be2b0909985551eb17d767694a6e64e31ebfa/edlib/src/edlib.cpp#L165
                    alignment.end_locations = Some(vec![target.len() as isize - 1])
                }
                AlignMode::SHW | AlignMode::HW => {
                    alignment.edit_distance = Some(query.len());
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
        let max_num_blocks = query.len() / usize::try_from(WORD_SIZE)?;
        // Number of redundant cells in last level blocks.
        let w = max_num_blocks * usize::try_from(WORD_SIZE)? - query.len();
        let equality_def = EqualityDefinition::new(&alphabet, Some(&config.added_equalities));
        let peq = build_peq_table(alphabet.len(), query, &equality_def)?;

        // Main Calculation
        let mut align_data = AlignmentData::new(max_num_blocks, target.len());
        let mut k = config.k;
        let mut dynamic_k = false;
        if k < 0 {
            dynamic_k = true;
            k = isize::try_from(WORD_SIZE)?;
        }
        loop {
            match config.mode {
                AlignMode::NW => alignment.myers_calc_edit_dst_nw(
                    &peq,
                    w,
                    max_num_blocks,
                    query.len(),
                    target,
                    k,
                    &config.mode,
                ),
                AlignMode::SHW | AlignMode::HW => todo!(),
            }
            k *= 2;

            if !(dynamic_k && alignment.edit_distance.is_none()) {
                break;
            }
        }

        if let Some(dst) = alignment.edit_distance {}

        Ok(alignment)
    }

    pub fn as_cigar(format: CigarFormat) {}
}

/// Build Peq (query profile) table for given query and alphabet.
/// * Peq is table of dimensions `alphabetLength+1 x maxNumBlocks`.
/// * Bit `i` of `Peq[s * maxNumBlocks + b]` is `1` if `i`-th symbol from block `b` of query equals symbol `s`, otherwise it is `0`.
pub fn build_peq_table(
    alphabet_length: usize,
    query: &str,
    equality_def: &EqualityDefinition,
) -> anyhow::Result<Vec<Word>> {
    /*
        query_len = 2000
        max_num_blocks = 31.25 -> 31
        alphabet_len = 4
        alphabet = "ATGC" (0123)

              A T C G
          b1  0 0 0 0 0
          b2  0 0 0 0 0
              ... x 1996
          b4  0 0 0 0 0
          b5  0 0 0 0 0

        symbol = 0 (A)
        block = 0 (b1)
        r = (0 + 1) * (64 - 1) -> 63
        idx = 0 * 31 + 0 -> 0
    */
    let word_size = usize::try_from(WORD_SIZE)?;
    let max_num_blocks = query.len() / word_size;
    // Table of dimension alphabet length + 1 x max_num_blocks.
    // Last symbol is a wildcard?
    let mut peq_table: Vec<Word> = vec![0; (alphabet_length + 1) * max_num_blocks];

    for symbol in 0..=alphabet_length {
        let Some(symbol_char) = equality_def.symbol(symbol) else {
            continue;
        };
        for block in 0..max_num_blocks {
            let idx = symbol * max_num_blocks + block;

            if symbol < alphabet_length {
                peq_table[idx] = 0;
                let r = (block + 1) * word_size - 1;

                for r in (0..r).rev() {
                    if r < block * word_size {
                        break;
                    }
                    // Cast query u8 byte at r into char.
                    // NOTE: Can panic.
                    let Some(r_char) = query.as_bytes().get(r).and_then(|c| char::from_u32(*c as u32)) else {
                        eprintln!("Cannot get query elem at {r}.");
                        continue;
                    };

                    // Default should be 0.
                    // Bitshift to set to 0? Why not just set to 0?
                    peq_table[idx] = 0;

                    // If position is greater than query len, treat as wildcard and pad with W wildcard symbols?
                    // - OR -
                    // Set to 1 if i-th symbol from block b of query equals symbol.
                    if r >= query.len() || equality_def.are_equal(r_char, symbol_char)? {
                        peq_table[idx] += 1
                    }
                }
            } else {
                // Last symbol. Wildcard.
                peq_table[idx] = 1
            }
        }
    }

    Ok(peq_table)
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

    #[test]
    fn test_build_peq_table() {
        const ALPHABET: &str = "ATGC";
        const QUERY: &str = "AGGATACA";
        let eq_def = EqualityDefinition::new(ALPHABET, None);
        let table = build_peq_table(ALPHABET.len(), &str::repeat(QUERY, 200), &eq_def).unwrap();

        dbg!(table);
    }
}
