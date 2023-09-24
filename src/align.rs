use crate::{
    cigar::{CigarFormat, EditOp},
    config::AlignConfig,
    mode::AlignMode,
};

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
    pub end_locations: Option<Vec<usize>>,
    /// Zero-based positions in target sequence where optimal alignment paths start.
    /// * Correspond to [`AlignResult::end_locations`].
    /// * If gap before query is penalized, gap counts as part of query ([`AlignMode::NW`](crate::mode::AlignMode::NW)), otherwise not.
    pub start_locations: Option<Vec<usize>>,
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
                    alignment.end_locations = Some(vec![target.len() - 1])
                }
                AlignMode::SHW | AlignMode::HW => todo!(),
            }
        }

        Ok(alignment)
    }

    pub fn as_cigar(format: CigarFormat) {}
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
