use crate::{
    align::{Word, WORD_SIZE},
    ceil_div,
    equal::EqualityDefinition,
};

/// Build Peq (query profile) table for given transformed query and alphabet.
/// * Peq is table of dimensions `alphabetLength+1 x maxNumBlocks`.
/// * Bit `i` of `Peq[s * maxNumBlocks + b]` is `1` if `i`-th symbol from block `b` of query equals symbol `s`, otherwise it is `0`.
pub fn build_peq_table(
    alphabet_length: usize,
    query: &[usize],
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
    let max_num_blocks = ceil_div!(query.len(), word_size);
    // Table of dimensions:
    // * (alphabet_length + 1) x max_num_blocks.
    // Last symbol is a wildcard?
    let mut peq_table: Vec<Word> = vec![0; (alphabet_length + 1) * max_num_blocks];

    for symbol in 0..=alphabet_length {
        // Get character for symbol idx.
        // If None, on last wildcard symbol. Set wildcard.
        if let Some(symbol_char) = equality_def.symbol(symbol) {
            for block in 0..max_num_blocks {
                let idx = symbol * max_num_blocks + block;
                peq_table[idx] = 0;

                let r = (block + 1) * word_size;

                for r in (0..r).rev().filter(|r| *r >= block * word_size) {
                    peq_table[idx] <<= 1;

                    // This is C++ pointer indexing behavior (query[r])? wtf just defaults to 0?
                    let r_char_idx = query.get(r).cloned().unwrap_or(0);
                    // Cast query u8 byte at r into char.
                    let r_char = equality_def.symbol(r_char_idx).unwrap();

                    // If position is greater than query len, treat as wildcard and pad with W wildcard symbols?
                    // - OR -
                    // Set to 1 if i-th symbol from block b of query equals symbol.
                    if r >= query.len() || equality_def.are_equal(r_char, symbol_char)? {
                        peq_table[idx] += 1
                    }
                }
            }
        } else {
            for block in 0..max_num_blocks {
                let idx = symbol * max_num_blocks + block;
                // All 1s. i.e. 111..11
                peq_table[idx] = Word::MAX
            }
        };
    }

    Ok(peq_table)
}

mod test {
    use super::*;
    use crate::align::transform_sequences;

    #[test]
    fn test_build_peq_table() {
        let query = str::repeat("AGGATACA", 10);
        let (alphabet, transformed_query, _) = transform_sequences(&query, &query);

        let eq_def = EqualityDefinition::new(&alphabet, None);
        let table = build_peq_table(alphabet.len(), &transformed_query, &eq_def).unwrap();

        assert_eq!(
            table,
            [
                3002117172780181929,
                434041037028460038,
                1157442765409226768,
                4629771061636907072,
                18446744073709551615
            ]
        )
    }
}
