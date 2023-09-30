use crate::align::{Word, HIGH_BIT_MASK, WORD_SIZE};

#[derive(Debug, Clone, Default)]
pub struct Block {
    /// Bit vector Pvin
    pub p: Word,
    /// Bit vector Mvin
    pub m: Word,
    /// Score of last cell in block.
    pub score: isize,
}

impl Block {
    /// Corresponds to `Advance_Block` function from Myers 1999.
    /// * Keep in mind. `v` for vertical and `h` for horizontal.
    /// * Calculates the horizontal output delta of one word(block), which is part of a column.
    /// * Highest bit of word (one most to the left) is most bottom cell of block from column.
    /// * `Pv[i]` and `Mv[i]` define `vin` of `cell[i]: vin = cell[i] - cell[i-1]`.
    ///
    /// # Arguments
    /// * `eq` -  Bitset, `Eq[i] == 1` if match, `0` if mismatch.
    /// * `hin` - Will be `+1`, `0` or `-1`.
    ///
    /// # Returns
    /// * `hout` -s Will be `+1`, `0` or `-1`.
    pub fn calculate_hout_delta(&mut self, mut eq: Word, hin: isize) -> anyhow::Result<isize> {
        // * Pv  Bitset, Pv[i] == 1 if vin is +1, otherwise Pv[i] == 0.
        // * Mv  Bitset, Mv[i] == 1 if vin is -1, otherwise Mv[i] == 0.

        let hin_is_neg: Word = hin.is_negative().into();
        let xv: Word = eq | self.m;

        // This is instruction below written using 'if': if (hin < 0) Eq |= (Word)1;
        eq |= hin_is_neg;
        let xh: Word = (((eq & self.p) + self.p) ^ self.p) | eq;

        let mut ph = self.m | !(xh | self.p);
        let mut mh = self.p & xh;

        let mut hout: isize = 0;
        // This is instruction below written using 'if': if (Ph & HIGH_BIT_MASK) hout = 1;
        hout = ((ph & HIGH_BIT_MASK) >> (WORD_SIZE - 1)).try_into()?;
        // This is instruction below written using 'if': if (Mh & HIGH_BIT_MASK) hout = -1;
        hout -= TryInto::<isize>::try_into((mh & HIGH_BIT_MASK) >> (WORD_SIZE - 1))?;

        ph <<= 1;
        mh <<= 1;

        mh |= hin_is_neg;
        ph |= Word::try_from((hin + 1) >> 1)?;

        // * PvOut  Bitset, PvOut[i] == 1 if vout is +1, otherwise PvOut[i] == 0.
        // * MvOut  Bitset, MvOut[i] == 1 if vout is -1, otherwise MvOut[i] == 0.
        self.p = mh | !(xv | ph);
        self.m = ph & xv;
        Ok(hout)
    }

    /// Get score values for `Block`.
    pub fn get_cell_values(&self) -> Vec<isize> {
        let mut scores = vec![];
        let mut score = self.score;
        // Mask shifted to go through block p and m by each bit modifying score.
        let mut mask = HIGH_BIT_MASK;

        // Iterate through word calculating score.
        for _ in 0..WORD_SIZE - 1 {
            scores.push(score);

            // Look at the pvin and mvin of the bottom cell of block from column. i.e. 100..00.
            if self.p & HIGH_BIT_MASK != 0 {
                score -= 1
            }
            if self.m & HIGH_BIT_MASK != 0 {
                score += 1
            }
            // Then keep moving mask up the column.
            // 100..00 -> 010..00
            mask >>= 1;
        }
        if let Some(old_score) = scores.get_mut(WORD_SIZE as usize - 1) {
            *old_score = score
        } else {
            scores.push(score)
        }

        scores
    }
}
