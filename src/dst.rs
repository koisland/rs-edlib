use anyhow::bail;

use crate::{
    align::{Alignment, AlignmentData, Word, WORD_1, WORD_SIZE},
    block::Block,
    mode::AlignMode,
};

// Each column is reduced in more expensive way.
pub const STRONG_REDUCE_NUM: usize = 2048;

impl Alignment {
    /// Uses Myers' bit-vector algorithm to find edit distance for one of semi-global alignment methods.
    /// * `Peq`: Query profile.
    /// * `w`: Size of padding in last block.
    /// * `maxNumBlocks`: Number of blocks needed to cover the whole query.
    /// * `queryLength`
    /// * `target`
    /// * `k`
    /// * `mode`: EDLIB_MODE_HW or EDLIB_MODE_SHW
    pub fn calc_edit_dst_semi_global(
        &mut self,
        peq: &[Word],
        w: usize,
        max_num_blocks: usize,
        query_len: usize,
        target: &[usize],
        mut k: usize,
        mode: &AlignMode,
    ) -> anyhow::Result<()> {
        let mut best_score: Option<usize> = None;
        let mut positions: Vec<isize> = vec![];

        let word_size = usize::try_from(WORD_SIZE)?;
        let int_word_size = isize::try_from(word_size)?;

        let mut first_block_idx = 0;
        let mut last_block_idx = std::cmp::min(k + 1 / word_size, max_num_blocks);

        // Initialize blocks to some score and default bit vecs.
        let mut blocks: Vec<Block> = (0..=last_block_idx)
            .map(|blk_n| Block {
                p: Word::MAX,
                m: 0,
                score: ((blk_n + 1) * word_size) as isize,
            })
            .collect();

        // If HW , gap before query is not penalized (hout == 0).
        let mut start_hout = 1;
        if *mode == AlignMode::HW {
            k = std::cmp::min(query_len, k);
            start_hout = 0;
        }
        let mut k = isize::try_from(k)?;

        for (idx, c) in target.iter().enumerate() {
            let col_idx = idx + *c * max_num_blocks;

            // calculate column
            let mut hout = start_hout;
            for block_idx in first_block_idx..=last_block_idx {
                let block = &mut blocks[block_idx];
                let peq_c = &peq[col_idx + block_idx];
                hout = block.calculate_hout_delta(*peq_c, hout)?;
                block.score += hout;
            }

            // let last_block = &mut blocks[last_block_idx - 1];
            // let peq_c = &peq[col_idx + last_block_idx - 1];
            let next_peq_c = &peq[col_idx + last_block_idx];

            if last_block_idx < max_num_blocks - 1
                && blocks[last_block_idx - 1].score - hout <= k
                && ((next_peq_c & WORD_1) != 0 || hout < 0)
            {
                let last_block = &mut blocks[last_block_idx - 1];
                last_block_idx += 1;
                last_block.p = Word::MAX;
                last_block.m = 0;
                last_block.score = last_block.score - hout
                    + int_word_size
                    + last_block.calculate_hout_delta(*next_peq_c, hout)?
            } else {
                loop {
                    let last_block = &blocks[last_block_idx - 1];
                    // CHANGE: Made > instead of >=. Why is an idx negative?
                    if !(last_block_idx > first_block_idx && last_block.score >= k + int_word_size)
                    {
                        break;
                    }
                    last_block_idx -= 1;
                }
            }

            // Every some columns, do some expensive but also more efficient block reducing.
            // This is important!
            //
            // Reduce the band by decreasing last block if possible.
            if c % STRONG_REDUCE_NUM == 0 {
                loop {
                    let last_block = &blocks[last_block_idx];
                    // CHANGE: Remove inline function.
                    // CHANGE: Made > instead of >=. Should prevent overflow but might have downstream effects.
                    // See https://github.com/Martinsos/edlib/blob/master/edlib/src/edlib.cpp#L624C3-L630
                    if !(last_block_idx > first_block_idx && last_block.all_block_cells_larger(k)) {
                        break;
                    }
                    last_block_idx -= 1
                }
            }
            // Reduce band by increasing first block if possible. Not applicable to HW.
            if *mode == AlignMode::SHW {
                while first_block_idx <= last_block_idx
                    && blocks[first_block_idx].score >= k + int_word_size
                {
                    first_block_idx += 1;
                }

                // Do strong reduction every some blcoks.
                if c % STRONG_REDUCE_NUM == 0 {
                    let first_block = &blocks[first_block_idx];
                    while first_block_idx <= last_block_idx && first_block.all_block_cells_larger(k)
                    {
                        first_block_idx += 1
                    }
                }
            }

            // If band stops to exist finish.
            if last_block_idx < first_block_idx {
                self.edit_distance = best_score;
                if best_score.is_some() {
                    self.end_locations = Some(positions)
                }

                return Ok(());
            }

            // Update best score.
            if last_block_idx == max_num_blocks - 1 {
                let col_score = blocks[last_block_idx].score;

                // Scores > k don't have correct values (so we cannot use them), but are certainly > k.
                if col_score <= k {
                    // NOTE: Score that I find in column c == score from column c - W
                    let int_best_score = best_score.map(|score| score as isize);
                    if best_score.is_none() || Some(col_score) <= int_best_score {
                        if Some(col_score) != int_best_score {
                            positions.clear();
                            best_score = Some(col_score.try_into()?);
                            // Change k so we will only look for equal or better scores than the best found so far.
                            k = col_score;
                        }
                        positions.push(isize::try_from(*c)? - isize::try_from(w)?)
                    }
                }
            }
        }

        // Obtain results for last W columns from last column.
        if last_block_idx == max_num_blocks - 1 {
            let block_scores = blocks[last_block_idx].get_cell_values();
            for block_score_idx in 0..w {
                let col_score = block_scores[block_score_idx + 1];
                let int_best_score = best_score.map(|score| score as isize);
                if col_score <= k
                    && (best_score.is_none()
                        || int_best_score.map(|score| score >= col_score) == Some(true))
                {
                    if Some(col_score) != int_best_score {
                        positions.clear();
                        k = col_score;
                        best_score = Some(col_score.try_into()?)
                    }
                    positions.push(
                        isize::try_from(target.len())? - isize::try_from(w)?
                            + isize::try_from(block_score_idx)?,
                    )
                }
            }
        }

        self.edit_distance = best_score;
        if best_score.is_some() {
            self.end_locations = Some(positions)
        }

        Ok(())
    }

    /// Use Myers' bit-vector algorithm to find edit distance for global (NW) alignment method.
    /// * `Peq`: Query profile
    /// * `w`: Size of padding in last block.
    /// * `max_num_blocks`: Number of blocks needed to cover the whole query.
    /// * `query_len`
    /// * `target`: Transformed target sequence.
    /// * `k`
    /// * `position`: 0-indexed position in target at which best score was found.
    /// * `align_data`: Data needed for alignment traceback (for reconstruction of alignment).
    ///     * Omitted original implementation `bool`, `find_alignment`.
    ///     * If `Some<AlignmentData>`, whole matrix is remembered and alignment data is returned.
    ///     * Quadratic amount of memory is consumed.
    /// * `target_stop_position`:
    ///     * If set to `None`, whole calculation is performed normally, as expected. Originally sentinel value of `-1`.
    ///     * If set to `p`, calculation is performed up to position `p` in target (inclusive) and column p is returned as the only column in align data.
    pub fn calc_edit_dst_nw(
        &mut self,
        peq: &[Word],
        w: usize,
        max_num_blocks: usize,
        query_len: usize,
        target: &[usize],
        k: usize,
        mut position: Option<&mut usize>,
        mut align_data: Option<&mut AlignmentData>,
        target_stop_position: Option<usize>,
    ) -> anyhow::Result<()> {
        let mut best_score: Option<&mut usize> = self.edit_distance.as_mut();
        let word_size = usize::try_from(WORD_SIZE)?;

        if target_stop_position.is_some() && align_data.is_some() {
            bail!("Cannot set alignment and target_stop_position at same time.");
        };

        if k < target.len().abs_diff(query_len) {
            best_score.take();
            position.take();
            return Ok(());
        }

        let mut k = std::cmp::min(k, std::cmp::max(query_len, target.len()));

        // 0-based index of first block of Ukkonen band.
        let mut first_block_idx = 0;

        /*
        Cleaned up from https://github.com/Martinsos/edlib/blob/master/edlib/src/edlib.cpp#L755C21-L755C112

        min(
            maxNumBlocks,
            ceilDiv(
                min(
                    k,
                    (k + queryLength - targetLength) / 2
                ) + 1,
                WORD_SIZE
            )
        ) - 1
        */
        // 0-based index of last block of Ukkonen band.
        // The cells below the band.
        let mut last_block_idx = std::cmp::min(
            max_num_blocks,
            std::cmp::min(k, (k + query_len - target.len() / 2) + 1) / word_size,
        ) - 1;

        // Initialize blocks to some score and default bit vecs.
        let mut blocks: Vec<Block> = (0..=last_block_idx)
            .map(|blk_n| Block {
                p: Word::MAX,
                m: 0,
                score: ((blk_n + 1) * word_size) as isize,
            })
            .collect();

        // if let Some(aln_data) = align_data.as_mut() {
        //     if target_stop_position.is_none() {
        //         // TODO: Need to resize vectors.
        //     }
        // } else {
        //     // No alignment.
        // }

        // Iterate thru columns.
        for (idx, c) in target.iter().enumerate() {
            let col_idx = idx + *c * max_num_blocks;

            // Init to hout of 1. Can range from 0, 1, and -1.
            let mut hout: isize = 1;
            // Adapted from Myers 1999.
            // Computes the output score of the block given:
            // * current level-b vertical delta vector (peq_c) given symbol c, as input
            // * horizontal delta hin (hout, init to 1)
            // Makes the resulting vertical delta the current one
            // Returning the horizontal output delta, hout.
            for block_idx in first_block_idx..=last_block_idx {
                let block = &mut blocks[block_idx];
                let peq_c = peq[col_idx + block_idx];

                hout = block.calculate_hout_delta(peq_c, hout)?;
                block.score += hout;
            }

            let mut last_block = &mut blocks[last_block_idx - 1];
            // Update k. I do it only on end of column because it would slow calculation too much otherwise. (Martinsos)
            // NOTICE: I add W when in last block because it is actually result from W cells to the left and W cells up. (Martinsos)
            /*
            https://github.com/Martinsos/edlib/blob/master/edlib/src/edlib.cpp#L795C20-L797C58
            bl->score
                + max(
                    targetLength - c - 1,
                    queryLength - (
                        (1 + lastBlock) * WORD_SIZE - 1
                    ) - 1
                )
                + (lastBlock == maxNumBlocks - 1 ? W : 0)
            */
            let cmp = last_block.score
                + isize::try_from(
                    std::cmp::max(
                        target.len() - *c - 1,
                        query_len - ((1 + last_block_idx) * word_size - 1) - 1,
                    ) + if last_block_idx == max_num_blocks - 1 {
                        w
                    } else {
                        0
                    },
                )?;
            k = std::cmp::min(k, cmp.try_into()?);

            // Adjust number of blocks accoring to Ukkonen (Martinsos)
            // Adjust last block (Martinsos)
            // If block is not beneath band, calculate next block. Only next because others are certainly beneath band. (Martinsos)
            if last_block_idx + 1 < max_num_blocks
                && (last_block_idx + 1) * word_size - 1
                    <= k - usize::try_from(last_block.score)? + 2 * word_size - 2 - target.len()
                        + c
                        + query_len
            {
                let prev_block_score = blocks[last_block_idx].score;
                // Next block.
                last_block_idx += 1;
                last_block = &mut blocks[last_block_idx];
                last_block.p = Word::MAX;
                last_block.m = 0;
                // Get eq for column.
                let peq_c = peq[col_idx + last_block_idx];

                let new_hout = last_block.calculate_hout_delta(peq_c, hout)?;
                last_block.score = prev_block_score - hout + isize::try_from(word_size)? + new_hout;

                // // Original impl never uses? Why?
                // hout = new_hout
            }

            // While block is out of band, move one block up. (Martinsos)
            // NOTE: Condition used here is more loose than the one from the article, since I simplified the max() part of it. (Martinsos)
            // I could consider adding that max part, for optimal performance. (Martinsos)
            // TODO: This might overflow.
            while last_block_idx >= first_block_idx
                && last_block.score >= (k + word_size).try_into()?
                || (last_block_idx + 1) * word_size - 1
                    > k - usize::try_from(last_block.score)? + 2 * word_size - 2 - target.len()
                        + c
                        + query_len
                        + 1
            {
                last_block_idx -= 1;
                last_block = &mut blocks[last_block_idx];
            }
            // Adjust first block (Martinsos)
            // While outside of band, advance block (Martinsos)
            while first_block_idx <= last_block_idx
                && (blocks[first_block_idx].score >= (k + word_size).try_into()?
                    || ((first_block_idx + 1) * word_size - 1
                        < (blocks[first_block_idx].score
                            - isize::try_from(k - target.len() + query_len + c)?)
                        .try_into()?))
            {
                first_block_idx += 1;
            }

            // TODO: consider if this part is useful, it does not seem to help much
            let k: isize = k.try_into()?;

            if c % STRONG_REDUCE_NUM == 0 {
                // Every some columns do more expensive but more efficient reduction
                'outer: while last_block_idx >= first_block_idx {
                    // If all cells outside of band remove block.
                    let scores = blocks[last_block_idx].get_cell_values();
                    let num_cells = if last_block_idx == max_num_blocks - 1 {
                        word_size - w
                    } else {
                        word_size
                    };
                    let mut r = last_block_idx * word_size + num_cells - 1;

                    for score in scores.iter().take(word_size).skip(word_size - num_cells) {
                        // TODO: Does not work if do not put +1! Why???
                        if (*score <= k)
                            && (isize::try_from(r)?
                                <= k - score - isize::try_from(target.len() + c + query_len)? + 1)
                        {
                            break 'outer;
                        }
                        r -= 1;
                    }
                    // TODO: This might overflow.
                    last_block_idx -= 1;
                }

                'outer: while first_block_idx <= last_block_idx {
                    let scores = blocks[first_block_idx].get_cell_values();
                    let num_cells = if first_block_idx == max_num_blocks - 1 {
                        word_size - w
                    } else {
                        word_size
                    };
                    let mut r = first_block_idx * word_size + num_cells - 1;

                    for score in scores.iter().take(word_size).skip(word_size - num_cells) {
                        // TODO: Does not work if do not put +1! Why???
                        if (*score <= k)
                            && (isize::try_from(r)?
                                <= k - score - isize::try_from(target.len() + c + query_len)? + 1)
                        {
                            break 'outer;
                        }
                        r -= 1;
                    }
                    first_block_idx += 1;
                }
            }

            // If band stops to exist finish.
            if last_block_idx < first_block_idx {
                best_score.take();
                position.take();
                return Ok(());
            }

            // Save column so it can be used for reconstruction
            if let Some(align_data) = align_data.as_mut() {
                if *c < target.len() {
                    for (b, block) in blocks[first_block_idx..=last_block_idx].iter().enumerate() {
                        let idx = max_num_blocks + c + b;
                        align_data.ps[idx] = Some(block.p);
                        align_data.ms[idx] = Some(block.m);
                        align_data.scores[idx] = Some(block.score);
                    }
                    align_data.first_blocks[*c] = Some(first_block_idx);
                    align_data.last_blocks[*c] = Some(last_block_idx);
                }

                // If this is stop column, save it and finish. (M)
                if Some(*c) == target_stop_position {
                    for (b, block) in blocks[first_block_idx..=last_block_idx].iter().enumerate() {
                        align_data.ps[b] = Some(block.p);
                        align_data.ms[b] = Some(block.m);
                        align_data.scores[b] = Some(block.score);
                    }
                    align_data.first_blocks[0] = Some(first_block_idx);
                    align_data.last_blocks[0] = Some(last_block_idx);
                    best_score.take();

                    if let (Some(position), Some(stop_position)) =
                        (position.as_mut(), target_stop_position)
                    {
                        **position = stop_position;
                    }
                    return Ok(());
                }
            }
        }

        // If last block of last column was calculated.
        if last_block_idx == max_num_blocks - 1 {
            let last_best_score = &blocks[last_block_idx].get_cell_values()[w];
            if last_best_score <= &isize::try_from(k)? {
                if let (Some(score), Some(pos)) = (best_score, position) {
                    *score = usize::try_from(*last_best_score)?;
                    *pos = target.len() - 1;
                }
                return Ok(());
            }
        }

        best_score.take();
        position.take();
        Ok(())
    }
}
