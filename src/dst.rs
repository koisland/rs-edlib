use anyhow::{bail, Context};

use crate::{
    align::{Alignment, AlignmentData, Word, WORD_SIZE},
    block::Block,
    mode::AlignMode,
};

impl Alignment {
    pub fn calc_edit_dst_semi_global(
        &mut self,
        peq: &[Word],
        w: usize,
        max_num_blocks: usize,
        query_len: usize,
        target: &[usize],
        k: usize,
        mode: &AlignMode,
    ) {
        let best_score = self.edit_distance.as_mut();
        let positions = &self.end_locations;
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
        align_data: Option<&mut AlignmentData>,
        target_stop_position: Option<usize>,
    ) -> anyhow::Result<()> {
        let mut best_score: Option<&mut usize> = self.edit_distance.as_mut();
        let word_size = usize::try_from(WORD_SIZE)?;

        if target_stop_position.is_some() && align_data.is_some() {
            bail!("Cannot set alignment and target_stop_position at same time.");
        };

        // Each column is reduced in more expensive way.
        const STRONG_REDUCE_NUM: usize = 2048;

        if k < target.len().abs_diff(query_len) {
            best_score.take();
            position.take();
            return Ok(());
        }

        let upper_k_bound = std::cmp::min(k, std::cmp::max(query_len, target.len()));

        // 0-based index of first block of Ukkonen band.
        let first_block_idx = 0;

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
        let last_block_idx = std::cmp::min(
            max_num_blocks,
            std::cmp::min(k, (k + query_len - target.len() / 2) + 1) / usize::try_from(WORD_SIZE)?,
        ) - 1;

        // Initialize blocks to some score and default bit vecs.
        let mut blocks: Vec<Block> = (0..=last_block_idx)
            .map(|blk_n| Block {
                p: Word::MAX,
                m: 0,
                score: ((blk_n + 1) * word_size) as isize,
            })
            .collect();

        if let Some(aln_data) = align_data {
            if target_stop_position.is_none() {
                // TODO: Need to resize vectors.
            }
        } else {
            // No alignment.
        }

        // Iterate thru columns.
        for (idx, c) in target.iter().enumerate() {
            let peq_c = peq[idx + *c * max_num_blocks];

            // Init to hout of 1. Can range from 0, 1, and -1.
            let mut hout: isize = 1;
            // Adapted from Myers 1999.
            // Computes the output score of the block given:
            // * current level-b vertical delta vector (peq_c) given symbol c, as input
            // * horizontal delta hin (hout, init to 1)
            // Makes the resulting vertical delta the current one
            // Returning the horizontal output delta, hout.
            for block in blocks
                .get_mut(first_block_idx..=last_block_idx)
                .with_context(|| {
                    format!("Invalid first {first_block_idx} or last {last_block_idx} block for blocks.")
                })?
                .iter_mut()
            {
                hout = block.calculate_hout_delta(peq_c, hout)?;
                block.score += hout;
            }

            let Some(last_block) = blocks.get_mut(last_block_idx) else {
                bail!("Missing last block.")
            };
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
            let cmp = if last_block_idx == max_num_blocks - 1 {
                w
            } else {
                0
            };
            let cmp = std::cmp::max(
                target.len() - *c - 1,
                query_len - ((1 + last_block_idx) * word_size - 1) - 1,
            ) + cmp;
            let cmp = last_block.score + isize::try_from(cmp)?;
            let new_k = std::cmp::min(k.try_into()?, cmp);
        }
        Ok(())
    }
}
