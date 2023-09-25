use crate::{
    align::{Alignment, Word},
    mode::AlignMode,
};

impl Alignment {
    pub fn myers_calc_edit_dst_nw(
        &mut self,
        peq: &[Word],
        w: usize,
        max_num_blocks: usize,
        query_len: usize,
        target: &str,
        k: isize,
        mode: &AlignMode,
    ) {
    }
}
