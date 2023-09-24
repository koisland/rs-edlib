#[derive(Debug, Clone, Default)]
/// How to treat gaps before and after query?
pub enum AlignMode {
    #[default]
    /// Global method.
    /// * Useful when you want to find out how similar is first sequence (target) to second sequence (query).
    NW,
    /// Prefix method.
    /// * **Doesn't penalize gaps at end of query**. Deleting elements at end of second sequence is free.
    /// * Useful when you want to find out **how well first sequence fits at the beginning of the second sequence**.
    ///
    /// ### Example
    /// `AACT` and `AACTGGC`
    /// * Edit distance would be 0.
    /// * Removing `GGC` from the end of second sequence is "free" and does not count into total edit distance.
    SHW,
    /// Infix method.
    /// * **Gaps at query end and start are not penalized.** Deleting elements at start and end of second sequence is free.
    /// ### Example
    /// `ACT` and `CGACTGAC`
    /// * Edit distance would be 0
    /// * Removing `CG` from the start and `GAC` from the end of second sequence is "free" and does not count into total edit distance.
    /// * In bioinformatics, this method is appropriate for aligning read to a sequence.
    HW,
}
