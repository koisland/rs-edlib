use crate::{equal::EqualityPair, mode::AlignMode, task::AlignTask};

#[derive(Debug, Clone, Default)]
/// Alignment configuration.
/// * `additionalEqualitiesLength` not necessary because of Rust's storage of container metadata.
pub struct AlignConfig {
    /// Threshold number of differences between target and query. Limits search space of possible solutions.
    ///
    /// **If**:
    /// * Non-negative:
    ///     * Edit distance is not larget than `k`.
    /// * None:
    ///     * `k` will be auto-adjusted until a score is found.
    /// * Small:
    ///     * Improve speed of computation.
    /// * Smaller than edit distance:
    ///     * Edit distance will be reduced to `-1`.
    ///
    /// https://dl.acm.org/doi/abs/10.1145/316542.316550 \[1\]
    ///
    pub k: Option<usize>,
    /// Alignment method, [`AlignMode`].
    pub mode: AlignMode,
    /// Alignment task, [`AlignTask`].
    pub task: AlignTask,
    /// List of pairs of characters as an [`EqualityPair`], where each pair defines two characters as equal.
    /// * Allows extension of the lib's definition of equality.
    pub added_equalities: Vec<EqualityPair>,
}
