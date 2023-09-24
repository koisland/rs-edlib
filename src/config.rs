use std::collections::HashSet;

use crate::{cigar::EqualityPair, mode::AlignMode, task::AlignTask};

#[derive(Debug, Clone)]
/// Alignment configuration.
/// * `additionalEqualitiesLength` not necessary because of Rust's storage of container metadata.
pub struct AlignConfig {
    /// If non-negative:
    /// * Edit distance is not larget than `k`.
    ///
    /// If small:
    /// * Improve speed of computation.
    ///
    /// If smaller than edit distance:
    /// * Edit distance will be reduced to `-1`.
    ///
    /// If negative:
    /// * `k` will be auto-adjusted until a score is found.
    pub k: isize,
    /// Alignment method, [`AlignMode`].
    pub mode: AlignMode,
    /// Alignment task, [`AlignTask`].
    pub task: AlignTask,
    /// List of pairs of characters as an [`EqualityPair`], where each pair defines two characters as equal.
    /// * Allows extension of the lib's definition of equality.
    pub added_equalities: HashSet<EqualityPair>,
}

impl Default for AlignConfig {
    fn default() -> Self {
        Self {
            k: -1,
            mode: Default::default(),
            task: Default::default(),
            added_equalities: Default::default(),
        }
    }
}
