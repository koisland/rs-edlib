#[derive(Debug, Clone, Default, PartialEq, Eq)]
/// What should rs_edlib do?
pub enum AlignTask {
    #[default]
    /// Find edit distance and end locations.
    Distance,
    /// Find edit distance, end locations, and start locations.
    Loc,
    /// Find edit distance, end locations, and start locations and alignment path.
    Path,
}
