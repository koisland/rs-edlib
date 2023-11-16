/// Describes CIGAR format.
/// * See http://samtools.github.io/hts-specs/SAMv1.pdf
/// * See http://drive5.com/usearch/manual/cigar.html
pub enum CigarFormat {
    /// Match: 'M', Insertion: 'I', Deletion: 'D', Mismatch: 'M'.
    Standard,
    /// Match: '=', Insertion: 'I', Deletion: 'D', Mismatch: 'X'.
    Extended,
}


#[derive(Debug, Clone, Copy)]
pub enum EditOp {
    /// Match
    Match,
    /// Insertion to target = deletion from query.
    ///
    /// ### Example:
    /// * Target: `ATCG`
    /// * Query:  `A-CG`
    Insert,
    /// Deletion from target = insertion to query.
    ///
    /// ### Example:
    /// * Target: `A-CG`
    /// * Query:  `ATCG`
    Delete,
    /// Mismatch
    Mismatch,
}
