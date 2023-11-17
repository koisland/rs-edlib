#![warn(missing_docs)]

//! Rust port of edlib.

pub mod align;
pub mod block;
pub mod cigar;
pub mod config;
pub mod dst;
pub mod equal;
pub mod mode;
pub mod peq;
pub mod task;

#[macro_export]
/// Ceiling division.
///
/// ```ignore
/// assert_eq!(3, ceil_div!(5, 2))
/// ```
macro_rules! ceil_div {
    ($x:expr, $y:expr) => {
        ($x + $y - 1) / $y
    };
}

pub use align::Alignment;
pub use config::AlignConfig;
