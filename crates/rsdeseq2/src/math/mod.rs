//! Small numerical helpers used by the core algorithms.

pub mod distributions;
pub mod median;
pub mod optim;
pub mod qr;
pub mod special;

pub use median::{median, median_finite};
pub use special::trigamma;
