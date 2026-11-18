//! Binary64 stand-in for the x87 kernels on hosts without them.
//!
//! [`crate::NumericProfile::RLinuxX86`] cannot be selected on these hosts, so
//! these functions are never reached; they exist so the profile branches
//! type-check everywhere.

#[derive(Clone, Copy)]
pub(super) struct Extended(f64);
impl Extended {
    pub(super) fn divide(&mut self, value: f64) {
        self.0 /= value;
    }
    pub(super) fn value(self) -> f64 {
        self.0
    }
}
pub(super) fn sum(values: impl Iterator<Item = f64>) -> Extended {
    Extended(values.sum())
}
pub(super) fn mean(values: &[f64]) -> f64 {
    values.iter().sum::<f64>() / values.len() as f64
}
pub(super) fn norm(values: &[f64]) -> f64 {
    values.iter().map(|value| value * value).sum::<f64>().sqrt()
}
