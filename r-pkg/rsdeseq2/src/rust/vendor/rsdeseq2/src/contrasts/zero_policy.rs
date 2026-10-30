//! Shared row overrides, independent of numerical storage representation.
use crate::errors::{DeseqError, invalid_dimensions};

/// Validated borrowed contrast/global-zero flags.
///
/// Globally zero genes are excluded from contrast overrides. Wald overrides
/// effect/statistic/p-value with 0/0/1 and retains the SE. LRT overrides only
/// the displayed effect, retaining the stored likelihood-ratio inference.
#[derive(Clone, Copy, Debug)]
pub struct ContrastZeroPolicy<'a> {
    contrast: &'a [bool],
    global: &'a [bool],
}

impl<'a> ContrastZeroPolicy<'a> {
    /// Validate aligned flags before any result field can be mutated.
    pub fn new(contrast: &'a [bool], global: &'a [bool]) -> Result<Self, DeseqError> {
        if global.len() != contrast.len() {
            return Err(invalid_dimensions(
                "allZero rows",
                contrast.len(),
                global.len(),
            ));
        }
        Ok(Self { contrast, global })
    }

    fn validate_rows(&self, rows: usize) -> Result<(), DeseqError> {
        if rows != self.contrast.len() {
            return Err(invalid_dimensions(
                "contrastAllZero rows",
                rows,
                self.contrast.len(),
            ));
        }
        Ok(())
    }

    /// Apply the shared Wald override to aligned numerical columns.
    /// `Option<f64>` and plain f64 storage follow the same policy; conversion
    /// only represents the assigned constants in the destination value type.
    pub fn apply_wald<T: From<f64>>(
        &self,
        effect: &mut [T],
        statistic: &mut [T],
        pvalue: &mut [T],
    ) -> Result<(), DeseqError> {
        self.validate_rows(effect.len())?;
        self.validate_rows(statistic.len())?;
        self.validate_rows(pvalue.len())?;
        for gene in 0..effect.len() {
            if self.contrast[gene] && !self.global[gene] {
                effect[gene] = T::from(0.0);
                statistic[gene] = T::from(0.0);
                pvalue[gene] = T::from(1.0);
            }
        }
        Ok(())
    }

    /// Apply the LRT display-effect override without touching SEs or inference.
    /// Accepts either contiguous effects or borrowed fields of result rows.
    pub fn apply_effects<'b, T: From<f64> + 'b>(
        &self,
        effects: impl ExactSizeIterator<Item = &'b mut T>,
    ) -> Result<(), DeseqError> {
        self.validate_rows(effects.len())?;
        for (gene, effect) in effects.enumerate() {
            if self.contrast[gene] && !self.global[gene] {
                *effect = T::from(0.0);
            }
        }
        Ok(())
    }
}
