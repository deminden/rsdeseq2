//! Resource policy independent of the likelihood and frontend.
use crate::Error;
use std::sync::{
    Arc,
    atomic::{AtomicBool, Ordering},
};

/// Shared cooperative cancellation. Cancelling never changes accepted arithmetic.
#[derive(Clone, Debug, Default)]
pub struct PosteriorCancellation(Arc<AtomicBool>);
impl PosteriorCancellation {
    /// Request cancellation; running fits stop at their next check.
    pub fn cancel(&self) {
        self.0.store(true, Ordering::Relaxed);
    }
    /// Whether cancellation was requested.
    pub fn is_cancelled(&self) -> bool {
        self.0.load(Ordering::Relaxed)
    }
}

/// Batch-wide grid budget and a cooperative cancellation source.
#[derive(Clone, Debug)]
pub struct PosteriorWorkControl {
    /// Maximum grid likelihood evaluations across the batch; None is unlimited.
    pub max_grid_evaluations: Option<usize>,
    /// Cancellation source shared with the caller.
    pub cancellation: PosteriorCancellation,
    interrupt: Option<(std::thread::ThreadId, fn() -> bool)>,
}
impl Default for PosteriorWorkControl {
    fn default() -> Self {
        Self {
            max_grid_evaluations: Some(10_000_000),
            cancellation: Default::default(),
            interrupt: None,
        }
    }
}
impl PosteriorWorkControl {
    /// Attach a polling function confined to the calling thread. It must return
    /// normally, including on interruption; it must never unwind across an FFI.
    pub fn with_interrupt_check(mut self, check: fn() -> bool) -> Self {
        self.interrupt = Some((std::thread::current().id(), check));
        self
    }
    #[doc(hidden)]
    pub fn check(&self) -> Result<(), Error> {
        if self.cancellation.is_cancelled() {
            return Err(Error::PosteriorCancelled);
        }
        if let Some((thread, check)) = self.interrupt
            && thread == std::thread::current().id()
            && check()
        {
            self.cancellation.cancel();
            return Err(Error::PosteriorCancelled);
        }
        Ok(())
    }
    pub(super) fn validate_estimate(&self, evaluations: usize) -> Result<(), Error> {
        self.check()?;
        if let Some(limit) = self.max_grid_evaluations
            && evaluations > limit
        {
            return Err(Error::PosteriorGridBudgetExceeded { evaluations, limit });
        }
        Ok(())
    }
}
