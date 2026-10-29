//! Per-builder worker pools, shared by stage/refit clones.
use crate::errors::DeseqError;
use std::sync::{Arc, OnceLock};

#[derive(Clone, Debug, Default)]
pub(super) struct ExecutionContext {
    pool: Arc<OnceLock<Result<Arc<rayon::ThreadPool>, String>>>,
}

impl ExecutionContext {
    pub(super) fn pool_for_entry(
        &self,
        requested: Option<usize>,
    ) -> Result<Option<Arc<rayon::ThreadPool>>, DeseqError> {
        let Some(workers) = requested else {
            return Ok(None);
        };
        if workers == 0 || workers > rayon::max_num_threads() {
            return Err(DeseqError::InvalidOptions {
                reason: "threads must be positive and within Rayon's supported range".into(),
            });
        }
        let pool = self.pool.get_or_init(|| {
            rayon::ThreadPoolBuilder::new()
                .num_threads(workers)
                .build()
                .map(Arc::new)
                .map_err(|error| error.to_string())
        });
        let pool = pool.as_ref().map_err(|error| DeseqError::InvalidOptions {
            reason: format!("cannot create worker pool: {error}"),
        })?;
        if pool.current_thread_index().is_some() {
            Ok(None)
        } else {
            Ok(Some(Arc::clone(pool)))
        }
    }
}
