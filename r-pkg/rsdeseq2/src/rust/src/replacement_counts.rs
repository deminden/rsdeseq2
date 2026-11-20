//! Short-lived owned output; never contains a pointer into R storage.
use super::{copy_counts_to_r, guard_value, output_slice, run_backend};
use rsdeseq2::prelude::CountMatrix;
use std::ptr;

/// Opaque across the ABI. C may query/copy it, then must consume its owner slot.
pub struct ReplacementCounts {
    counts: CountMatrix,
}
impl ReplacementCounts {
    pub(super) fn new(counts: CountMatrix) -> Self {
        Self { counts }
    }
}

/// # Safety
/// `owner` is a live owner returned by fit. Both output slots are distinct,
/// aligned writable usize values, and no pointer is retained.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rsdeseq2_r_replacement_shape(
    owner: *const ReplacementCounts,
    genes: *mut usize,
    samples: *mut usize,
) -> i32 {
    run_backend(|| {
        // SAFETY: C passes null (handled here) or a pointer to a live, aligned value it owns for this call.
        let owner = unsafe { owner.as_ref() }.ok_or("replacement owner is null")?;
        // SAFETY: the C caller passes a writable buffer of this length, uniquely borrowed for this call (see the export's safety contract).
        let genes_out = unsafe { output_slice(genes, 1, "replacement genes")? };
        // SAFETY: the C caller passes a writable buffer of this length, uniquely borrowed for this call (see the export's safety contract).
        let samples_out = unsafe { output_slice(samples, 1, "replacement samples")? };
        genes_out[0].write(owner.counts.n_genes());
        samples_out[0].write(owner.counts.n_samples());
        Ok(())
    })
}

/// # Safety
/// `owner` is live, unaliased with the output, and remains owned by the caller.
/// C validates output allocation for genes*samples doubles before this call.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rsdeseq2_r_replacement_copy(
    owner: *const ReplacementCounts,
    genes: usize,
    samples: usize,
    output: *mut f64,
) -> i32 {
    run_backend(|| {
        // SAFETY: C passes null (handled here) or a pointer to a live, aligned value it owns for this call.
        let owner = unsafe { owner.as_ref() }.ok_or("replacement owner is null")?;
        if (genes, samples) != (owner.counts.n_genes(), owner.counts.n_samples()) {
            return Err("replacement output dimensions do not match owner".into());
        }
        let size = genes
            .checked_mul(samples)
            .ok_or("replacement shape overflow")?;
        // SAFETY: the C caller passes a writable buffer of this length, uniquely borrowed for this call (see the export's safety contract).
        unsafe {
            copy_counts_to_r(
                &owner.counts,
                output_slice(output, size, "replacement counts")?,
            )?
        };
        Ok(())
    })
}

/// # Safety
/// A non-null slot is writable and contains NULL or the unique live owner
/// returned by fit. This clears the slot before releasing its allocation.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rsdeseq2_r_replacement_drop(slot: *mut *mut ReplacementCounts) {
    guard_value((), || {
        if slot.is_null() {
            return;
        }
        // SAFETY: the caller passes the owner slot it received from this backend.
        let owner = unsafe { ptr::replace(&mut *slot, ptr::null_mut()) };
        if !owner.is_null() {
            // SAFETY: a non-null owner was produced by Box::into_raw here.
            drop(unsafe { Box::from_raw(owner) });
        }
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::mem::MaybeUninit;

    #[test]
    fn owned_counts_validate_copy_shape_and_release_the_slot() {
        let counts = CountMatrix::from_row_major_u32(2, 3, vec![1, 2, 3, 4, 5, 6]).unwrap();
        let mut owner = Box::into_raw(Box::new(ReplacementCounts::new(counts)));
        let (mut genes, mut samples) = (0, 0);
        let mut untouched = [-1.0; 6];
        let mut output = MaybeUninit::<[f64; 6]>::uninit();
        // SAFETY: test-owned buffers and pointers outlive every call in this block.
        unsafe {
            assert_eq!(
                rsdeseq2_r_replacement_shape(owner, &mut genes, &mut samples),
                0
            );
            assert_eq!((genes, samples), (2, 3));
            assert_ne!(
                rsdeseq2_r_replacement_copy(owner, 3, 2, untouched.as_mut_ptr()),
                0
            );
            assert_eq!(untouched, [-1.0; 6]);
            assert_ne!(rsdeseq2_r_replacement_copy(owner, 2, 3, ptr::null_mut()), 0);
            assert_eq!(
                rsdeseq2_r_replacement_copy(owner, 2, 3, output.as_mut_ptr().cast()),
                0
            );
            // A successful copy writes the complete 2×3 destination.
            assert_eq!(output.assume_init(), [1.0, 4.0, 2.0, 5.0, 3.0, 6.0]);
            rsdeseq2_r_replacement_drop(&mut owner);
            assert!(owner.is_null());
            rsdeseq2_r_replacement_drop(&mut owner);
            rsdeseq2_r_replacement_drop(ptr::null_mut());
            assert_ne!(
                rsdeseq2_r_replacement_shape(owner, &mut genes, &mut samples),
                0
            );
        }
    }
}
