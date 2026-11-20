//! Synchronous count/row preparation and offset composition; no R numerical APIs.
use super::{input_slice, output_slice, run_backend};
use rsdeseq2::glm::NormalShrinkageRows;
use rsdeseq2::matrix::{CountRows, CountView, MatrixOrder, MatrixView};
use rsdeseq2::normalization::visit_offset_factors_column_major;
use rsdeseq2::results::RowFlags;

unsafe fn counts<'a>(
    genes: usize,
    samples: usize,
    values: *const f64,
) -> Result<MatrixView<'a, f64>, String> {
    let length = genes.checked_mul(samples).ok_or("count shape overflow")?;
    MatrixView::new(
        genes,
        samples,
        // SAFETY: the C caller passes a readable buffer of these dimensions that stays valid for this call (see the export's safety contract).
        unsafe { input_slice(values, length, "counts")? },
        MatrixOrder::ColumnMajor,
    )
    .map_err(|e| e.to_string())
}
/// # Safety
/// Counts and optional factors reference G×N doubles, optional size factors N
/// doubles and output 999 distinct writable doubles. All borrows are synchronous.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rsdeseq2_r_frozen_local_grid(
    genes: usize,
    samples: usize,
    values: *const f64,
    sf: *const f64,
    nf: *const f64,
    output: *mut f64,
) -> i32 {
    run_backend(|| {
        // SAFETY: these pointers come from the C caller under this export's contract, which satisfies the helper's `# Safety` requirements.
        let counts = CountView::from_doubles(unsafe { counts(genes, samples, values)? })
            .map_err(|e| e.to_string())?;
        let sf = if sf.is_null() {
            None
        } else {
            // SAFETY: the C caller passes a readable buffer of these dimensions that stays valid for this call (see the export's safety contract).
            Some(unsafe { input_slice(sf, samples, "size factors")? })
        };
        let nf = if nf.is_null() {
            None
        } else {
            // SAFETY: these pointers come from the C caller under this export's contract, which satisfies the helper's `# Safety` requirements.
            Some(unsafe { self::counts(genes, samples, nf)? })
        };
        let grid = rsdeseq2::transform::frozen_local_vst_grid_from_counts(counts, sf, nf)
            .map_err(|e| e.to_string())?;
        // SAFETY: the C caller passes a writable buffer of this length, uniquely borrowed for this call (see the export's safety contract).
        unsafe { output_slice(output, 999, "frozen local grid")? }.write_copy_of_slice(&grid);
        Ok(())
    })
}
/// # Safety
/// The input references G×N readable doubles; C validates complete matrix storage.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rsdeseq2_r_validate_counts(
    genes: usize,
    samples: usize,
    values: *const f64,
) -> i32 {
    run_backend(|| {
        // SAFETY: these pointers come from the C caller under this export's contract, which satisfies the helper's `# Safety` requirements.
        CountView::from_doubles(unsafe { counts(genes, samples, values)? })
            .map_err(|e| e.to_string())?;
        Ok(())
    })
}
/// # Safety
/// Counts reference G×N doubles, dispersions/trends G doubles, optional excluded
/// flags G integers. Output indices reference G distinct writable integers and
/// used one writable usize. No borrowed pointers are retained after return.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rsdeseq2_r_normal_shrinkage_rows(
    genes: usize,
    samples: usize,
    values: *const f64,
    dispersions: *const f64,
    fitted: *const f64,
    excluded: *const i32,
    indices: *mut i32,
    used: *mut usize,
) -> i32 {
    run_backend(|| {
        // SAFETY: these pointers come from the C caller under this export's contract, which satisfies the helper's `# Safety` requirements.
        let rows = CountRows::from_doubles(unsafe { counts(genes, samples, values)? })
            .map_err(|e| e.to_string())?;
        let flags = if excluded.is_null() {
            None
        } else {
            Some(
                // SAFETY: the C caller passes a readable buffer of these dimensions that stays valid for this call (see the export's safety contract).
                RowFlags::from_integers(unsafe { input_slice(excluded, genes, "excluded rows")? })
                    .map_err(|e| e.to_string())?,
            )
        };
        let selected = NormalShrinkageRows::new(
            &rows,
            // SAFETY: the C caller passes a readable buffer of these dimensions that stays valid for this call (see the export's safety contract).
            unsafe { input_slice(dispersions, genes, "dispersions")? },
            // SAFETY: the C caller passes a readable buffer of these dimensions that stays valid for this call (see the export's safety contract).
            unsafe { input_slice(fitted, genes, "fitted dispersions")? },
            flags,
        )
        .map_err(|e| e.to_string())?;
        // SAFETY: the C caller passes a writable buffer of this length, uniquely borrowed for this call (see the export's safety contract).
        let out = unsafe { output_slice(indices, genes, "selected rows")? };
        for (slot, &index) in out.iter_mut().zip(selected.indices()) {
            slot.write(i32::try_from(index + 1).map_err(|e| e.to_string())?);
        }
        // SAFETY: the C caller passes a writable buffer of this length, uniquely borrowed for this call (see the export's safety contract).
        let length = unsafe { output_slice(used, 1, "selected row count")? };
        length[0].write(selected.indices().len());
        Ok(())
    })
}
/// # Safety
/// Optional base factors reference G×N doubles, size factors N doubles, offsets
/// N doubles; the distinct writable output references G×N doubles. C validates
/// all lengths and keeps inputs alive throughout this synchronous call.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rsdeseq2_r_offset_factors(
    genes: usize,
    samples: usize,
    sf: *const f64,
    nf: *const f64,
    offsets: *const f64,
    output: *mut f64,
) -> i32 {
    run_backend(|| {
        let sf = if sf.is_null() {
            None
        } else {
            // SAFETY: the C caller passes a readable buffer of these dimensions that stays valid for this call (see the export's safety contract).
            Some(unsafe { input_slice(sf, samples, "size factors")? })
        };
        let nf = if nf.is_null() {
            None
        } else {
            // SAFETY: these pointers come from the C caller under this export's contract, which satisfies the helper's `# Safety` requirements.
            Some(unsafe { counts(genes, samples, nf)? })
        };
        // SAFETY: the C caller passes a readable buffer of these dimensions that stays valid for this call (see the export's safety contract).
        let offsets = unsafe { input_slice(offsets, samples, "formula offsets")? };
        let length = genes
            .checked_mul(samples)
            .ok_or("normalization offset shape overflow")?;
        // SAFETY: the C caller passes a writable buffer of this length, uniquely borrowed for this call (see the export's safety contract).
        let out = unsafe { output_slice(output, length, "effective factors")? };
        visit_offset_factors_column_major(genes, samples, sf, nf, offsets, |index, value| {
            out[index].write(value);
        })
        .map_err(|error| error.to_string())
    })
}
