//! Synchronous borrowed fitted-result policy. No R APIs or retained pointers.
use super::{borrowed_cooks_heuristic, input_slice, output_slice, run_backend};
use rsdeseq2::matrix::{MatrixOrder, MatrixView};
use rsdeseq2::results::{CooksPolicy, ReplacementRows, RowFlags, two_level_design_matches};

#[repr(C)]
pub(super) struct FactorCodes {
    codes: *const i32,
    levels: usize,
}
#[repr(C)]
pub(super) struct ResultContext {
    max_cooks: *const f64,
    refitted: *const i32,
    new_zero: *const i32,
    samples: usize,
    coefficients: usize,
    design: *const f64,
    counts: *const f64,
    cooks: *const f64,
    factors_len: usize,
    factors: *const FactorCodes,
}

unsafe fn replacement_rows<'a>(
    genes: usize,
    refitted: *const i32,
    new_zero: *const i32,
) -> Result<ReplacementRows<'a>, String> {
    let flags = |pointer: *const i32| -> Result<_, String> {
        if pointer.is_null() {
            Ok(None)
        } else {
            RowFlags::from_integers(unsafe { input_slice(pointer, genes, "replacement flags")? })
                .map(Some)
                .map_err(|e| e.to_string())
        }
    };
    ReplacementRows::new(genes, flags(refitted)?, flags(new_zero)?).map_err(|e| e.to_string())
}

/// # Safety
/// Optional row indicators reference G readable integers; storage points to one
/// writable byte. C keeps every buffer alive during this synchronous query.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rsdeseq2_r_result_replacement_storage(
    genes: usize,
    refitted: *const i32,
    new_zero: *const i32,
    storage: *mut u8,
) -> i32 {
    run_backend(|| {
        let rows = unsafe { replacement_rows(genes, refitted, new_zero)? };
        let output = unsafe { output_slice(storage, 1, "replacement storage flag")? };
        output[0] = u8::from(rows.has_new_zero());
        Ok(())
    })
}

/// # Safety
/// C validates the context descriptor, optional G row flags/maxima, N×P design,
/// G×N count/Cook's matrices, F factor descriptors and N code buffers. When
/// replacement zero rows exist, all four G output columns are distinct initialized
/// copies disjoint from inputs. Otherwise effect/SE/stat may be null; p-values
/// may be null only when filtering is disabled. G output flags are writable.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rsdeseq2_r_present_results(
    genes: usize,
    effect: *mut f64,
    se: *mut f64,
    stat: *mut f64,
    pvalue: *mut f64,
    active: u8,
    cutoff: f64,
    context: *const ResultContext,
    outlier: *mut u8,
) -> i32 {
    run_backend(|| {
        if active > 1 {
            return Err("Cook's mode must be zero or one".into());
        }
        let context = unsafe { context.as_ref() }.ok_or("result context is null")?;
        let rows = unsafe { replacement_rows(genes, context.refitted, context.new_zero)? };
        let needs_zero = rows.has_new_zero();
        // C only copies columns that Rust may modify. Untouched columns remain
        // shared R inputs and are never borrowed mutably here.
        let effect =
            unsafe { output_slice(effect, if needs_zero { genes } else { 0 }, "result effects")? };
        let se = unsafe { output_slice(se, if needs_zero { genes } else { 0 }, "result SEs")? };
        let stat = unsafe {
            output_slice(
                stat,
                if needs_zero { genes } else { 0 },
                "result statistics",
            )?
        };
        let pvalue = unsafe {
            output_slice(
                pvalue,
                if needs_zero || active == 1 { genes } else { 0 },
                "result p-values",
            )?
        };
        let outlier = unsafe { output_slice(outlier, genes, "result outlier flags")? };
        // Disabled filtering preserves untouched IEEE NaN/NA payloads.
        if active == 0 {
            if needs_zero {
                rows.apply_columns(effect, se, stat, pvalue)
                    .map_err(|e| e.to_string())?;
            }
            outlier.fill(u8::MAX);
            return Ok(());
        }
        let factors =
            unsafe { input_slice(context.factors, context.factors_len, "Cook's factors")? };
        let mut low_count = false;
        if !context.design.is_null() && context.coefficients == 2 {
            let length = context
                .samples
                .checked_mul(context.coefficients)
                .ok_or("Cook's design shape overflow")?;
            let design = MatrixView::new(
                context.samples,
                context.coefficients,
                unsafe { input_slice(context.design, length, "Cook's design")? },
                MatrixOrder::ColumnMajor,
            )
            .map_err(|e| e.to_string())?;
            for factor in factors {
                if factor.levels != 2 {
                    continue;
                }
                let codes =
                    unsafe { input_slice(factor.codes, context.samples, "Cook's factor codes")? };
                // Missing factor observations cannot equal a complete finite
                // design partition; valid factor metadata uses one-based codes.
                if codes.contains(&i32::MIN) {
                    continue;
                }
                if codes.iter().any(|&v| v < 1 || v as usize > factor.levels) {
                    return Err("invalid Cook's factor code".into());
                }
                if two_level_design_matches(design, codes, factor.levels)
                    .map_err(|e| e.to_string())?
                {
                    low_count = true;
                    break;
                }
            }
        }
        let heuristic = if low_count {
            let length = genes
                .checked_mul(context.samples)
                .ok_or("Cook's shape overflow")?;
            let values = unsafe { input_slice(context.cooks, length, "Cook's distances")? };
            if values.iter().any(|v| v.is_infinite()) {
                return Err("cooks must contain finite values or NA".into());
            }
            Some(unsafe {
                borrowed_cooks_heuristic(genes, context.samples, context.counts, context.cooks)?
            })
        } else {
            None
        };
        let maximum = unsafe { input_slice(context.max_cooks, genes, "maximum Cook's distances")? };
        let policy = CooksPolicy::new(Some(cutoff), heuristic)
            .map_err(|e| e.to_string())?
            .with_replacement(rows);
        if needs_zero {
            rows.apply_columns(effect, se, stat, pvalue)
                .map_err(|e| e.to_string())?;
        }
        if pvalue
            .iter()
            .any(|v| !v.is_nan() && (!v.is_finite() || !(0.0..=1.0).contains(v)))
        {
            return Err("pvalue must contain values in [0, 1] or NA".into());
        }
        let option = |value: f64| if value.is_nan() { None } else { Some(value) };
        let mask = policy
            .evaluate(
                pvalue.iter().copied().map(option),
                maximum.iter().copied().map(option),
            )
            .map_err(|e| e.to_string())?;
        for gene in 0..genes {
            // Retain the legacy active-mask missing-value encoding.
            pvalue[gene] = mask.pvalue[gene].unwrap_or(f64::NAN);
            outlier[gene] = mask.outlier[gene].map_or(u8::MAX, u8::from);
        }
        Ok(())
    })
}
