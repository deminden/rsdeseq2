//! Formula term expansion edge cases.

use crate::support::assert_name_not_found;
use rsdeseq2::errors::NameScope;
use rsdeseq2::prelude::*;

fn frame() -> FormulaModelFrame {
    FormulaModelFrame {
        factors: vec![FormulaFactorColumn {
            name: "condition".to_string(),
            sample_levels: ["A", "A", "B", "B"].map(String::from).to_vec(),
            levels: None,
            reference: None,
        }],
        numeric_covariates: vec![FormulaNumericColumn {
            name: "age".to_string(),
            values: vec![30.0, 41.0, 52.0, 63.0],
        }],
    }
}

/// Unknown calls used to bounce between two expansion functions until the
/// stack overflowed, aborting the process (found by the formula fuzz target).
/// Run on a small stack so a regression fails instead of aborting.
#[test]
fn unknown_calls_are_errors_not_stack_overflows() {
    std::thread::Builder::new()
        .stack_size(1 << 20)
        .spawn(|| {
            for (formula, missing) in [
                ("~ q(age)", "q(age)"),
                ("~q()", "q()"),
                ("~ f() + condition", "f()"),
                ("~ q(age):condition", "q(age)"),
            ] {
                let error =
                    expanded_formula_design_with_offsets_from_model_frame(formula, &frame())
                        .err()
                        .unwrap_or_else(|| panic!("{formula} should be rejected"));
                assert_name_not_found(
                    error,
                    ErrorKind::Formula,
                    NameScope::FormulaVariable,
                    missing,
                );
            }
            assert!(
                expanded_formula_design_with_offsets_from_model_frame("~ log(age)", &frame())
                    .is_ok()
            );
        })
        .unwrap()
        .join()
        .unwrap();
}
