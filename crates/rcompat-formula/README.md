# rcompat-formula

Build R-style statistical design matrices in Rust. Turn a model formula and
sample metadata into named columns for model fitting, including factor
contrasts, interactions, numeric transforms, and offsets.

## Install

Requires Rust 1.96 or newer.

```sh
cargo add rcompat-formula@0.1.1
```

## Example

Build a two-group design with `control` as the reference level. The resulting
matrix has an intercept and an indicator for the treated samples.

```rust
use rcompat_formula::design::{
    FormulaFactorColumn, FormulaModelFrame, formula_design_for_fitting_from_model_frame,
};

fn main() -> Result<(), rcompat_formula::Error> {
    let frame = FormulaModelFrame {
        factors: vec![FormulaFactorColumn {
            name: "condition".into(),
            sample_levels: ["control", "control", "treated", "treated"]
                .map(String::from)
                .to_vec(),
            levels: Some(vec!["control".into(), "treated".into()]),
            reference: Some("control".into()),
        }],
        numeric_covariates: vec![],
    };

    let fit_design = formula_design_for_fitting_from_model_frame("~ condition", &frame)?;
    let design = fit_design.design.standard_design;
    println!("{:?}", design.coefficient_names());
    println!("{:?}", design.matrix().as_slice());
    Ok(())
}
```

## What it provides

- Standard treatment-coded and expanded factor designs.
- Formula terms and interactions, with built-in transforms such as
  `log()`, `scale()`, and `poly()`.
- Separate per-sample offsets and model-frame metadata for downstream fitting.
- C-locale collation by default and locale-aware collation with the optional
  `icu` feature.

Formula expressions use supported built-ins and caller-supplied values. The
crate runs without an R interpreter.

[API documentation](https://docs.rs/rcompat-formula) ·
[License](https://github.com/deminden/rsdeseq2/blob/rcompat-v0.1.1/crates/rcompat-formula/LICENSE) ·
[Third-party notices](https://github.com/deminden/rsdeseq2/blob/rcompat-v0.1.1/crates/rcompat-formula/THIRD_PARTY_NOTICES.md)
