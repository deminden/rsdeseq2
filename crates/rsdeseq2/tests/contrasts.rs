use rsdeseq2::prelude::*;

fn named_design() -> DesignMatrix {
    DesignMatrix::from_row_major(
        4,
        3,
        vec![
            1.0, 0.0, 0.0, //
            1.0, 0.0, 1.0, //
            1.0, 1.0, 0.0, //
            1.0, 1.0, 1.0,
        ],
        Some(vec![
            "Intercept".into(),
            "condition_B_vs_A".into(),
            "batch_Y_vs_X".into(),
        ]),
    )
    .unwrap()
}

#[path = "contrasts/all_zero.rs"]
mod all_zero;
#[path = "contrasts/factor_levels.rs"]
mod factor_levels;
#[path = "contrasts/model_frame.rs"]
mod model_frame;
#[path = "contrasts/result_semantics.rs"]
mod result_semantics;
#[path = "contrasts/selectors.rs"]
mod selectors;
