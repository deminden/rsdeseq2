mod common;

use common::*;
use rsdeseq2::prelude::*;

use common::assert_reference_value as assert_f64_or_missing;

#[path = "dispersion_reference/fixed.rs"]
mod fixed;

#[path = "dispersion_reference/weighted_map.rs"]
mod weighted_map;

#[path = "dispersion_reference/map.rs"]
mod map;

#[path = "dispersion_reference/gene_wise.rs"]
mod gene_wise;

#[path = "dispersion_reference/normalization_factors.rs"]
mod normalization_factors;

#[path = "dispersion_reference/trends.rs"]
mod trends;

#[path = "dispersion_reference/priors.rs"]
mod priors;
