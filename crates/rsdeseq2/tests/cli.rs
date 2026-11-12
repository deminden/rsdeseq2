mod common;
use rsdeseq2::prelude::*;

use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicUsize, Ordering};

#[path = "cli/support.rs"]
mod support;
use support::*;

#[path = "cli/contrast_options.rs"]
mod contrast_options;
#[path = "cli/normalization_options.rs"]
mod normalization_options;
#[path = "cli/numerical_assertions.rs"]
mod numerical_assertions;

#[path = "cli/beta_prior_options.rs"]
mod beta_prior_options;
#[path = "cli/filter_options.rs"]
mod filter_options;
#[path = "cli/statistical_options.rs"]
mod statistical_options;

#[path = "cli/outputs.rs"]
mod outputs;

#[path = "cli/beta_prior_outputs.rs"]
mod beta_prior_outputs;

#[path = "cli/beta_prior_validation.rs"]
mod beta_prior_validation;

#[path = "cli/contrast_validation.rs"]
mod contrast_validation;

#[path = "cli/transforms.rs"]
mod transforms;
