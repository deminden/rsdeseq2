//! Standalone probability primitives shared by inference and row policy.

mod f_quantile;
mod gamma;
mod normal;
mod normal_log;
mod quadrature;
mod student_t;

pub(crate) use f_quantile::f_99_quantile;
pub(crate) use gamma::large_shape_chi_square_tail;
pub(crate) use normal::{normal_two_sided_tail, normal_upper_tail};
pub(crate) use normal_log::normal_log_upper_tail;
pub(crate) use student_t::student_t_upper_tail;
