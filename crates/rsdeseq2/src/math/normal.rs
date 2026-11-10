/// Log probability that a standard normal variable exceeds `q`.
///
/// Evaluates far tails directly in log space, including tails whose ordinary
/// probabilities underflow. Infinite arguments have their limiting values;
/// NaN propagates. This does not change the existing Wald probability path.
pub fn normal_log_upper_tail(q: f64) -> f64 {
    crate::probability::normal_log_upper_tail(q)
}

/// Standard normal log CDF, retaining small log probabilities near one.
pub fn normal_log_cdf(q: f64) -> f64 {
    normal_log_upper_tail(-q)
}
