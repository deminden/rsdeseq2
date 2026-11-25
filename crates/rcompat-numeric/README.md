# rcompat-numeric

Numerical building blocks for R-compatible statistical software in Rust.
Use it to compute probability tails and quantiles, evaluate special
functions, or share checked arithmetic and dense matrices across an analysis.

## Install

Requires Rust 1.96 or newer.

```sh
cargo add rcompat-numeric@0.2.0
```

## Example

Compute the upper-tail probability for a standard normal statistic with
`z = 8`. Computing the tail directly preserves small probabilities.

```rust
use rcompat_numeric::probability::{normal_log_upper_tail, normal_upper_tail};

fn main() {
    let z = 8.0;
    println!("P(Z > {z}) = {:.3e}", normal_upper_tail(z));
    println!("log P(Z > {z}) = {:.3}", normal_log_upper_tail(z));
}
```

## What it provides

- Normal and Student's t tails, normal quantiles, and special functions.
- Profile-aware elementary functions and reductions. The default portable
  profile uses Rust math. The R-host profile follows system math and the
  reduction schedule for the host. Both run without an R installation.
- A checked row-major matrix, borrowed views, and structured errors.

[API documentation](https://docs.rs/rcompat-numeric) ·
[MIT license](LICENSE) ·
[Third-party notices](THIRD_PARTY_NOTICES.md)
