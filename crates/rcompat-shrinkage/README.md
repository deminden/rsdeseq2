# rcompat-shrinkage

Empirical Bayes shrinkage in Rust. Fit ashr-style adaptive mixtures from
effect estimates and their standard errors, or use apeglm-style heavy-tailed
priors for supported regression likelihoods.

## Install

Requires Rust 1.96 or newer.

```sh
cargo add rcompat-shrinkage@0.1.2
```

## Example

Compute a posterior for an estimate of 1.2 with a standard error of 0.5.
This example supplies a fixed prior with a point mass at zero and a
unit-scale normal component.

```rust
use rcompat_shrinkage::shrinkage::NormalMixture;

fn main() -> Result<(), rcompat_shrinkage::Error> {
    let prior = NormalMixture::new(vec![0.3, 0.7], vec![0.0, 0.0], vec![0.0, 1.0])?;
    let posterior = prior.posterior(1.2, 0.5)?;

    println!("Posterior mean: {:.3}", posterior.mean());
    println!("Posterior SD: {:.3}", posterior.standard_deviation());
    Ok(())
}
```

## What it provides

- Fixed and adaptive mixture priors, with posterior means, standard
  deviations, and false-sign summaries for supported models.
- Heavy-tailed coefficient priors and posterior fitting for supported
  regression likelihoods, including negative-binomial counts.
- Native Rust execution and structured errors, with no R runtime dependency.

[API documentation](https://docs.rs/rcompat-shrinkage) ·
[MIT license](LICENSE) ·
[Third-party notices](THIRD_PARTY_NOTICES.md)
