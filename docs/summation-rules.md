# Summation rules

Which reduction each R primitive imitation uses. Changing a rule is a
NUMERIC task: it moves results and needs reference evidence.

| R primitive | Used for | Rust rule | Location |
| --- | --- | --- | --- |
| `rowMeans(log(counts))` | Size-factor geometric means | Compensated (Neumaier) sum of logs, then divide by n | `normalization::checked_compensated_sum` |
| `rowMeans(normalized)` | `baseMean` | Scale-normalised mean `((Σ x/s)/n)·s`, `s = max|x|` | `numeric::checked::mean` |
| `rowVars` / `var` | `baseVar` | Two-pass: scaled centred squares, `(Σ((x−m)/s)²/(n−1))` rescaled | `normalization::checked_sample_variance` |
| `median` | Size factors | Sort, middle value or `a/2 + b/2` | `math::median` |
| `quantile(type = 7)` | Beta-prior variance, local VST, filtering thresholds | `lo + f·(hi − lo)`; beta-prior keeps equal endpoints unchanged | `quantile_type7_*` |
| `mean(x, trim)` | Cook's replacement counts | Plain sum of the central values / count; median for `trim >= 0.5` | `cooks::r_trimmed_mean` |
| `sum` | Likelihoods, Wald contrasts | Sequential sum, or scale-normalised sum where overflow is possible | `numeric::checked::{sum, scaled_sum}` |

## Known mismatch

`baseMean` imitates `rowMeans` with the scale-normalised mean, while the
geometric means imitate `rowMeans` with a compensated sum. R uses the same
extended-precision accumulator for both. The scaled mean is kept because it
matches the committed DESeq2 references within the documented bounds; moving
`baseMean` to the compensated rule is a NUMERIC candidate that needs the
Phase 6.4 parity metrics (log10 differences and significance flips) before and
after.
