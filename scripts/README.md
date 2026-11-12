# Repository tooling

Run commands from the repository root. Python tools require Python 3.10 or
newer. Reference generation and R interface checks use the pinned
`deseq2-r462-native` environment; production numerical code remains in Rust.

| Folder | Responsibility |
| --- | --- |
| `references/` | Independent R/mathematical reference generators and deterministic fixture/contrast generation. The DESeq2 panels live together in `references/deseq2/`. |
| `validation/` | Numerical/interface comparisons, boundary audits, result scoring and target checks. Shared R check helpers live beside their callers. |
| `benchmarks/` | Matched time/RSS measurements and workflow profiling. |
| `environment/` | Hash-verified installation and validation of the pinned R environment. |
| `release/` | Vendoring, checked-source manifests, package staging and binary provenance. These helpers do not authorize publication. |
| `tests/` | Tests of tooling behavior and acceptance gates. |

Examples of maintainer commands:

```sh
scripts/release/sync_r_vendor.sh --check
python3 -m unittest discover -s scripts/tests -p 'test_*.py'
python3 scripts/environment/install_r_reference.py --help
python3 scripts/references/generate_gtex_mode_contrasts.py --help
```

Use fresh output directories below `results/`. Reference generation is distinct
from acceptance: normal validation replays frozen expected values and preserves
their hashes. See [reproducibility](../docs/reproducibility.md).
