# Documentation

These guides describe the 0.6.2 source and its supported behavior.

| Task | Guide |
| --- | --- |
| Run file-based analysis | [CLI](cli.md) |
| Use R fits, results, transforms and shrinkage | [R wrapper](r-wrapper.md) |
| Check supported methods and limits | [Compatibility](compatibility.md), [numerical ownership contract](native-supported-contract.md) |
| Understand statistical implementation | [Algorithms](algorithms.md), [platform arithmetic](platform_arithmetic.md), [summation rules](summation-rules.md), [DESeq2 quirks](deseq2-quirks.md) |
| Review precision, speed and memory | [Current validation](testing-parity.md#ci-status), [benchmark procedure](benchmarks.md) |
| Reproduce data and independent references | [Reproducibility](reproducibility.md) |
| Run tests or review validation | [Testing and parity](testing-parity.md) |
| Read development history | [Hardening progress](reports/hardening-progress.md), [0.5.0 baseline](reports/baseline-0.5.0.md), [GTEx comparison](reports/gtex-native-wald-benchmark.md) |
| Change or package the code | [Development](development.md), [release procedure](release.md), [0.6.2 notes](releases/0.6.2.md) |

`data/` contains the raw evidence used by these guides and active validation
tools. Operational R environment locks live in [`config/r-reference/`](../config/r-reference/).
Reference metadata and measurement source hashes identify the exact environments
and binaries that produced their values.
