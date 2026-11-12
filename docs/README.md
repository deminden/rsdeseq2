# Documentation

These guides describe the current 0.5.0 source and its supported behavior.

| Task | Guide |
| --- | --- |
| Run file-based analysis | [CLI](cli.md) |
| Use R fits, results, transforms and shrinkage | [R wrapper](r-wrapper.md) |
| Check supported methods and limits | [Compatibility](compatibility.md), [numerical ownership contract](native-supported-contract.md) |
| Understand statistical implementation | [Algorithms](algorithms.md), [platform arithmetic](linux_numerical_parity.md) |
| Review precision, speed and memory | [GTEx comparison](reports/gtex-native-wald-benchmark.md), [benchmark procedure](benchmarks.md) |
| Reproduce data and independent references | [Reproducibility](reproducibility.md) |
| Run tests or inspect passing CI | [Testing and parity](testing-parity.md) |
| Change or package the code | [Development](development.md), [release procedure](release.md), [0.5.0 notes](releases/0.5.0.md) |

`data/` contains the raw evidence used by these guides and active validation
tools. Operational R environment locks live in [`config/r-reference/`](../config/r-reference/).
Reference metadata and measurement source hashes identify the exact environments
and binaries that produced their values.
