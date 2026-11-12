# Documentation map

The files in this directory are maintained user and contributor guides.

- [CLI](cli.md), [R wrapper](r-wrapper.md), [supported contract](native-supported-contract.md)
- [Compatibility](compatibility.md), [gap analysis](deseq2-gap-analysis.md), [algorithms](algorithms.md)
- [Testing and parity](testing-parity.md), [reproducibility](reproducibility.md), [benchmarks](benchmarks.md)
- [Development](development.md), [release procedure](release.md)

`reports/` holds the current empirical comparison and latest correction evidence.
`history/` preserves superseded development plans, implementation checkpoints and
older benchmark reports. Their source revisions, counts and measurements retain
their original scope. `releases/` contains versioned release notes. `data/` holds
supporting measurement and reference evidence; operational environment locks live
separately in [`config/r-reference/`](../config/r-reference/).

Historical JSON/TSV provenance records retain original command paths and hashes.
Current executable commands and documentation links use the reorganised paths.
