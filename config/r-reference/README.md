# Pinned independent R environment

The environment JSON specifies R, Bioconductor and reference-package versions.
The source TSV pins package source URLs and hashes. Keep their original dated
identities when reproducing or auditing saved numerical references.

Install/check through `scripts/environment/install_r_reference.py` and its R
helper. CI caches include both lock files and the installer implementations.
Benchmark outputs and validation measurements belong in `docs/data/`, not here.
