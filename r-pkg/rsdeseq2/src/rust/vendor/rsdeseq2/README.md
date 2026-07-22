# Rust Core Packaging Directory

The R build currently excludes the full numerical Rust core and provides
selected primitive `.Call` bridges elsewhere in the package. This directory is
reserved for the core sources needed by a CRAN/Bioconductor-compatible build.
