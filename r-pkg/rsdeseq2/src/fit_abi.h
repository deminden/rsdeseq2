#ifndef RSDESEQ2_FIT_ABI_H
#define RSDESEQ2_FIT_ABI_H
#include <stddef.h>
#include <stdint.h>

/* Opaque, Rust-owned replacement counts. No R pointer is retained. A fit
 * publishes at most one owner, only after successful computation. NULL means
 * there is no replacement-count output to retain (nRefit == 0). C queries and
 * copies a live owner once, then drops it on ordinary return or R unwinding. */
typedef struct Rsdeseq2ReplacementCounts Rsdeseq2ReplacementCounts;

extern int rsdeseq2_r_fit(
    size_t genes, size_t samples, const double *counts,
    const double *size_factors, const double *normalization, const double *weights,
    size_t coefficients, const double *design,
    size_t reduced_coefficients, const double *reduced,
    int fit_type, int test_type, size_t min_replicates, double replacement_cutoff,
    size_t threads, int execution_mode,
    double *results, double *beta, double *beta_se, double *covariance,
    double *mu, double *cooks, Rsdeseq2ReplacementCounts **replacement,
    uint8_t *replace, uint8_t *replaceable, uint8_t *refit_rows,
    uint8_t *new_all_zero, size_t *n_refit,
    double *factors_out, double *dispersions, double *disp_fit,
    double *base_var, int *all_zero);

/* Shape/copy require a live owner returned by fit. Output allocations must be
 * validated by C; copy also checks that the requested shape matches the owner.
 * Drop consumes the owner and clears its slot; a NULL slot/owner is harmless. */
extern int rsdeseq2_r_replacement_shape(const Rsdeseq2ReplacementCounts *, size_t *, size_t *);
extern int rsdeseq2_r_replacement_copy(const Rsdeseq2ReplacementCounts *, size_t, size_t, double *);
extern void rsdeseq2_r_replacement_drop(Rsdeseq2ReplacementCounts **);
#endif
