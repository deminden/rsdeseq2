#ifndef RSDESEQ2_CONTRAST_ZERO_H
#define RSDESEQ2_CONTRAST_ZERO_H
#include <Rinternals.h>
#include <stddef.h>
#include <stdint.h>

/* All pointed-to buffers are validated and borrowed for one synchronous call.
 * The counts/design use column-major storage. Rules: 1 numeric, 2 factor. */
typedef struct {
    size_t samples;
    const double *counts;
    const double *design;
    uint8_t has_design;
    int32_t rule;
    const char *const *sample_levels;
    const char *numerator;
    const char *denominator;
} Rsdeseq2ContrastZero;

/* Called entirely before entering Rust, including UTF-8 translation. */
void rsdeseq2_validate_contrast_zero(SEXP, R_xlen_t, R_xlen_t, Rsdeseq2ContrastZero *);
#endif
