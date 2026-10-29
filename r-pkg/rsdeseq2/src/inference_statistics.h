#ifndef RSDESEQ2_INFERENCE_STATISTICS_H
#define RSDESEQ2_INFERENCE_STATISTICS_H
#include <stddef.h>
#include <stdint.h>
/* Borrowed primitive buffers and caller-allocated outputs; no R calculations. */
extern int rsdeseq2_r_wald_statistics(size_t, const double *, const double *, double,
    int, uint8_t, size_t, const double *, double, double *, double *);
extern int rsdeseq2_r_lrt_pvalues(size_t, const double *, size_t, const double *, double, double *);
extern int rsdeseq2_r_default_cooks_cutoff(size_t, size_t, double *, uint8_t *);
#endif
