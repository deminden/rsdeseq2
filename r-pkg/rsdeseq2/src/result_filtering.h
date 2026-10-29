#ifndef RSDESEQ2_RESULT_FILTERING_H
#define RSDESEQ2_RESULT_FILTERING_H
#include <stddef.h>
#include <stdint.h>
/* Synchronous primitive-array bridge; no callback, ownership transfer or R
 * arithmetic. Caller allocates G output cells and K metadata cells as stated. */
extern int rsdeseq2_r_bh_adjust(size_t, const double *, double, double *);
extern int rsdeseq2_r_independent_filtering(size_t, const double *, const double *,
    size_t, const double *, uint8_t, double, double, double *, uint8_t *,
    double *, int32_t *, double *, int32_t *, double *);
#endif
