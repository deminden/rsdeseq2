#ifndef RSDESEQ2_RESULT_POLICY_H
#define RSDESEQ2_RESULT_POLICY_H
#include <stddef.h>
#include <stdint.h>

typedef struct { const int32_t *codes; size_t levels; } Rsdeseq2FactorCodes;
typedef struct {
    const double *max_cooks;
    const int32_t *refitted;
    const int32_t *new_zero;
    size_t samples;
    size_t coefficients;
    const double *design;
    const double *counts;
    const double *cooks;
    size_t factors_len;
    const Rsdeseq2FactorCodes *factors;
} Rsdeseq2ResultContext;

extern int rsdeseq2_r_present_results(size_t, double *, double *, double *, double *,
    uint8_t, double, const Rsdeseq2ResultContext *, uint8_t *);
extern int rsdeseq2_r_result_replacement_storage(size_t, const int32_t *, const int32_t *, uint8_t *);
#endif
