#ifndef RSDESEQ2_CONTRAST_INFERENCE_H
#define RSDESEQ2_CONTRAST_INFERENCE_H
#include <Rinternals.h>
#include <stddef.h>
#include <stdint.h>
#include "contrast_zero.h"
#define RSDESEQ2_MAX_CONTRAST_BATCH 64

typedef struct {
    const double *contrast;
    double threshold;
    int32_t alternative;
    uint8_t use_wald;
    uint8_t use_t;
    size_t df_len;
    const double *df;
    const Rsdeseq2ContrastZero *zero;
    double *lfc;
    double *lfc_se;
    double *stat;
    double *pvalue;
} Rsdeseq2ContrastRequest;

extern int rsdeseq2_r_contrast_batch(size_t, size_t, const double *, uint8_t, const double *,
    uint8_t, const double *, size_t, const Rsdeseq2ContrastRequest *, double);
void rsdeseq2_validate_inference_matrices(SEXP, SEXP, SEXP);
void rsdeseq2_validate_wald_transport(SEXP, SEXP, SEXP, SEXP, R_xlen_t, Rsdeseq2ContrastRequest *);
SEXP rsdeseq2_allocate_contrast_output(R_xlen_t, Rsdeseq2ContrastRequest *);
#endif
