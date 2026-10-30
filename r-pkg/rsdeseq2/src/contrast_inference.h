#ifndef RSDESEQ2_CONTRAST_INFERENCE_H
#define RSDESEQ2_CONTRAST_INFERENCE_H
#include <stddef.h>
#include <stdint.h>
#include "contrast_zero.h"
extern int rsdeseq2_r_contrast_inference(size_t, size_t, const double *,
    uint8_t, const double *, uint8_t, const double *, const double *, uint8_t,
    double, int, uint8_t, size_t, const double *, double,
    double *, double *, double *, double *, const Rsdeseq2ContrastZero *);
#endif
