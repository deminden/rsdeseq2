#include <R.h>
#include <Rinternals.h>
#include <math.h>
#include <stdint.h>
#include <stddef.h>
#include <limits.h>
#include <string.h>
#include "validation.h"
#include "contrast_inference.h"

extern const char *rsdeseq2_r_last_error(void);

SEXP rsdeseq2_contrast_batch(SEXP beta, SEXP se, SEXP covariance, SEXP requests) {
    rsdeseq2_validate_inference_matrices(beta, se, covariance);
    if (TYPEOF(requests) != VECSXP || XLENGTH(requests) > RSDESEQ2_MAX_CONTRAST_BATCH)
        Rf_error("contrast batch must be a list of at most 64 requests");
    R_xlen_t genes = matrix_rows(beta), coefficients = matrix_cols(beta), count = XLENGTH(requests);
    if (count == 0) return Rf_allocVector(VECSXP, 0);
    checked_product(checked_product(genes, count), 4);
    checked_product(count, sizeof(Rsdeseq2ContrastRequest));
    checked_product(count, sizeof(Rsdeseq2ContrastZero));
    Rsdeseq2ContrastRequest *native = (Rsdeseq2ContrastRequest *)R_alloc((size_t)count, sizeof(*native));
    Rsdeseq2ContrastZero *zero = (Rsdeseq2ContrastZero *)R_alloc((size_t)count, sizeof(*zero));
    memset(native, 0, (size_t)count * sizeof(*native));
    memset(zero, 0, (size_t)count * sizeof(*zero));
    for (R_xlen_t i = 0; i < count; ++i) {
        SEXP item = VECTOR_ELT(requests, i);
        SEXP contrast = named_list_field(item, "contrast", "contrast batch request");
        SEXP wald = named_list_field(item, "wald", "contrast batch request");
        SEXP context = named_list_field(item, "zero", "contrast batch request");
        validate_vector(contrast, REALSXP, coefficients, "contrast");
        native[i].contrast = REAL(contrast);
        if (wald != R_NilValue) {
            rsdeseq2_validate_wald_transport(
                named_list_field(wald, "threshold", "Wald request"),
                named_list_field(wald, "alternative", "Wald request"),
                named_list_field(wald, "useT", "Wald request"),
                named_list_field(wald, "df", "Wald request"), genes, &native[i]);
        }
        if (context != R_NilValue) {
            rsdeseq2_validate_contrast_zero(context, genes, coefficients, &zero[i]);
            native[i].zero = &zero[i];
        }
    }
    SEXP result = PROTECT(Rf_allocVector(VECSXP, count));
    for (R_xlen_t i = 0; i < count; ++i) {
        SEXP item = PROTECT(rsdeseq2_allocate_contrast_output(genes, &native[i]));
        SET_VECTOR_ELT(result, i, item);
        UNPROTECT(1);
    }
    if (count) {
        int status = rsdeseq2_r_contrast_batch((size_t)genes, (size_t)coefficients, REAL(beta),
            se != R_NilValue, se == R_NilValue ? NULL : REAL(se), covariance != R_NilValue,
            covariance == R_NilValue ? NULL : REAL(covariance), (size_t)count, native, NA_REAL);
        if (status) Rf_error("rsdeseq2 Rust backend: %s", rsdeseq2_r_last_error());
    }
    UNPROTECT(1);
    return result;
}
