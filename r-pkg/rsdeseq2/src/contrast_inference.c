#include <R.h>
#include <Rinternals.h>
#include <math.h>
#include <stdint.h>
#include <stddef.h>
#include <limits.h>
#include "validation.h"
#include "contrast_inference.h"

extern const char *rsdeseq2_r_last_error(void);

void rsdeseq2_validate_inference_matrices(SEXP beta, SEXP se, SEXP covariance) {
    validate_matrix(beta, "beta");
    R_xlen_t genes = matrix_rows(beta), coefficients = matrix_cols(beta);
    if (se != R_NilValue) {
        validate_matrix(se, "betaSE");
        if (matrix_rows(se) != genes || matrix_cols(se) != coefficients)
            Rf_error("betaSE dimensions must match beta");
    }
    if (covariance != R_NilValue) {
        SEXP dims = Rf_getAttrib(covariance, R_DimSymbol);
        if (TYPEOF(dims) != INTSXP || XLENGTH(dims) != 3 ||
            INTEGER(dims)[0] != genes || INTEGER(dims)[1] != coefficients ||
            INTEGER(dims)[2] != coefficients)
            Rf_error("beta covariance dimensions must be genes x coefficients x coefficients");
        validate_vector(covariance, REALSXP,
            checked_product(genes, checked_product(coefficients, coefficients)), "beta covariance");
    }
}

void rsdeseq2_validate_wald_transport(SEXP threshold, SEXP alternative, SEXP use_t, SEXP df,
                                     R_xlen_t genes, Rsdeseq2ContrastRequest *request) {
    validate_real_scalar(threshold, "lfcThreshold");
    request->alternative = validate_integer_scalar(alternative, 0, 5, "alternative");
    validate_boolean(use_t, "useT");
    request->use_t = (uint8_t)LOGICAL(use_t)[0];
    request->use_wald = 1;
    request->threshold = Rf_asReal(threshold);
    if (request->use_t) {
        validate_vector(df, REALSXP, XLENGTH(df), "degrees of freedom");
        if (XLENGTH(df) != 1 && XLENGTH(df) != genes)
            Rf_error("degrees of freedom must be scalar or one per gene");
        checked_product(XLENGTH(df), 1);
        request->df_len = (size_t)XLENGTH(df);
        request->df = REAL(df);
    } else if (df != R_NilValue) Rf_error("degrees of freedom require useT=TRUE");
}

SEXP rsdeseq2_allocate_contrast_output(R_xlen_t genes, Rsdeseq2ContrastRequest *request) {
    const char *wald_names[] = {"log2FoldChange", "lfcSE", "stat", "pvalue", ""};
    const char *projection_names[] = {"log2FoldChange", "lfcSE", ""};
    SEXP result = PROTECT(Rf_mkNamed(VECSXP, request->use_wald ? wald_names : projection_names));
    for (int column = 0; column < (request->use_wald ? 4 : 2); ++column) {
        SEXP values = PROTECT(Rf_allocVector(REALSXP, genes));
        SET_VECTOR_ELT(result, column, values);
        switch (column) {
            case 0: request->lfc = REAL(values); break;
            case 1: request->lfc_se = REAL(values); break;
            case 2: request->stat = REAL(values); break;
            case 3: request->pvalue = REAL(values); break;
        }
        UNPROTECT(1);
    }
    UNPROTECT(1);
    return result;
}

static SEXP contrast_inference(SEXP beta, SEXP se, SEXP covariance, SEXP contrast,
                              Rsdeseq2ContrastRequest *request, SEXP context) {
    rsdeseq2_validate_inference_matrices(beta, se, covariance);
    R_xlen_t genes = matrix_rows(beta), coefficients = matrix_cols(beta);
    validate_vector(contrast, REALSXP, coefficients, "contrast");
    request->contrast = REAL(contrast);
    Rsdeseq2ContrastZero zero = {0};
    if (context != R_NilValue) {
        rsdeseq2_validate_contrast_zero(context, genes, coefficients, &zero);
        request->zero = &zero;
    }
    SEXP result = PROTECT(rsdeseq2_allocate_contrast_output(genes, request));
    int status = rsdeseq2_r_contrast_batch((size_t)genes, (size_t)coefficients, REAL(beta),
        se != R_NilValue, se == R_NilValue ? NULL : REAL(se), covariance != R_NilValue,
        covariance == R_NilValue ? NULL : REAL(covariance), 1, request, NA_REAL);
    if (status) Rf_error("rsdeseq2 Rust backend: %s", rsdeseq2_r_last_error());
    UNPROTECT(1);
    return result;
}

SEXP rsdeseq2_project_contrast_zero(SEXP beta, SEXP se, SEXP covariance, SEXP contrast, SEXP context) {
    Rsdeseq2ContrastRequest request = {0};
    return contrast_inference(beta, se, covariance, contrast, &request, context);
}

SEXP rsdeseq2_project_contrast(SEXP beta, SEXP se, SEXP covariance, SEXP contrast) {
    return rsdeseq2_project_contrast_zero(beta, se, covariance, contrast, R_NilValue);
}

SEXP rsdeseq2_wald_contrast_zero(SEXP beta, SEXP se, SEXP covariance, SEXP contrast,
                          SEXP threshold, SEXP alternative, SEXP use_t, SEXP df, SEXP context) {
    rsdeseq2_validate_inference_matrices(beta, se, covariance);
    Rsdeseq2ContrastRequest request = {0};
    rsdeseq2_validate_wald_transport(threshold, alternative, use_t, df, matrix_rows(beta), &request);
    return contrast_inference(beta, se, covariance, contrast, &request, context);
}

SEXP rsdeseq2_wald_contrast(SEXP beta, SEXP se, SEXP covariance, SEXP contrast,
                          SEXP threshold, SEXP alternative, SEXP use_t, SEXP df) {
    return rsdeseq2_wald_contrast_zero(beta, se, covariance, contrast, threshold, alternative, use_t, df, R_NilValue);
}
