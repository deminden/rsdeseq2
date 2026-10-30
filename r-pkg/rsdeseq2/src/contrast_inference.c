#include <R.h>
#include <Rinternals.h>
#include <math.h>
#include <stdint.h>
#include <stddef.h>
#include <limits.h>
#include "validation.h"
#include "contrast_inference.h"

extern const char *rsdeseq2_r_last_error(void);

static void validate_estimates(SEXP beta, SEXP se, SEXP covariance, SEXP contrast) {
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
    validate_vector(contrast, REALSXP, coefficients, "contrast");
}

static SEXP contrast_inference(SEXP beta, SEXP se, SEXP covariance, SEXP contrast,
                              int wald, double threshold, int alternative, int use_t, SEXP df, SEXP context) {
    R_xlen_t genes = matrix_rows(beta);
    Rsdeseq2ContrastZero zero = {0};
    if (context != R_NilValue)
        rsdeseq2_validate_contrast_zero(context, genes, matrix_cols(beta), &zero);
    SEXP lfc = PROTECT(Rf_allocVector(REALSXP, genes));
    SEXP lfc_se = PROTECT(Rf_allocVector(REALSXP, genes));
    SEXP stat = PROTECT(wald ? Rf_allocVector(REALSXP, genes) : R_NilValue);
    SEXP pvalue = PROTECT(wald ? Rf_allocVector(REALSXP, genes) : R_NilValue);
    int status = rsdeseq2_r_contrast_inference((size_t)genes, (size_t)matrix_cols(beta), REAL(beta),
        se != R_NilValue, se == R_NilValue ? NULL : REAL(se),
        covariance != R_NilValue, covariance == R_NilValue ? NULL : REAL(covariance),
        REAL(contrast), (uint8_t)wald, threshold, alternative, (uint8_t)use_t,
        df == R_NilValue ? 0 : (size_t)XLENGTH(df), df == R_NilValue ? NULL : REAL(df),
        NA_REAL, REAL(lfc), REAL(lfc_se), wald ? REAL(stat) : NULL, wald ? REAL(pvalue) : NULL,
        context == R_NilValue ? NULL : &zero);
    if (status) Rf_error("rsdeseq2 Rust backend: %s", rsdeseq2_r_last_error());
    const char *wald_names[] = {"log2FoldChange", "lfcSE", "stat", "pvalue", ""};
    const char *projection_names[] = {"log2FoldChange", "lfcSE", ""};
    SEXP output = PROTECT(Rf_mkNamed(VECSXP, wald ? wald_names : projection_names));
    SET_VECTOR_ELT(output, 0, lfc); SET_VECTOR_ELT(output, 1, lfc_se);
    if (wald) { SET_VECTOR_ELT(output, 2, stat); SET_VECTOR_ELT(output, 3, pvalue); }
    UNPROTECT(5);
    return output;
}

SEXP rsdeseq2_project_contrast_zero(SEXP beta, SEXP se, SEXP covariance, SEXP contrast, SEXP context) {
    validate_estimates(beta, se, covariance, contrast);
    return contrast_inference(beta, se, covariance, contrast, 0, 0, 0, 0, R_NilValue, context);
}

SEXP rsdeseq2_project_contrast(SEXP beta, SEXP se, SEXP covariance, SEXP contrast) {
    return rsdeseq2_project_contrast_zero(beta, se, covariance, contrast, R_NilValue);
}

SEXP rsdeseq2_wald_contrast_zero(SEXP beta, SEXP se, SEXP covariance, SEXP contrast,
                          SEXP threshold, SEXP alternative, SEXP use_t, SEXP df, SEXP context) {
    validate_estimates(beta, se, covariance, contrast);
    validate_real_scalar(threshold, "lfcThreshold");
    int alt = validate_integer_scalar(alternative, 0, 5, "alternative");
    validate_boolean(use_t, "useT");
    if (LOGICAL(use_t)[0]) {
        validate_vector(df, REALSXP, XLENGTH(df), "degrees of freedom");
        if (XLENGTH(df) != 1 && XLENGTH(df) != matrix_rows(beta))
            Rf_error("degrees of freedom must be scalar or one per gene");
        checked_product(XLENGTH(df), 1);
    } else if (df != R_NilValue) Rf_error("degrees of freedom require useT=TRUE");
    return contrast_inference(beta, se, covariance, contrast, 1,
        Rf_asReal(threshold), alt, LOGICAL(use_t)[0], df, context);
}

SEXP rsdeseq2_wald_contrast(SEXP beta, SEXP se, SEXP covariance, SEXP contrast,
                          SEXP threshold, SEXP alternative, SEXP use_t, SEXP df) {
    return rsdeseq2_wald_contrast_zero(beta, se, covariance, contrast, threshold, alternative, use_t, df, R_NilValue);
}
