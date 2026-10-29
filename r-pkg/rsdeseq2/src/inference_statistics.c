#include <R.h>
#include <Rinternals.h>
#include <math.h>
#include <stdint.h>
#include <stddef.h>
#include <limits.h>
#include "validation.h"
#include "inference_statistics.h"

extern const char *rsdeseq2_r_last_error(void);
static void check_statistics_status(int status) {
    if (status) Rf_error("rsdeseq2 Rust backend: %s", rsdeseq2_r_last_error());
}
static void validate_statistic_df(SEXP df, R_xlen_t genes) {
    validate_vector(df, REALSXP, XLENGTH(df), "degrees of freedom");
    if (XLENGTH(df) != 1 && XLENGTH(df) != genes)
        Rf_error("degrees of freedom must be scalar or one per gene");
    checked_product(XLENGTH(df), 1);
}
SEXP rsdeseq2_wald_statistics(SEXP beta, SEXP se, SEXP threshold, SEXP alternative, SEXP use_t, SEXP df) {
    validate_vector(beta, REALSXP, XLENGTH(beta), "beta");
    R_xlen_t genes = checked_product(XLENGTH(beta), 1);
    validate_vector(se, REALSXP, genes, "betaSE");
    validate_real_scalar(threshold, "lfcThreshold");
    int alt = validate_integer_scalar(alternative, 0, 5, "alternative");
    validate_boolean(use_t, "useT");
    if (LOGICAL(use_t)[0]) validate_statistic_df(df, genes);
    else if (df != R_NilValue) Rf_error("degrees of freedom require useT=TRUE");
    SEXP stat = PROTECT(Rf_allocVector(REALSXP, genes));
    SEXP pvalue = PROTECT(Rf_allocVector(REALSXP, genes));
    check_statistics_status(rsdeseq2_r_wald_statistics((size_t)genes, REAL(beta), REAL(se),
        Rf_asReal(threshold), alt, (uint8_t)LOGICAL(use_t)[0], df == R_NilValue ? 0 : (size_t)XLENGTH(df),
        df == R_NilValue ? NULL : REAL(df), NA_REAL, REAL(stat), REAL(pvalue)));
    const char *names[] = {"stat", "pvalue", ""};
    SEXP output = PROTECT(Rf_mkNamed(VECSXP, names));
    SET_VECTOR_ELT(output, 0, stat); SET_VECTOR_ELT(output, 1, pvalue);
    UNPROTECT(3);
    return output;
}
SEXP rsdeseq2_lrt_pvalues(SEXP statistics, SEXP df) {
    validate_vector(statistics, REALSXP, XLENGTH(statistics), "LRT statistics");
    R_xlen_t genes = checked_product(XLENGTH(statistics), 1);
    validate_statistic_df(df, genes);
    SEXP output = PROTECT(Rf_allocVector(REALSXP, genes));
    check_statistics_status(rsdeseq2_r_lrt_pvalues((size_t)genes, REAL(statistics),
        (size_t)XLENGTH(df), REAL(df), NA_REAL, REAL(output)));
    Rf_setAttrib(output, R_NamesSymbol, Rf_getAttrib(statistics, R_NamesSymbol));
    UNPROTECT(1);
    return output;
}
SEXP rsdeseq2_default_cooks_cutoff(SEXP samples, SEXP coefficients) {
    int m = validate_integer_scalar(samples, 0, INT_MAX, "samples");
    int p = validate_integer_scalar(coefficients, 0, INT_MAX, "coefficients");
    double cutoff = 0;
    uint8_t present = 0;
    check_statistics_status(rsdeseq2_r_default_cooks_cutoff((size_t)m, (size_t)p, &cutoff, &present));
    return present ? Rf_ScalarReal(cutoff) : R_NilValue;
}
