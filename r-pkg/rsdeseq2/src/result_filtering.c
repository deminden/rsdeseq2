#include <R.h>
#include <Rinternals.h>
#include <math.h>
#include <stdint.h>
#include <stddef.h>
#include <limits.h>
#include "validation.h"
#include "result_filtering.h"

extern const char *rsdeseq2_r_last_error(void);
static void check_filter_status(int status) {
    if (status) Rf_error("rsdeseq2 Rust backend: %s", rsdeseq2_r_last_error());
}
static void validate_filter_pvalues(SEXP values) {
    validate_vector(values, REALSXP, XLENGTH(values), "p-values");
    checked_product(XLENGTH(values), 1);
    for (R_xlen_t i = 0; i < XLENGTH(values); i++) {
        double value = REAL(values)[i];
        if (!ISNAN(value) && (!R_FINITE(value) || value < 0 || value > 1))
            Rf_error("p-values must be in [0, 1] or missing");
    }
}
SEXP rsdeseq2_bh_adjust(SEXP values) {
    validate_filter_pvalues(values);
    SEXP output = PROTECT(Rf_allocVector(REALSXP, XLENGTH(values)));
    check_filter_status(rsdeseq2_r_bh_adjust((size_t)XLENGTH(values), REAL(values), NA_REAL, REAL(output)));
    Rf_setAttrib(output, R_NamesSymbol, Rf_getAttrib(values, R_NamesSymbol));
    UNPROTECT(1);
    return output;
}
SEXP rsdeseq2_independent_filtering(SEXP filter, SEXP values, SEXP theta, SEXP alpha, SEXP enabled) {
    validate_filter_pvalues(values);
    validate_vector(filter, REALSXP, XLENGTH(values), "baseMean");
    validate_real_scalar(alpha, "alpha");
    validate_boolean(enabled, "enabled");
    if (!R_FINITE(Rf_asReal(alpha)) || Rf_asReal(alpha) <= 0 || Rf_asReal(alpha) >= 1)
        Rf_error("alpha must be in (0, 1)");
    int active = LOGICAL(enabled)[0];
    if (active && theta != R_NilValue) {
        validate_vector(theta, REALSXP, XLENGTH(theta), "theta");
        if (XLENGTH(theta) <= 1 || XLENGTH(theta) > INT_MAX) Rf_error("theta length is outside the supported integer range");
    }
    R_xlen_t genes = XLENGTH(values), k = active ? (theta == R_NilValue ? 50 : XLENGTH(theta)) : 0;
    checked_product(k, 1);
    SEXP padj = PROTECT(Rf_allocVector(REALSXP, genes));
    SEXP flags = PROTECT(Rf_allocVector(LGLSXP, genes));
    SEXP theta_output = PROTECT(Rf_allocVector(REALSXP, k));
    SEXP rejections = PROTECT(Rf_allocVector(INTSXP, k));
    SEXP smooth = PROTECT(Rf_allocVector(REALSXP, k));
    SEXP selected = PROTECT(Rf_allocVector(INTSXP, 1));
    SEXP scalars = PROTECT(Rf_allocVector(REALSXP, 2));
    uint8_t *native_flags = (uint8_t *)R_alloc((size_t)genes, sizeof(uint8_t));
    check_filter_status(rsdeseq2_r_independent_filtering((size_t)genes, REAL(filter), REAL(values),
        active && theta != R_NilValue ? (size_t)XLENGTH(theta) : 0,
        active && theta != R_NilValue ? REAL(theta) : NULL, (uint8_t)active,
        Rf_asReal(alpha), NA_REAL, REAL(padj), native_flags, REAL(theta_output),
        INTEGER(rejections), REAL(smooth), INTEGER(selected), REAL(scalars)));
    for (R_xlen_t i = 0; i < genes; i++) LOGICAL(flags)[i] = native_flags[i] == UINT8_MAX ? NA_LOGICAL : native_flags[i];
    if (INTEGER(selected)[0] < 0) INTEGER(selected)[0] = NA_INTEGER;
    const char *names[] = {"padj", "filtered", "theta", "numRej", "selected", "scalars", "lo.fit", ""};
    SEXP output = PROTECT(Rf_mkNamed(VECSXP, names));
    SET_VECTOR_ELT(output, 0, padj); SET_VECTOR_ELT(output, 1, flags);
    SET_VECTOR_ELT(output, 2, theta_output); SET_VECTOR_ELT(output, 3, rejections);
    SET_VECTOR_ELT(output, 4, selected); SET_VECTOR_ELT(output, 5, scalars);
    if (active) {
        const char *curve_names[] = {"x", "y", ""};
        SEXP curve = PROTECT(Rf_mkNamed(VECSXP, curve_names));
        SET_VECTOR_ELT(curve, 0, theta_output); SET_VECTOR_ELT(curve, 1, smooth);
        SET_VECTOR_ELT(output, 6, curve);
        UNPROTECT(1);
    }
    UNPROTECT(8);
    return output;
}
