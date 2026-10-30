#ifndef RSDESEQ2_VALIDATION_H
#define RSDESEQ2_VALIDATION_H
#include <string.h>

/* Validate allocated R storage before constructing any native slice. */
static R_xlen_t checked_product(R_xlen_t a, R_xlen_t b) {
    if (a < 0 || b < 0 || (b && (a > R_XLEN_T_MAX / b ||
        (uint64_t)a > SIZE_MAX / sizeof(double) / (uint64_t)b ||
        (uint64_t)a > PTRDIFF_MAX / sizeof(double) / (uint64_t)b)))
        Rf_error("native matrix dimensions overflow");
    return a * b;
}
static void validate_vector(SEXP x, int type, R_xlen_t length, const char *name) {
    if (TYPEOF(x) != type || XLENGTH(x) != length)
        Rf_error("%s must have the required type and length (%lld)", name, (long long)length);
}
static void validate_optional_vector(SEXP x, int type, R_xlen_t length, const char *name) {
    if (x != R_NilValue) validate_vector(x, type, length, name);
}
static void validate_matrix(SEXP x, const char *name) {
    SEXP dims = Rf_getAttrib(x, R_DimSymbol);
    if (TYPEOF(x) != REALSXP || TYPEOF(dims) != INTSXP || XLENGTH(dims) != 2 ||
        INTEGER(dims)[0] < 1 || INTEGER(dims)[1] < 1)
        Rf_error("%s must be a nonempty double matrix", name);
    validate_vector(x, REALSXP, checked_product(INTEGER(dims)[0], INTEGER(dims)[1]), name);
}
static R_xlen_t matrix_rows(SEXP x) { return INTEGER(Rf_getAttrib(x, R_DimSymbol))[0]; }
static R_xlen_t matrix_cols(SEXP x) { return INTEGER(Rf_getAttrib(x, R_DimSymbol))[1]; }
static void validate_design(SEXP x, R_xlen_t rows, const char *name) {
    validate_matrix(x, name);
    if (matrix_rows(x) != rows) Rf_error("%s has incompatible rows", name);
}
static void validate_normalization(SEXP counts, SEXP sf, SEXP nf, SEXP weights) {
    validate_optional_vector(sf, REALSXP, matrix_cols(counts), "size_factors");
    if (nf != R_NilValue) {
        validate_matrix(nf, "normalization_factors");
        if (matrix_rows(nf) != matrix_rows(counts) || matrix_cols(nf) != matrix_cols(counts))
            Rf_error("normalization_factors dimensions must match counts");
    }
    if (weights != R_NilValue) {
        validate_matrix(weights, "weights");
        if (matrix_rows(weights) != matrix_rows(counts) || matrix_cols(weights) != matrix_cols(counts))
            Rf_error("weights dimensions must match counts");
    }
}
static int validate_integer_scalar(SEXP x, double minimum, double maximum, const char *name) {
    if ((TYPEOF(x) != INTSXP && TYPEOF(x) != REALSXP) || XLENGTH(x) != 1)
        Rf_error("%s must be an integer scalar", name);
    double value = Rf_asReal(x);
    if (!R_FINITE(value) || value < minimum || value > maximum || value != floor(value))
        Rf_error("%s is outside its supported integer range", name);
    return (int)value;
}
static void validate_real_scalar(SEXP x, const char *name) {
    if ((TYPEOF(x) != REALSXP && TYPEOF(x) != INTSXP) || XLENGTH(x) != 1)
        Rf_error("%s must be a numeric scalar", name);
}
static void validate_boolean(SEXP x, const char *name) {
    if (TYPEOF(x) != LGLSXP || XLENGTH(x) != 1 || LOGICAL(x)[0] == NA_LOGICAL)
        Rf_error("%s must be TRUE or FALSE", name);
}
/* Metadata lookup only; no numerical processing. */
static SEXP named_list_field(SEXP x, const char *name, const char *context) {
    SEXP names = Rf_getAttrib(x, R_NamesSymbol);
    if (TYPEOF(x) != VECSXP || TYPEOF(names) != STRSXP || XLENGTH(names) != XLENGTH(x))
        Rf_error("%s must be a named list", context);
    SEXP value = R_NilValue;
    int found = 0;
    for (R_xlen_t i = 0; i < XLENGTH(x); ++i) {
        if (STRING_ELT(names, i) == NA_STRING) Rf_error("%s names cannot be missing", context);
        if (!strcmp(CHAR(STRING_ELT(names, i)), name)) {
            if (found) Rf_error("duplicate %s field: %s", context, name);
            value = VECTOR_ELT(x, i); found = 1;
        }
    }
    return value;
}
#endif
