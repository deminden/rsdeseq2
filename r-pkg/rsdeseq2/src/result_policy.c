#include <R.h>
#include <Rinternals.h>
#include <math.h>
#include <stdint.h>
#include <stddef.h>
#include <limits.h>
#include <string.h>
#include "validation.h"
#include "result_policy.h"

extern const char *rsdeseq2_r_last_error(void);

static const int32_t *row_flags(SEXP values, R_xlen_t genes, const char *name) {
    if (values == R_NilValue) return NULL;
    validate_vector(values, LGLSXP, genes, name);
    return (const int32_t *)LOGICAL(values);
}

static const char *origin_string(SEXP values, R_xlen_t index) {
    if (TYPEOF(values) != STRSXP || index >= XLENGTH(values) || STRING_ELT(values,index) == NA_STRING)
        Rf_error("design origin identifiers must be nonmissing strings");
    const char *value = Rf_translateCharUTF8(STRING_ELT(values,index));
    if (!*value) Rf_error("design origin identifiers must be nonempty strings");
    return value;
}

static void result_origin(SEXP origin, Rsdeseq2ResultContext *native) {
    native->origin_kind = 2; /* compatibility for old direct-call contexts */
    if (origin == R_NilValue) return;
    SEXP kind = named_list_field(origin,"kind","design origin");
    validate_vector(kind,STRSXP,1,"design origin kind");
    const char *name = origin_string(kind,0);
    if (!strcmp(name,"matrix")) native->origin_kind = 0;
    else if (!strcmp(name,"formula")) native->origin_kind = 1;
    else if (strcmp(name,"legacy")) Rf_error("unsupported design origin kind");
    SEXP variables = named_list_field(origin,"variables","design origin");
    if (variables == R_NilValue && native->origin_kind != 1) return;
    validate_vector(variables,STRSXP,XLENGTH(variables),"design origin variables");
    if (native->origin_kind != 1 && XLENGTH(variables)) Rf_error("only formula origin can carry variables");
    checked_product(XLENGTH(variables),sizeof(const char *));
    const char **names = XLENGTH(variables) ? (const char **)R_alloc((size_t)XLENGTH(variables),sizeof(*names)) : NULL;
    for (R_xlen_t i=0;i<XLENGTH(variables);++i) {
        names[i]=origin_string(variables,i);
        for (R_xlen_t j=0;j<i;++j) if (!strcmp(names[i],names[j]))
            Rf_error("duplicate design origin variable");
    }
    native->variables_len=(size_t)XLENGTH(variables);native->variables=names;
}

SEXP rsdeseq2_present_results(SEXP columns, SEXP context, SEXP cutoff) {
    const char *names[] = {"log2FoldChange", "lfcSE", "stat", "pvalue", "cooksOutlier", ""};
    SEXP source[4];
    for (int i = 0; i < 4; ++i) source[i] = named_list_field(columns, names[i], "result columns");
    R_xlen_t genes = XLENGTH(source[0]);
    for (int i = 0; i < 4; ++i) validate_vector(source[i], REALSXP, genes, names[i]);
    checked_product(genes, 4);
    const int active = cutoff != R_NilValue;
    if (active) validate_real_scalar(cutoff, "cooksCutoff");
    Rsdeseq2ResultContext native = {0};
    native.refitted = row_flags(named_list_field(context, "refitRows", "result context"), genes, "refitRows");
    native.new_zero = row_flags(named_list_field(context, "newAllZeroRows", "result context"), genes, "newAllZeroRows");
    result_origin(named_list_field(context,"origin","result context"),&native);
    if (active) {
        SEXP maximum = named_list_field(context, "maxCooks", "result context");
        if (maximum == R_NilValue) Rf_error("The fit lacks Cook's context; refit");
        validate_vector(maximum, REALSXP, genes, "maxCooks"); native.max_cooks = REAL(maximum);
        SEXP design = named_list_field(context, "design", "result context");
        if (design != R_NilValue) {
            validate_matrix(design, "Cook's design");
            native.samples = (size_t)matrix_rows(design); native.coefficients = (size_t)matrix_cols(design);
            native.design = REAL(design);
        }
        SEXP counts = named_list_field(context, "counts", "result context");
        SEXP cooks = named_list_field(context, "cooks", "result context");
        if (counts != R_NilValue) {
            validate_matrix(counts, "Cook's counts");
            if (matrix_rows(counts) != genes || (native.samples && matrix_cols(counts) != (R_xlen_t)native.samples))
                Rf_error("Cook's counts dimensions must match result rows and samples");
            native.samples = (size_t)matrix_cols(counts); native.counts = REAL(counts);
        }
        if (cooks != R_NilValue) {
            validate_matrix(cooks, "Cook's distances");
            if (matrix_rows(cooks) != genes || (native.samples && matrix_cols(cooks) != (R_xlen_t)native.samples))
                Rf_error("Cook's distances dimensions must match result rows and samples");
            native.samples = (size_t)matrix_cols(cooks); native.cooks = REAL(cooks);
        }
        SEXP factors = named_list_field(context, "factors", "result context");
        if (factors != R_NilValue) {
            if (TYPEOF(factors) != VECSXP) Rf_error("Cook's factors must be a list");
            R_xlen_t count = XLENGTH(factors); checked_product(count, sizeof(Rsdeseq2FactorCodes));
            SEXP factor_names = Rf_getAttrib(factors,R_NamesSymbol);
            if (native.origin_kind == 1 && count && (TYPEOF(factor_names)!=STRSXP || XLENGTH(factor_names)!=count))
                Rf_error("formula Cook's factors require names");
            Rsdeseq2FactorCodes *codes = count ? (Rsdeseq2FactorCodes *)R_alloc((size_t)count, sizeof(*codes)) : NULL;
            for (R_xlen_t i = 0; i < count; ++i) {
                SEXP factor = VECTOR_ELT(factors, i);
                codes[i].name = TYPEOF(factor_names)==STRSXP ? origin_string(factor_names,i) : NULL;
                SEXP values = named_list_field(factor, "codes", "Cook's factor");
                validate_vector(values, INTSXP, (R_xlen_t)native.samples, "Cook's factor codes");
                codes[i].codes = (const int32_t *)INTEGER(values);
                codes[i].levels = (size_t)validate_integer_scalar(named_list_field(factor, "levels", "Cook's factor"), 0, INT_MAX, "Cook's factor levels");
            }
            native.factors_len = (size_t)count; native.factors = codes;
        }
    }
    uint8_t replacement_storage = 0;
    int status = rsdeseq2_r_result_replacement_storage((size_t)genes, native.refitted, native.new_zero, &replacement_storage);
    if (status) Rf_error("rsdeseq2 Rust backend: %s", rsdeseq2_r_last_error());
    SEXP output = PROTECT(Rf_mkNamed(VECSXP, names));
    double *values[4] = {NULL, NULL, NULL, NULL};
    for (int i = 0; i < 4; ++i) {
        if (replacement_storage || (active && i == 3)) {
            SEXP copy = PROTECT(Rf_duplicate(source[i])); SET_VECTOR_ELT(output, i, copy);
            values[i] = REAL(copy); UNPROTECT(1);
        } else SET_VECTOR_ELT(output, i, source[i]);
    }
    SEXP outlier = PROTECT(Rf_allocVector(LGLSXP, genes)); SET_VECTOR_ELT(output, 4, outlier);
    uint8_t *flags = (uint8_t *)R_alloc((size_t)genes, sizeof(uint8_t));
    status = rsdeseq2_r_present_results((size_t)genes, values[0], values[1], values[2], values[3],
        (uint8_t)active, active ? Rf_asReal(cutoff) : 0.0, &native, flags);
    if (status) Rf_error("rsdeseq2 Rust backend: %s", rsdeseq2_r_last_error());
    for (R_xlen_t gene = 0; gene < genes; ++gene)
        LOGICAL(outlier)[gene] = flags[gene] == UINT8_MAX ? NA_LOGICAL : (flags[gene] ? TRUE : FALSE);
    UNPROTECT(2);
    return output;
}
