#include <R.h>
#include <Rinternals.h>
#include <math.h>
#include <stdint.h>
#include <stddef.h>
#include <limits.h>
#include <string.h>
#include "validation.h"
#include "contrast_zero.h"

static SEXP context_field(SEXP context, const char *name) {
    return named_list_field(context, name, "contrast zero context");
}

static const char *level_string(SEXP values, R_xlen_t index, const char *name) {
    if (TYPEOF(values) != STRSXP || index >= XLENGTH(values) || STRING_ELT(values, index) == NA_STRING)
        Rf_error("%s must contain nonmissing strings", name);
    const char *value = Rf_translateCharUTF8(STRING_ELT(values, index));
    if (!*value) Rf_error("%s must contain nonempty strings", name);
    return value;
}

void rsdeseq2_validate_contrast_zero(SEXP context, R_xlen_t genes, R_xlen_t coefficients,
                                    Rsdeseq2ContrastZero *output) {
    SEXP counts = context_field(context, "counts");
    SEXP rule = context_field(context, "rule");
    validate_matrix(counts, "contrast zero counts");
    if (matrix_rows(counts) != genes) Rf_error("contrast zero counts dimensions must match beta rows");
    output->samples = (size_t)matrix_cols(counts);
    output->counts = REAL(counts);
    output->rule = validate_integer_scalar(rule, 1, 2, "contrast zero rule");
    if (output->rule == 1) {
        SEXP design = context_field(context, "design");
        if (design != R_NilValue) {
            validate_matrix(design, "contrast zero design");
            if (matrix_rows(design) != (R_xlen_t)output->samples || matrix_cols(design) != coefficients)
                Rf_error("contrast zero design dimensions must match samples and coefficients");
            output->design = REAL(design);
            output->has_design = 1;
        }
    } else {
        SEXP levels = context_field(context, "sampleLevels");
        SEXP numerator = context_field(context, "numerator");
        SEXP denominator = context_field(context, "denominator");
        validate_vector(levels, STRSXP, (R_xlen_t)output->samples, "sampleLevels");
        validate_vector(numerator, STRSXP, 1, "numerator");
        validate_vector(denominator, STRSXP, 1, "denominator");
        const char **strings = (const char **)R_alloc(output->samples, sizeof(const char *));
        for (R_xlen_t i = 0; i < (R_xlen_t)output->samples; ++i)
            strings[i] = level_string(levels, i, "sampleLevels");
        output->sample_levels = strings;
        output->numerator = level_string(numerator, 0, "numerator");
        output->denominator = level_string(denominator, 0, "denominator");
    }
}
