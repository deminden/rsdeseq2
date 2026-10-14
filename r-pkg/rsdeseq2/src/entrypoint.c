#include <R.h>
#include <R_ext/Rdynload.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <stddef.h>

extern int rsdeseq2_r_estimate_size_factors(size_t, size_t, const double *, int, const double *, const size_t *, size_t, double *);
extern int rsdeseq2_r_normalized_counts(size_t, size_t, const double *, const double *, const double *, double *);
extern int rsdeseq2_r_base_metadata(size_t, size_t, const double *, const double *, const double *, const double *, double *, double *, uint8_t *);
extern int rsdeseq2_r_fit(size_t, size_t, const double *, const double *, const double *, const double *, size_t, const double *, size_t, const double *, int, int, size_t, int, double *, double *, double *, double *, double *, double *, double *);
extern int rsdeseq2_r_lfc_shrink_normal(size_t, size_t, const double *, const double *, const double *, const double *, size_t, const double *, const double *, const double *, const double *, size_t, ptrdiff_t, double, double *, double *, double *, double *, double *);
extern int rsdeseq2_r_lfc_shrink_normal_contrast(size_t, size_t, const double *, const double *, const double *, const double *, size_t, const double *, size_t, const double *, const int32_t *, const int32_t *, const double *, const double *, const double *, const int32_t *, const double *, ptrdiff_t, ptrdiff_t, double, double *, double *, double *, double *, double *);
extern int rsdeseq2_r_transform(size_t, size_t, const double *, const double *, const double *, const double *, size_t, const double *, int, int, uint8_t, size_t, size_t, const double *, double, int, const double *, const double *, const double *, double *, double *, double *, double *);
extern int rsdeseq2_r_cooks_mask(size_t, const double *, const double *, double, const double *, size_t, const double *, uint8_t, double *, uint8_t *);
extern const char *rsdeseq2_r_last_error(void);
static void check_rust_backend(int status) {
    if (status != 0) Rf_error("rsdeseq2 Rust backend: %s", rsdeseq2_r_last_error());
}

SEXP rsdeseq2_placeholder(void) {
    Rf_error("rsdeseq2 native bridge is not implemented yet");
    return R_NilValue;
}

SEXP rsdeseq2_diagnostic_schema(void) {
    const char *names[] = {
        "betaConv",
        "fullBetaConv",
        "reducedBetaConv",
        "betaIter",
        "reducedBetaIter",
        "deviance",
        "maxCooks",
    };
    const R_xlen_t n_names = (R_xlen_t)(sizeof(names) / sizeof(names[0]));
    SEXP out = PROTECT(Rf_allocVector(STRSXP, n_names));
    for (R_xlen_t i = 0; i < n_names; ++i) {
        SET_STRING_ELT(out, i, Rf_mkChar(names[i]));
    }
    UNPROTECT(1);
    return out;
}

SEXP rsdeseq2_estimate_size_factors(SEXP counts, SEXP geo_means, SEXP control_genes, SEXP method) {
    SEXP dims = Rf_getAttrib(counts, R_DimSymbol);
    if (Rf_length(dims) != 2) Rf_error("counts must be a matrix");
    const size_t n_genes = (size_t)INTEGER(dims)[0], n_samples = (size_t)INTEGER(dims)[1];
    const size_t n_control = control_genes == R_NilValue ? n_genes : (size_t)Rf_xlength(control_genes);
    size_t *controls = (size_t *)R_alloc(n_control, sizeof(size_t));
    if (control_genes == R_NilValue) for (size_t i = 0; i < n_control; ++i) controls[i] = i;
    else for (size_t i = 0; i < n_control; ++i) controls[i] = (size_t)(INTEGER(control_genes)[i] - 1);
    SEXP out = PROTECT(Rf_allocVector(REALSXP, n_samples));
    check_rust_backend(rsdeseq2_r_estimate_size_factors(n_genes, n_samples, REAL(counts), Rf_asInteger(method),
        geo_means == R_NilValue ? NULL : REAL(geo_means), controls, n_control, REAL(out)));
    SEXP dimnames = Rf_getAttrib(counts, R_DimNamesSymbol);
    if (dimnames != R_NilValue && Rf_length(dimnames) >= 2 && VECTOR_ELT(dimnames, 1) != R_NilValue)
        Rf_setAttrib(out, R_NamesSymbol, VECTOR_ELT(dimnames, 1));
    UNPROTECT(1); return out;
}
SEXP rsdeseq2_normalized_counts(SEXP counts, SEXP size_factors, SEXP normalization_factors) {
    SEXP dims = Rf_getAttrib(counts, R_DimSymbol);
    if (Rf_length(dims) != 2) Rf_error("counts must be a matrix");
    const size_t n_genes = (size_t)INTEGER(dims)[0], n_samples = (size_t)INTEGER(dims)[1];
    SEXP out = PROTECT(Rf_allocMatrix(REALSXP, (int)n_genes, (int)n_samples));
    check_rust_backend(rsdeseq2_r_normalized_counts(n_genes, n_samples, REAL(counts),
        size_factors == R_NilValue ? NULL : REAL(size_factors),
        normalization_factors == R_NilValue ? NULL : REAL(normalization_factors), REAL(out)));
    SEXP dimnames = Rf_getAttrib(counts, R_DimNamesSymbol);
    if (dimnames != R_NilValue) Rf_setAttrib(out, R_DimNamesSymbol, dimnames);
    UNPROTECT(1); return out;
}
SEXP rsdeseq2_base_mean(SEXP counts, SEXP size_factors, SEXP normalization_factors) {
    SEXP dims = Rf_getAttrib(counts, R_DimSymbol);
    if (Rf_length(dims) != 2) Rf_error("counts must be a matrix");
    const size_t n_genes = (size_t)INTEGER(dims)[0], n_samples = (size_t)INTEGER(dims)[1];
    SEXP out = PROTECT(Rf_allocVector(REALSXP, n_genes));
    SEXP variance = PROTECT(Rf_allocVector(REALSXP, n_genes));
    SEXP zero = PROTECT(Rf_allocVector(RAWSXP, n_genes));
    check_rust_backend(rsdeseq2_r_base_metadata(n_genes, n_samples, REAL(counts),
        size_factors == R_NilValue ? NULL : REAL(size_factors),
        normalization_factors == R_NilValue ? NULL : REAL(normalization_factors), NULL,
        REAL(out), REAL(variance), RAW(zero)));
    SEXP dimnames = Rf_getAttrib(counts, R_DimNamesSymbol);
    if (dimnames != R_NilValue && Rf_length(dimnames) >= 1 && VECTOR_ELT(dimnames, 0) != R_NilValue)
        Rf_setAttrib(out, R_NamesSymbol, VECTOR_ELT(dimnames, 0));
    UNPROTECT(3); return out;
}
SEXP rsdeseq2_base_metadata(SEXP counts, SEXP size_factors, SEXP normalization_factors, SEXP weights) {
    SEXP dims = Rf_getAttrib(counts, R_DimSymbol);
    if (Rf_length(dims) != 2) Rf_error("counts must be a matrix");
    const size_t n_genes = (size_t)INTEGER(dims)[0], n_samples = (size_t)INTEGER(dims)[1];
    SEXP base_mean = PROTECT(Rf_allocVector(REALSXP, n_genes));
    SEXP base_var = PROTECT(Rf_allocVector(REALSXP, n_genes));
    SEXP all_zero = PROTECT(Rf_allocVector(LGLSXP, n_genes));
    uint8_t *zero = (uint8_t *)R_alloc(n_genes, sizeof(uint8_t));
    check_rust_backend(rsdeseq2_r_base_metadata(n_genes, n_samples, REAL(counts),
        size_factors == R_NilValue ? NULL : REAL(size_factors),
        normalization_factors == R_NilValue ? NULL : REAL(normalization_factors),
        weights == R_NilValue ? NULL : REAL(weights), REAL(base_mean), REAL(base_var), zero));
    for (size_t i = 0; i < n_genes; ++i) LOGICAL(all_zero)[i] = zero[i] ? TRUE : FALSE;
    SEXP out = PROTECT(Rf_allocVector(VECSXP, 3));
    SET_VECTOR_ELT(out, 0, base_mean); SET_VECTOR_ELT(out, 1, base_var); SET_VECTOR_ELT(out, 2, all_zero);
    SEXP names = PROTECT(Rf_allocVector(STRSXP, 3));
    SET_STRING_ELT(names, 0, Rf_mkChar("baseMean")); SET_STRING_ELT(names, 1, Rf_mkChar("baseVar")); SET_STRING_ELT(names, 2, Rf_mkChar("allZero"));
    Rf_setAttrib(out, R_NamesSymbol, names);
    SEXP dimnames = Rf_getAttrib(counts, R_DimNamesSymbol);
    SEXP row_names = dimnames == R_NilValue ? R_NilValue : VECTOR_ELT(dimnames, 0);
    if (row_names == R_NilValue) {
        row_names = PROTECT(Rf_allocVector(INTSXP, n_genes));
        for (size_t i = 0; i < n_genes; ++i) INTEGER(row_names)[i] = (int)i + 1;
        Rf_setAttrib(out, R_RowNamesSymbol, row_names); UNPROTECT(1);
    } else Rf_setAttrib(out, R_RowNamesSymbol, row_names);
    SEXP classes = PROTECT(Rf_allocVector(STRSXP, 1)); SET_STRING_ELT(classes, 0, Rf_mkChar("data.frame")); Rf_setAttrib(out, R_ClassSymbol, classes);
    UNPROTECT(6); return out;
}

SEXP rsdeseq2_fit(SEXP counts, SEXP design, SEXP reduced_design, SEXP size_factors_input, SEXP normalization_factors, SEXP observation_weights, SEXP fit_type, SEXP test_type, SEXP threads, SEXP execution_mode) {
    SEXP count_dims = Rf_getAttrib(counts, R_DimSymbol);
    SEXP design_dims = Rf_getAttrib(design, R_DimSymbol);
    if (Rf_length(count_dims) != 2 || Rf_length(design_dims) != 2)
        Rf_error("counts and design must be matrices");
    const size_t n_genes = (size_t)INTEGER(count_dims)[0];
    const size_t n_samples = (size_t)INTEGER(count_dims)[1];
    if ((size_t)INTEGER(design_dims)[0] != n_samples)
        Rf_error("design must contain one row per count sample");
    size_t reduced_coefficients = 0;
    const double *reduced_values = NULL;
    if (reduced_design != R_NilValue) {
        SEXP reduced_dims = Rf_getAttrib(reduced_design, R_DimSymbol);
        if (Rf_length(reduced_dims) != 2 || (size_t)INTEGER(reduced_dims)[0] != n_samples)
            Rf_error("reduced design must contain one row per count sample");
        reduced_coefficients = (size_t)INTEGER(reduced_dims)[1];
        reduced_values = REAL(reduced_design);
    }
    const char *columns[] = {"baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "dispersion", "converged", "maxCooks", "cooksOutlier", "filtered"};
    const size_t n_columns = sizeof(columns) / sizeof(columns[0]);
    const size_t n_coefficients = (size_t)INTEGER(design_dims)[1];
    SEXP values = PROTECT(Rf_allocMatrix(REALSXP, (int)n_genes, (int)n_columns));
    SEXP beta = PROTECT(Rf_allocMatrix(REALSXP, (int)n_genes, (int)n_coefficients));
    SEXP beta_se = PROTECT(Rf_allocMatrix(REALSXP, (int)n_genes, (int)n_coefficients));
    SEXP beta_covariance = PROTECT(Rf_allocVector(REALSXP, n_genes * n_coefficients * n_coefficients));
    SEXP covariance_dims = PROTECT(Rf_allocVector(INTSXP, 3));
    INTEGER(covariance_dims)[0] = (int)n_genes;
    INTEGER(covariance_dims)[1] = (int)n_coefficients;
    INTEGER(covariance_dims)[2] = (int)n_coefficients;
    Rf_setAttrib(beta_covariance, R_DimSymbol, covariance_dims);
    SEXP size_factors_out = PROTECT(Rf_allocVector(REALSXP, n_samples));
    SEXP dispersions_out = PROTECT(Rf_allocVector(REALSXP, n_genes));
    SEXP disp_fit_out = PROTECT(Rf_allocVector(REALSXP, n_genes));
    check_rust_backend(rsdeseq2_r_fit(n_genes, n_samples, REAL(counts),
        size_factors_input == R_NilValue ? NULL : REAL(size_factors_input),
        normalization_factors == R_NilValue ? NULL : REAL(normalization_factors),
        observation_weights == R_NilValue ? NULL : REAL(observation_weights),
        n_coefficients, REAL(design), reduced_coefficients, reduced_values,
        Rf_asInteger(fit_type), Rf_asInteger(test_type), (size_t)Rf_asInteger(threads), Rf_asInteger(execution_mode), REAL(values),
        REAL(beta), REAL(beta_se), REAL(beta_covariance), REAL(size_factors_out),
        REAL(dispersions_out), REAL(disp_fit_out)));
    SEXP out = PROTECT(Rf_allocVector(VECSXP, n_columns));
    for (size_t column = 0; column < n_columns; ++column) {
        const int logical_column = column == 7 || column == 9 || column == 10;
        SEXP value = PROTECT(Rf_allocVector(logical_column ? LGLSXP : REALSXP, n_genes));
        if (logical_column) {
            for (size_t gene = 0; gene < n_genes; ++gene) {
                const double item = REAL(values)[gene + column * n_genes];
                LOGICAL(value)[gene] = ISNAN(item) ? NA_LOGICAL : (item != 0.0 ? TRUE : FALSE);
            }
        } else {
            memcpy(REAL(value), REAL(values) + column * n_genes, n_genes * sizeof(double));
        }
        SET_VECTOR_ELT(out, column, value);
        UNPROTECT(1);
    }
    SEXP names = PROTECT(Rf_allocVector(STRSXP, n_columns));
    for (size_t column = 0; column < n_columns; ++column) SET_STRING_ELT(names, column, Rf_mkChar(columns[column]));
    Rf_setAttrib(out, R_NamesSymbol, names);
    SEXP row_names = PROTECT(Rf_allocVector(STRSXP, n_genes));
    SEXP dimnames = Rf_getAttrib(counts, R_DimNamesSymbol);
    SEXP gene_names = dimnames == R_NilValue ? R_NilValue : VECTOR_ELT(dimnames, 0);
    if (gene_names == R_NilValue) {
        for (size_t gene = 0; gene < n_genes; ++gene) SET_STRING_ELT(row_names, gene, Rf_mkChar(""));
        Rf_setAttrib(out, R_RowNamesSymbol, R_NilValue);
    } else {
        for (size_t gene = 0; gene < n_genes; ++gene) SET_STRING_ELT(row_names, gene, STRING_ELT(gene_names, gene));
        Rf_setAttrib(out, R_RowNamesSymbol, row_names);
    }
    SEXP classes = PROTECT(Rf_allocVector(STRSXP, 1));
    SET_STRING_ELT(classes, 0, Rf_mkChar("data.frame"));
    Rf_setAttrib(out, R_ClassSymbol, classes);
    SEXP beta_dimnames = PROTECT(Rf_allocVector(VECSXP, 2));
    SEXP count_dimnames = Rf_getAttrib(counts, R_DimNamesSymbol);
    SEXP design_dimnames = Rf_getAttrib(design, R_DimNamesSymbol);
    SEXP coefficient_names = design_dimnames == R_NilValue ? R_NilValue : VECTOR_ELT(design_dimnames, 1);
    SET_VECTOR_ELT(beta_dimnames, 0, gene_names);
    SET_VECTOR_ELT(beta_dimnames, 1, coefficient_names);
    Rf_setAttrib(beta, R_DimNamesSymbol, beta_dimnames);
    Rf_setAttrib(beta_se, R_DimNamesSymbol, beta_dimnames);
    SEXP covariance_names = PROTECT(Rf_allocVector(VECSXP, 3));
    SET_VECTOR_ELT(covariance_names, 0, gene_names);
    SET_VECTOR_ELT(covariance_names, 1, coefficient_names);
    SET_VECTOR_ELT(covariance_names, 2, coefficient_names);
    Rf_setAttrib(beta_covariance, R_DimNamesSymbol, covariance_names);
    if (count_dimnames != R_NilValue && Rf_length(count_dimnames) >= 2)
        Rf_setAttrib(size_factors_out, R_NamesSymbol, VECTOR_ELT(count_dimnames, 1));
    SEXP fit = PROTECT(Rf_allocVector(VECSXP, 7));
    SET_VECTOR_ELT(fit, 0, out);
    SET_VECTOR_ELT(fit, 1, beta);
    SET_VECTOR_ELT(fit, 2, beta_se);
    SET_VECTOR_ELT(fit, 3, beta_covariance);
    SET_VECTOR_ELT(fit, 4, size_factors_out);
    SET_VECTOR_ELT(fit, 5, dispersions_out);
    SET_VECTOR_ELT(fit, 6, disp_fit_out);
    SEXP fit_names = PROTECT(Rf_allocVector(STRSXP, 7));
    SET_STRING_ELT(fit_names, 0, Rf_mkChar("results"));
    SET_STRING_ELT(fit_names, 1, Rf_mkChar("beta"));
    SET_STRING_ELT(fit_names, 2, Rf_mkChar("betaSE"));
    SET_STRING_ELT(fit_names, 3, Rf_mkChar("betaCovariance"));
    SET_STRING_ELT(fit_names, 4, Rf_mkChar("sizeFactors"));
    SET_STRING_ELT(fit_names, 5, Rf_mkChar("dispersions"));
    SET_STRING_ELT(fit_names, 6, Rf_mkChar("dispFit"));
    Rf_setAttrib(fit, R_NamesSymbol, fit_names);
    UNPROTECT(16);
    return fit;
}

SEXP rsdeseq2_transform(SEXP counts, SEXP design, SEXP size_factors_input,
                        SEXP normalization_factors, SEXP observation_weights,
                        SEXP fit_type, SEXP transform_type, SEXP blind,
                        SEXP nsub, SEXP threads, SEXP frozen_intercept,
                        SEXP sample_prior_variance, SEXP frozen_trend_type,
                        SEXP frozen_trend_parameters, SEXP frozen_disp_fit,
                        SEXP frozen_local_curve) {
    SEXP count_dims = Rf_getAttrib(counts, R_DimSymbol);
    SEXP design_dims = Rf_getAttrib(design, R_DimSymbol);
    if (Rf_length(count_dims) != 2 || Rf_length(design_dims) != 2)
        Rf_error("counts and design must be matrices");
    const size_t n_genes = (size_t)INTEGER(count_dims)[0];
    const size_t n_samples = (size_t)INTEGER(count_dims)[1];
    if ((size_t)INTEGER(design_dims)[0] != n_samples)
        Rf_error("design must contain one row per count sample");

    SEXP transformed = PROTECT(Rf_allocMatrix(REALSXP, (int)n_genes, (int)n_samples));
    SEXP size_factors = PROTECT(Rf_allocVector(REALSXP, n_samples));
    SEXP intercept = PROTECT(Rf_allocVector(REALSXP, n_genes));
    SEXP prior_variance = PROTECT(Rf_allocVector(REALSXP, 1));
    check_rust_backend(rsdeseq2_r_transform(
        n_genes, n_samples, REAL(counts),
        size_factors_input == R_NilValue ? NULL : REAL(size_factors_input),
        normalization_factors == R_NilValue ? NULL : REAL(normalization_factors),
        observation_weights == R_NilValue ? NULL : REAL(observation_weights),
        (size_t)INTEGER(design_dims)[1], REAL(design),
        Rf_asInteger(fit_type), Rf_asInteger(transform_type),
        (uint8_t)(LOGICAL(blind)[0] == TRUE), (size_t)Rf_asInteger(nsub),
        (size_t)Rf_asInteger(threads),
        frozen_intercept == R_NilValue ? NULL : REAL(frozen_intercept),
        sample_prior_variance == R_NilValue ? NAN : Rf_asReal(sample_prior_variance),
        Rf_asInteger(frozen_trend_type),
        frozen_trend_parameters == R_NilValue ? NULL : REAL(frozen_trend_parameters),
        frozen_disp_fit == R_NilValue ? NULL : REAL(frozen_disp_fit),
        frozen_local_curve == R_NilValue ? NULL : REAL(frozen_local_curve),
        REAL(transformed), REAL(size_factors),
        REAL(intercept), REAL(prior_variance)));
    SEXP dimnames = Rf_getAttrib(counts, R_DimNamesSymbol);
    if (dimnames != R_NilValue) Rf_setAttrib(transformed, R_DimNamesSymbol, dimnames);
    if (dimnames != R_NilValue && VECTOR_ELT(dimnames, 1) != R_NilValue)
        Rf_setAttrib(size_factors, R_NamesSymbol, VECTOR_ELT(dimnames, 1));

    SEXP out = PROTECT(Rf_allocVector(VECSXP, 4));
    SET_VECTOR_ELT(out, 0, transformed);
    SET_VECTOR_ELT(out, 1, size_factors);
    SET_VECTOR_ELT(out, 2, Rf_asInteger(transform_type) == 1 ? intercept : R_NilValue);
    SET_VECTOR_ELT(out, 3, Rf_asInteger(transform_type) == 1 ? prior_variance : R_NilValue);
    SEXP names = PROTECT(Rf_allocVector(STRSXP, 4));
    SET_STRING_ELT(names, 0, Rf_mkChar("transformed"));
    SET_STRING_ELT(names, 1, Rf_mkChar("sizeFactors"));
    SET_STRING_ELT(names, 2, Rf_mkChar("intercept"));
    SET_STRING_ELT(names, 3, Rf_mkChar("samplePriorVariance"));
    Rf_setAttrib(out, R_NamesSymbol, names);
    UNPROTECT(6);
    return out;
}

SEXP rsdeseq2_lfc_shrink_normal(SEXP counts, SEXP design,
                                SEXP size_factors, SEXP normalization_factors,
                                SEXP observation_weights, SEXP dispersions,
                                SEXP base_mean, SEXP disp_fit, SEXP coefficient,
                                SEXP intercept_index, SEXP lfc_threshold) {
    SEXP count_dims = Rf_getAttrib(counts, R_DimSymbol);
    SEXP design_dims = Rf_getAttrib(design, R_DimSymbol);
    if (Rf_length(count_dims) != 2 || Rf_length(design_dims) != 2)
        Rf_error("counts and design must be matrices");
    const size_t n_genes = (size_t)INTEGER(count_dims)[0];
    const size_t n_samples = (size_t)INTEGER(count_dims)[1];
    const size_t n_coefficients = (size_t)INTEGER(design_dims)[1];
    if ((size_t)INTEGER(design_dims)[0] != n_samples)
        Rf_error("design must contain one row per count sample");
    if (XLENGTH(dispersions) != n_genes || XLENGTH(base_mean) != n_genes ||
        XLENGTH(disp_fit) != n_genes)
        Rf_error("dispersion and base-mean vectors must have one value per gene");
    const int selected_coefficient = Rf_asInteger(coefficient);
    const int selected_intercept = Rf_asInteger(intercept_index);
    if (selected_coefficient < 1 || (size_t)selected_coefficient > n_coefficients)
        Rf_error("coef must select a design coefficient");
    SEXP beta = PROTECT(Rf_allocVector(REALSXP, n_genes));
    SEXP beta_se = PROTECT(Rf_allocVector(REALSXP, n_genes));
    SEXP prior_variance = PROTECT(Rf_allocVector(REALSXP, n_coefficients));
    SEXP stat = PROTECT(Rf_allocVector(REALSXP, n_genes));
    SEXP pvalue = PROTECT(Rf_allocVector(REALSXP, n_genes));
    check_rust_backend(rsdeseq2_r_lfc_shrink_normal(
        n_genes, n_samples, REAL(counts),
        size_factors == R_NilValue ? NULL : REAL(size_factors),
        normalization_factors == R_NilValue ? NULL : REAL(normalization_factors),
        observation_weights == R_NilValue ? NULL : REAL(observation_weights),
        n_coefficients, REAL(design), REAL(dispersions), REAL(base_mean), REAL(disp_fit),
        (size_t)(selected_coefficient - 1), (ptrdiff_t)selected_intercept,
        Rf_asReal(lfc_threshold), REAL(beta), REAL(beta_se), REAL(prior_variance),
        REAL(stat), REAL(pvalue)));
    SEXP out = PROTECT(Rf_allocVector(VECSXP, 5));
    SET_VECTOR_ELT(out, 0, beta);
    SET_VECTOR_ELT(out, 1, beta_se);
    SET_VECTOR_ELT(out, 2, prior_variance);
    SET_VECTOR_ELT(out, 3, stat);
    SET_VECTOR_ELT(out, 4, pvalue);
    SEXP names = PROTECT(Rf_allocVector(STRSXP, 5));
    SET_STRING_ELT(names, 0, Rf_mkChar("beta"));
    SET_STRING_ELT(names, 1, Rf_mkChar("betaSE"));
    SET_STRING_ELT(names, 2, Rf_mkChar("betaPriorVar"));
    SET_STRING_ELT(names, 3, Rf_mkChar("stat"));
    SET_STRING_ELT(names, 4, Rf_mkChar("pvalue"));
    Rf_setAttrib(out, R_NamesSymbol, names);
    UNPROTECT(7);
    return out;
}

SEXP rsdeseq2_lfc_shrink_normal_contrast(SEXP counts, SEXP design,
                                         SEXP size_factors, SEXP normalization_factors,
                                         SEXP observation_weights, SEXP dispersions,
                                         SEXP base_mean, SEXP disp_fit,
                                         SEXP standard_beta, SEXP standard_factor_groups,
                                         SEXP expanded_prior_map, SEXP factor_groups,
                                         SEXP contrast, SEXP intercept_index,
                                         SEXP standard_intercept_index,
                                         SEXP lfc_threshold) {
    SEXP count_dims = Rf_getAttrib(counts, R_DimSymbol);
    SEXP design_dims = Rf_getAttrib(design, R_DimSymbol);
    if (Rf_length(count_dims) != 2 || Rf_length(design_dims) != 2)
        Rf_error("counts and expanded design must be matrices");
    const size_t n_genes = (size_t)INTEGER(count_dims)[0];
    const size_t n_samples = (size_t)INTEGER(count_dims)[1];
    const size_t n_coefficients = (size_t)INTEGER(design_dims)[1];
    SEXP standard_dims = Rf_getAttrib(standard_beta, R_DimSymbol);
    if (Rf_length(standard_dims) != 2 || (size_t)INTEGER(standard_dims)[0] != n_genes)
        Rf_error("standard beta estimates must be a genes by coefficients matrix");
    const size_t n_standard_coefficients = (size_t)INTEGER(standard_dims)[1];
    if ((size_t)INTEGER(design_dims)[0] != n_samples)
        Rf_error("expanded design must contain one row per count sample");
    if (XLENGTH(dispersions) != n_genes || XLENGTH(base_mean) != n_genes ||
        XLENGTH(disp_fit) != n_genes || XLENGTH(standard_factor_groups) != n_standard_coefficients ||
        XLENGTH(expanded_prior_map) != n_coefficients || XLENGTH(factor_groups) != n_coefficients ||
        XLENGTH(contrast) != n_coefficients)
        Rf_error("dispersion, prior, factor-group, and contrast dimensions do not match");
    SEXP beta = PROTECT(Rf_allocVector(REALSXP, n_genes));
    SEXP beta_se = PROTECT(Rf_allocVector(REALSXP, n_genes));
    SEXP prior_variance = PROTECT(Rf_allocVector(REALSXP, n_coefficients));
    SEXP stat = PROTECT(Rf_allocVector(REALSXP, n_genes));
    SEXP pvalue = PROTECT(Rf_allocVector(REALSXP, n_genes));
    check_rust_backend(rsdeseq2_r_lfc_shrink_normal_contrast(
        n_genes, n_samples, REAL(counts),
        size_factors == R_NilValue ? NULL : REAL(size_factors),
        normalization_factors == R_NilValue ? NULL : REAL(normalization_factors),
        observation_weights == R_NilValue ? NULL : REAL(observation_weights),
        n_coefficients, REAL(design), n_standard_coefficients, REAL(standard_beta),
        INTEGER(standard_factor_groups), INTEGER(expanded_prior_map),
        REAL(dispersions), REAL(base_mean), REAL(disp_fit),
        INTEGER(factor_groups), REAL(contrast), (ptrdiff_t)(Rf_asInteger(intercept_index) - 1),
        (ptrdiff_t)(Rf_asInteger(standard_intercept_index) - 1),
        Rf_asReal(lfc_threshold), REAL(beta), REAL(beta_se), REAL(prior_variance),
        REAL(stat), REAL(pvalue)));
    SEXP out = PROTECT(Rf_allocVector(VECSXP, 5));
    SET_VECTOR_ELT(out, 0, beta);
    SET_VECTOR_ELT(out, 1, beta_se);
    SET_VECTOR_ELT(out, 2, prior_variance);
    SET_VECTOR_ELT(out, 3, stat);
    SET_VECTOR_ELT(out, 4, pvalue);
    SEXP names = PROTECT(Rf_allocVector(STRSXP, 5));
    SET_STRING_ELT(names, 0, Rf_mkChar("beta"));
    SET_STRING_ELT(names, 1, Rf_mkChar("betaSE"));
    SET_STRING_ELT(names, 2, Rf_mkChar("betaPriorVar"));
    SET_STRING_ELT(names, 3, Rf_mkChar("stat"));
    SET_STRING_ELT(names, 4, Rf_mkChar("pvalue"));
    Rf_setAttrib(out, R_NamesSymbol, names);
    UNPROTECT(7);
    return out;
}
SEXP rsdeseq2_apply_cooks_cutoff(SEXP pvalue,
                                 SEXP max_cooks,
                                 SEXP cooks_cutoff,
                                 SEXP counts,
                                 SEXP cooks,
                                 SEXP low_count_heuristic) {
    const R_xlen_t n_genes = Rf_xlength(pvalue);
    const int use_low_count = LOGICAL(low_count_heuristic)[0] == TRUE;
    size_t n_samples = 0;
    const double *count_values = NULL, *cooks_values = NULL;
    if (use_low_count) {
        SEXP dims = Rf_getAttrib(counts, R_DimSymbol);
        if (Rf_length(dims) != 2) {
            Rf_error("counts must be a matrix");
        }
        if ((R_xlen_t)INTEGER(dims)[0] != n_genes || (size_t)INTEGER(Rf_getAttrib(cooks, R_DimSymbol))[0] != (size_t)n_genes) {
            Rf_error("counts must have one row per p-value");
        }
        n_samples = (size_t)INTEGER(dims)[1];
        count_values = REAL(counts);
        cooks_values = REAL(cooks);
    }
    SEXP masked_pvalue = PROTECT(Rf_allocVector(REALSXP, n_genes));
    SEXP cooks_outlier = PROTECT(Rf_allocVector(LGLSXP, n_genes));
    uint8_t *outlier = (uint8_t *)R_alloc((size_t)n_genes, sizeof(uint8_t));
    check_rust_backend(rsdeseq2_r_cooks_mask((size_t)n_genes, REAL(pvalue), REAL(max_cooks),
        REAL(cooks_cutoff)[0], count_values, n_samples, cooks_values,
        (uint8_t)use_low_count, REAL(masked_pvalue), outlier));
    for (R_xlen_t gene = 0; gene < n_genes; ++gene) {
        LOGICAL(cooks_outlier)[gene] = outlier[gene] == UINT8_MAX ? NA_LOGICAL : (outlier[gene] ? TRUE : FALSE);
    }

    SEXP names = Rf_getAttrib(pvalue, R_NamesSymbol);
    if (names != R_NilValue) {
        Rf_setAttrib(masked_pvalue, R_NamesSymbol, names);
    }

    SEXP out = PROTECT(Rf_allocVector(VECSXP, 2));
    SET_VECTOR_ELT(out, 0, masked_pvalue);
    SET_VECTOR_ELT(out, 1, cooks_outlier);
    SEXP out_names = PROTECT(Rf_allocVector(STRSXP, 2));
    SET_STRING_ELT(out_names, 0, Rf_mkChar("pvalue"));
    SET_STRING_ELT(out_names, 1, Rf_mkChar("cooksOutlier"));
    Rf_setAttrib(out, R_NamesSymbol, out_names);

    UNPROTECT(4);
    return out;
}

static const R_CallMethodDef call_methods[] = {
    {"rsdeseq2_placeholder", (DL_FUNC)&rsdeseq2_placeholder, 0},
    {"rsdeseq2_diagnostic_schema", (DL_FUNC)&rsdeseq2_diagnostic_schema, 0},
    {"rsdeseq2_estimate_size_factors", (DL_FUNC)&rsdeseq2_estimate_size_factors, 4},
    {"rsdeseq2_normalized_counts", (DL_FUNC)&rsdeseq2_normalized_counts, 3},
    {"rsdeseq2_base_mean", (DL_FUNC)&rsdeseq2_base_mean, 3},
    {"rsdeseq2_base_metadata", (DL_FUNC)&rsdeseq2_base_metadata, 4},
    {"rsdeseq2_fit", (DL_FUNC)&rsdeseq2_fit, 10},
    {"rsdeseq2_transform", (DL_FUNC)&rsdeseq2_transform, 16},
    {"rsdeseq2_lfc_shrink_normal", (DL_FUNC)&rsdeseq2_lfc_shrink_normal, 11},
    {"rsdeseq2_lfc_shrink_normal_contrast", (DL_FUNC)&rsdeseq2_lfc_shrink_normal_contrast, 16},
    {"rsdeseq2_apply_cooks_cutoff", (DL_FUNC)&rsdeseq2_apply_cooks_cutoff, 6},
    {NULL, NULL, 0},
};

void R_init_rsdeseq2(DllInfo *dll) {
    R_registerRoutines(dll, NULL, call_methods, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
