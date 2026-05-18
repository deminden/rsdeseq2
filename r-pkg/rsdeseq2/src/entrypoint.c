#include <R.h>
#include <R_ext/Rdynload.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <math.h>
#include <stdlib.h>

static int compare_double_ascending(const void *lhs, const void *rhs) {
    const double left = *(const double *)lhs;
    const double right = *(const double *)rhs;
    return (left > right) - (left < right);
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

SEXP rsdeseq2_estimate_size_factors(SEXP counts,
                                    SEXP log_geo_means,
                                    SEXP control_genes,
                                    SEXP stabilize) {
    SEXP dims = Rf_getAttrib(counts, R_DimSymbol);
    if (Rf_length(dims) != 2) {
        Rf_error("counts must be a matrix");
    }
    const int n_genes = INTEGER(dims)[0];
    const int n_samples = INTEGER(dims)[1];
    const R_xlen_t n_control = Rf_xlength(control_genes);
    const double *count_values = REAL(counts);
    const double *log_geo_mean_values = REAL(log_geo_means);
    const int *control_values = INTEGER(control_genes);

    double *log_ratios = (double *)R_alloc((size_t)n_control, sizeof(double));
    SEXP out = PROTECT(Rf_allocVector(REALSXP, n_samples));
    for (int sample = 0; sample < n_samples; ++sample) {
        R_xlen_t n_usable = 0;
        for (R_xlen_t control = 0; control < n_control; ++control) {
            const int gene = control_values[control] - 1;
            if (gene < 0 || gene >= n_genes) {
                Rf_error("controlGenes contains an out-of-range row index");
            }
            const R_xlen_t index = gene + (R_xlen_t)sample * n_genes;
            const double count = count_values[index];
            const double log_geo_mean = log_geo_mean_values[gene];
            if (isfinite(log_geo_mean) && count > 0.0) {
                log_ratios[n_usable] = log(count) - log_geo_mean;
                ++n_usable;
            }
        }
        if (n_usable == 0) {
            Rf_error("sample %d has no usable positive count ratios", sample + 1);
        }
        qsort(log_ratios, (size_t)n_usable, sizeof(double), compare_double_ascending);
        const double median = n_usable % 2 == 1
                                  ? log_ratios[n_usable / 2]
                                  : (log_ratios[n_usable / 2 - 1] + log_ratios[n_usable / 2]) / 2.0;
        const double size_factor = exp(median);
        if (!isfinite(size_factor) || size_factor <= 0.0) {
            Rf_error("estimated size factors must be finite and positive");
        }
        REAL(out)[sample] = size_factor;
    }

    if (LOGICAL(stabilize)[0]) {
        double log_sum = 0.0;
        for (int sample = 0; sample < n_samples; ++sample) {
            log_sum += log(REAL(out)[sample]);
        }
        const double scale = exp(log_sum / (double)n_samples);
        if (!isfinite(scale) || scale <= 0.0) {
            Rf_error("cannot stabilize size factors to geometric mean one");
        }
        for (int sample = 0; sample < n_samples; ++sample) {
            REAL(out)[sample] /= scale;
        }
    }

    SEXP dimnames = Rf_getAttrib(counts, R_DimNamesSymbol);
    if (dimnames != R_NilValue && Rf_length(dimnames) >= 2) {
        SEXP sample_names = VECTOR_ELT(dimnames, 1);
        if (sample_names != R_NilValue) {
            Rf_setAttrib(out, R_NamesSymbol, sample_names);
        }
    }

    UNPROTECT(1);
    return out;
}

SEXP rsdeseq2_normalized_counts(SEXP counts,
                                SEXP size_factors,
                                SEXP normalization_factors) {
    SEXP dims = Rf_getAttrib(counts, R_DimSymbol);
    if (Rf_length(dims) != 2) {
        Rf_error("counts must be a matrix");
    }
    const int n_genes = INTEGER(dims)[0];
    const int n_samples = INTEGER(dims)[1];
    const R_xlen_t n_values = (R_xlen_t)n_genes * n_samples;
    const double *count_values = REAL(counts);
    const double *size_factor_values = size_factors == R_NilValue ? NULL : REAL(size_factors);
    const double *normalization_factor_values =
        normalization_factors == R_NilValue ? NULL : REAL(normalization_factors);

    SEXP out = PROTECT(Rf_allocVector(REALSXP, n_values));
    for (int sample = 0; sample < n_samples; ++sample) {
        for (int gene = 0; gene < n_genes; ++gene) {
            const R_xlen_t index = gene + (R_xlen_t)sample * n_genes;
            const double denom = normalization_factor_values == NULL
                                     ? size_factor_values[sample]
                                     : normalization_factor_values[index];
            REAL(out)[index] = count_values[index] / denom;
        }
    }

    Rf_setAttrib(out, R_DimSymbol, dims);
    SEXP dimnames = Rf_getAttrib(counts, R_DimNamesSymbol);
    if (dimnames != R_NilValue) {
        Rf_setAttrib(out, R_DimNamesSymbol, dimnames);
    }

    UNPROTECT(1);
    return out;
}

SEXP rsdeseq2_base_mean(SEXP counts,
                        SEXP size_factors,
                        SEXP normalization_factors) {
    SEXP dims = Rf_getAttrib(counts, R_DimSymbol);
    if (Rf_length(dims) != 2) {
        Rf_error("counts must be a matrix");
    }
    const int n_genes = INTEGER(dims)[0];
    const int n_samples = INTEGER(dims)[1];
    const double *count_values = REAL(counts);
    const double *size_factor_values = size_factors == R_NilValue ? NULL : REAL(size_factors);
    const double *normalization_factor_values =
        normalization_factors == R_NilValue ? NULL : REAL(normalization_factors);

    SEXP out = PROTECT(Rf_allocVector(REALSXP, n_genes));
    for (int gene = 0; gene < n_genes; ++gene) {
        double sum = 0.0;
        for (int sample = 0; sample < n_samples; ++sample) {
            const R_xlen_t index = gene + (R_xlen_t)sample * n_genes;
            const double denom = normalization_factor_values == NULL
                                     ? size_factor_values[sample]
                                     : normalization_factor_values[index];
            sum += count_values[index] / denom;
        }
        REAL(out)[gene] = sum / (double)n_samples;
    }

    SEXP dimnames = Rf_getAttrib(counts, R_DimNamesSymbol);
    if (dimnames != R_NilValue && Rf_length(dimnames) >= 1) {
        SEXP row_names = VECTOR_ELT(dimnames, 0);
        if (row_names != R_NilValue) {
            Rf_setAttrib(out, R_NamesSymbol, row_names);
        }
    }

    UNPROTECT(1);
    return out;
}

SEXP rsdeseq2_base_metadata(SEXP counts,
                            SEXP size_factors,
                            SEXP normalization_factors,
                            SEXP weights) {
    SEXP dims = Rf_getAttrib(counts, R_DimSymbol);
    if (Rf_length(dims) != 2) {
        Rf_error("counts must be a matrix");
    }
    const int n_genes = INTEGER(dims)[0];
    const int n_samples = INTEGER(dims)[1];
    const double *count_values = REAL(counts);
    const double *size_factor_values = size_factors == R_NilValue ? NULL : REAL(size_factors);
    const double *normalization_factor_values =
        normalization_factors == R_NilValue ? NULL : REAL(normalization_factors);
    const double *weight_values = weights == R_NilValue ? NULL : REAL(weights);

    SEXP base_mean = PROTECT(Rf_allocVector(REALSXP, n_genes));
    SEXP base_var = PROTECT(Rf_allocVector(REALSXP, n_genes));
    SEXP all_zero = PROTECT(Rf_allocVector(LGLSXP, n_genes));

    for (int gene = 0; gene < n_genes; ++gene) {
        double sum = 0.0;
        int zero_row = 1;
        for (int sample = 0; sample < n_samples; ++sample) {
            const R_xlen_t index = gene + (R_xlen_t)sample * n_genes;
            const double count = count_values[index];
            const double denom = normalization_factor_values == NULL
                                     ? size_factor_values[sample]
                                     : normalization_factor_values[index];
            double value = count / denom;
            if (weight_values != NULL) {
                value *= weight_values[index];
            }
            if (count != 0.0) {
                zero_row = 0;
            }
            sum += value;
        }

        const double mean = sum / (double)n_samples;
        REAL(base_mean)[gene] = mean;
        LOGICAL(all_zero)[gene] = zero_row;

        if (n_samples <= 1) {
            REAL(base_var)[gene] = R_NaN;
        } else {
            double sum_sq = 0.0;
            for (int sample = 0; sample < n_samples; ++sample) {
                const R_xlen_t index = gene + (R_xlen_t)sample * n_genes;
                const double denom = normalization_factor_values == NULL
                                         ? size_factor_values[sample]
                                         : normalization_factor_values[index];
                double value = count_values[index] / denom;
                if (weight_values != NULL) {
                    value *= weight_values[index];
                }
                const double delta = value - mean;
                sum_sq += delta * delta;
            }
            REAL(base_var)[gene] = sum_sq / (double)(n_samples - 1);
        }
    }

    SEXP out = PROTECT(Rf_allocVector(VECSXP, 3));
    SET_VECTOR_ELT(out, 0, base_mean);
    SET_VECTOR_ELT(out, 1, base_var);
    SET_VECTOR_ELT(out, 2, all_zero);

    SEXP names = PROTECT(Rf_allocVector(STRSXP, 3));
    SET_STRING_ELT(names, 0, Rf_mkChar("baseMean"));
    SET_STRING_ELT(names, 1, Rf_mkChar("baseVar"));
    SET_STRING_ELT(names, 2, Rf_mkChar("allZero"));
    Rf_setAttrib(out, R_NamesSymbol, names);

    SEXP dimnames = Rf_getAttrib(counts, R_DimNamesSymbol);
    SEXP row_names = R_NilValue;
    if (dimnames != R_NilValue && Rf_length(dimnames) >= 1) {
        row_names = VECTOR_ELT(dimnames, 0);
    }
    if (row_names == R_NilValue) {
        row_names = PROTECT(Rf_allocVector(INTSXP, n_genes));
        for (int gene = 0; gene < n_genes; ++gene) {
            INTEGER(row_names)[gene] = gene + 1;
        }
        Rf_setAttrib(out, R_RowNamesSymbol, row_names);
        UNPROTECT(1);
    } else {
        Rf_setAttrib(out, R_RowNamesSymbol, row_names);
    }

    SEXP classes = PROTECT(Rf_allocVector(STRSXP, 1));
    SET_STRING_ELT(classes, 0, Rf_mkChar("data.frame"));
    Rf_setAttrib(out, R_ClassSymbol, classes);

    UNPROTECT(6);
    return out;
}

SEXP rsdeseq2_apply_cooks_cutoff(SEXP pvalue,
                                 SEXP max_cooks,
                                 SEXP cooks_cutoff,
                                 SEXP counts,
                                 SEXP cooks,
                                 SEXP low_count_heuristic) {
    const R_xlen_t n_genes = Rf_xlength(pvalue);
    const double *pvalue_values = REAL(pvalue);
    const double *max_cooks_values = REAL(max_cooks);
    const int use_cutoff = cooks_cutoff != R_NilValue;
    const int use_low_count = LOGICAL(low_count_heuristic)[0] == TRUE;
    const double cutoff = use_cutoff ? REAL(cooks_cutoff)[0] : 0.0;

    int n_samples = 0;
    const double *count_values = NULL;
    const double *cooks_values = NULL;
    if (use_low_count) {
        SEXP dims = Rf_getAttrib(counts, R_DimSymbol);
        if (Rf_length(dims) != 2) {
            Rf_error("counts must be a matrix");
        }
        if ((R_xlen_t)INTEGER(dims)[0] != n_genes) {
            Rf_error("counts must have one row per p-value");
        }
        n_samples = INTEGER(dims)[1];
        count_values = REAL(counts);
        cooks_values = REAL(cooks);
    }

    SEXP masked_pvalue = PROTECT(Rf_allocVector(REALSXP, n_genes));
    SEXP cooks_outlier = PROTECT(Rf_allocVector(LGLSXP, n_genes));
    for (R_xlen_t gene = 0; gene < n_genes; ++gene) {
        REAL(masked_pvalue)[gene] = pvalue_values[gene];
        LOGICAL(cooks_outlier)[gene] = NA_LOGICAL;

        if (!use_cutoff || ISNAN(max_cooks_values[gene])) {
            continue;
        }
        int outlier = max_cooks_values[gene] > cutoff;

        if (outlier && use_low_count) {
            int max_sample = -1;
            double max_cook = R_NegInf;
            for (int sample = 0; sample < n_samples; ++sample) {
                const R_xlen_t index = gene + (R_xlen_t)sample * n_genes;
                const double cook = cooks_values[index];
                if (!ISNA(cook) && isfinite(cook) && (max_sample < 0 || cook > max_cook)) {
                    max_sample = sample;
                    max_cook = cook;
                }
            }
            if (max_sample >= 0) {
                const double out_count = count_values[gene + (R_xlen_t)max_sample * n_genes];
                int larger_count_samples = 0;
                for (int sample = 0; sample < n_samples; ++sample) {
                    const R_xlen_t index = gene + (R_xlen_t)sample * n_genes;
                    if (count_values[index] > out_count) {
                        ++larger_count_samples;
                    }
                }
                if (larger_count_samples >= 3) {
                    outlier = 0;
                }
            }
        }

        LOGICAL(cooks_outlier)[gene] = outlier ? TRUE : FALSE;
        if (outlier) {
            REAL(masked_pvalue)[gene] = NA_REAL;
        }
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
    {"rsdeseq2_apply_cooks_cutoff", (DL_FUNC)&rsdeseq2_apply_cooks_cutoff, 6},
    {NULL, NULL, 0},
};

void R_init_rsdeseq2(DllInfo *dll) {
    R_registerRoutines(dll, NULL, call_methods, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
