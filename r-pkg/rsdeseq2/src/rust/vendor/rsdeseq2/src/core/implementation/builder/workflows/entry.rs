//! Default and named top-level fit workflows.

use super::super::{
    CountMatrix, DeseqBuilder, DeseqError, DeseqFit, DeseqResults, DesignMatrix, TestType,
    default_results_coefficient, resolve_coefficient_index,
};

impl DeseqBuilder {
    /// Run the currently implemented DESeq-like workflow.
    ///
    /// For `test=Wald`, this selects DESeq2-compatible dispersion means and
    /// reports the last design coefficient, matching DESeq2's default
    /// coefficient selection shape. For `test=Lrt`, callers must first store a
    /// reduced design with [`DeseqBuilder::reduced_design`].
    pub fn fit(&self, counts: &CountMatrix, design: &DesignMatrix) -> Result<DeseqFit, DeseqError> {
        if let Some(pool) = self.execution_context.pool_for_entry(self.threads)? {
            return pool.install(|| self.fit(counts, design));
        }

        self.fit_with_results(counts, design)
            .map(|(fit, _results)| fit)
    }

    /// Build a supported formula design from stored model-frame metadata and
    /// run the currently implemented DESeq-like workflow, returning only the
    /// fit state.
    pub fn fit_formula(&self, counts: &CountMatrix, formula: &str) -> Result<DeseqFit, DeseqError> {
        if let Some(pool) = self.execution_context.pool_for_entry(self.threads)? {
            return pool.install(|| self.fit_formula(counts, formula));
        }

        self.fit_formula_with_results(counts, formula)
            .map(|(fit, _results)| fit)
    }

    /// Run the currently implemented DESeq-like workflow and return result rows.
    ///
    /// This is the result-table-producing companion to [`DeseqBuilder::fit`].
    pub fn fit_with_results(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        if let Some(pool) = self.execution_context.pool_for_entry(self.threads)? {
            return pool.install(|| self.fit_with_results(counts, design));
        }

        match self.test {
            TestType::Wald => {
                let coefficient = default_results_coefficient(design)?;
                self.with_automatic_dispersion_means()
                    .fit_wald_glm_mu(counts, design, coefficient)
            }
            TestType::Lrt => {
                let reduced_design = self.reduced_design_for_top_level_lrt()?;
                self.with_automatic_dispersion_means().fit_lrt_with_results(
                    counts,
                    design,
                    reduced_design,
                )
            }
        }
    }

    /// Build a supported formula design from stored model-frame metadata and
    /// run the currently implemented DESeq-like workflow with result rows.
    ///
    pub fn fit_formula_with_results(
        &self,
        counts: &CountMatrix,
        formula: &str,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        if let Some(pool) = self.execution_context.pool_for_entry(self.threads)? {
            return pool.install(|| self.fit_formula_with_results(counts, formula));
        }

        let formula_design = self.expanded_formula_design_with_offsets(formula)?;
        self.with_formula_offsets(counts, &formula_design)?
            .fit_with_results(counts, &formula_design.design.standard_design)
    }

    /// Run the top-level Wald workflow and report a named design coefficient.
    pub fn fit_with_results_name(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        coefficient_name: &str,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        if let Some(pool) = self.execution_context.pool_for_entry(self.threads)? {
            return pool.install(|| self.fit_with_results_name(counts, design, coefficient_name));
        }

        match self.test {
            TestType::Wald => {
                let coefficient = resolve_coefficient_index(design, coefficient_name)?;
                self.with_automatic_dispersion_means()
                    .fit_wald_glm_mu(counts, design, coefficient)
            }
            TestType::Lrt => {
                let reduced_design = self.reduced_design_for_top_level_lrt()?;
                self.with_automatic_dispersion_means()
                    .fit_lrt_with_results_name(counts, design, reduced_design, coefficient_name)
            }
        }
    }

    /// Build a supported formula design from stored model-frame metadata and
    /// report a named coefficient.
    pub fn fit_formula_with_results_name(
        &self,
        counts: &CountMatrix,
        formula: &str,
        coefficient_name: &str,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        if let Some(pool) = self.execution_context.pool_for_entry(self.threads)? {
            return pool
                .install(|| self.fit_formula_with_results_name(counts, formula, coefficient_name));
        }

        let formula_design = self.expanded_formula_design_with_offsets(formula)?;
        self.with_formula_offsets(counts, &formula_design)?
            .fit_with_results_name(
                counts,
                &formula_design.design.standard_design,
                coefficient_name,
            )
    }

    /// Build a supported formula design from stored model-frame metadata and
    /// report a named coefficient, returning only the fit state.
    pub fn fit_formula_name(
        &self,
        counts: &CountMatrix,
        formula: &str,
        coefficient_name: &str,
    ) -> Result<DeseqFit, DeseqError> {
        if let Some(pool) = self.execution_context.pool_for_entry(self.threads)? {
            return pool.install(|| self.fit_formula_name(counts, formula, coefficient_name));
        }

        self.fit_formula_with_results_name(counts, formula, coefficient_name)
            .map(|(fit, _results)| fit)
    }
}
