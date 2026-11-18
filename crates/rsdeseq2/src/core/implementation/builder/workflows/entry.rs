//! Default and named top-level fit workflows.

use super::super::*;

impl DeseqBuilder {
    /// Run the currently implemented DESeq-like workflow.
    ///
    /// For `test=Wald`, this selects DESeq2-compatible dispersion means and
    /// reports the last design coefficient, matching DESeq2's default
    /// coefficient selection shape. For `test=Lrt`, callers must first store a
    /// reduced design with [`DeseqBuilder::reduced_design`].
    #[deprecated(
        since = "0.6.0",
        note = "use DeseqBuilder::run(&Analysis::new(counts, design)), adding .lrt(reduced) when the builder test is LRT"
    )]
    pub fn fit(&self, counts: &CountMatrix, design: &DesignMatrix) -> Result<DeseqFit, DeseqError> {
        let test = self.top_level_test()?;
        self.run(&Analysis::new(counts, design).test(test))?
            .into_fit_only()
    }

    /// Build a supported formula design from stored model-frame metadata and
    /// run the currently implemented DESeq-like workflow, returning only the
    /// fit state.
    #[deprecated(
        since = "0.6.0",
        note = "use DeseqBuilder::run(&Analysis::formula(counts, formula)), adding .lrt(reduced) when the builder test is LRT"
    )]
    pub fn fit_formula(&self, counts: &CountMatrix, formula: &str) -> Result<DeseqFit, DeseqError> {
        let test = self.top_level_test()?;
        self.run(&Analysis::formula(counts, formula).test(test))?
            .into_fit_only()
    }

    /// Run the currently implemented DESeq-like workflow and return result rows.
    ///
    /// This is the result-table-producing companion to [`DeseqBuilder::fit`].
    #[deprecated(
        since = "0.6.0",
        note = "use DeseqBuilder::run(&Analysis::new(counts, design)), adding .lrt(reduced) when the builder test is LRT"
    )]
    pub fn fit_with_results(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        let test = self.top_level_test()?;
        self.run(&Analysis::new(counts, design).test(test))?
            .into_results_pair()
    }

    /// Build a supported formula design from stored model-frame metadata and
    /// run the currently implemented DESeq-like workflow with result rows.
    ///
    #[deprecated(
        since = "0.6.0",
        note = "use DeseqBuilder::run(&Analysis::formula(counts, formula)), adding .lrt(reduced) when the builder test is LRT"
    )]
    pub fn fit_formula_with_results(
        &self,
        counts: &CountMatrix,
        formula: &str,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        let test = self.top_level_test()?;
        self.run(&Analysis::formula(counts, formula).test(test))?
            .into_results_pair()
    }

    /// Run the top-level Wald workflow and report a named design coefficient.
    #[deprecated(
        since = "0.6.0",
        note = "use DeseqBuilder::run(&Analysis::new(counts, design).result(ResultColumn::Name(coefficient_name))), adding .lrt(reduced) when the builder test is LRT"
    )]
    pub fn fit_with_results_name(
        &self,
        counts: &CountMatrix,
        design: &DesignMatrix,
        coefficient_name: &str,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        let test = self.top_level_test()?;
        self.run(
            &Analysis::new(counts, design)
                .test(test)
                .result(ResultColumn::Name(coefficient_name)),
        )?
        .into_results_pair()
    }

    /// Build a supported formula design from stored model-frame metadata and
    /// report a named coefficient.
    #[deprecated(
        since = "0.6.0",
        note = "use DeseqBuilder::run(&Analysis::formula(counts, formula).result(ResultColumn::Name(coefficient_name))), adding .lrt(reduced) when the builder test is LRT"
    )]
    pub fn fit_formula_with_results_name(
        &self,
        counts: &CountMatrix,
        formula: &str,
        coefficient_name: &str,
    ) -> Result<(DeseqFit, DeseqResults), DeseqError> {
        let test = self.top_level_test()?;
        self.run(
            &Analysis::formula(counts, formula)
                .test(test)
                .result(ResultColumn::Name(coefficient_name)),
        )?
        .into_results_pair()
    }

    /// Build a supported formula design from stored model-frame metadata and
    /// report a named coefficient, returning only the fit state.
    #[deprecated(
        since = "0.6.0",
        note = "use DeseqBuilder::run(&Analysis::formula(counts, formula).result(ResultColumn::Name(coefficient_name))), adding .lrt(reduced) when the builder test is LRT"
    )]
    pub fn fit_formula_name(
        &self,
        counts: &CountMatrix,
        formula: &str,
        coefficient_name: &str,
    ) -> Result<DeseqFit, DeseqError> {
        let test = self.top_level_test()?;
        self.run(
            &Analysis::formula(counts, formula)
                .test(test)
                .result(ResultColumn::Name(coefficient_name)),
        )?
        .into_fit_only()
    }
}
