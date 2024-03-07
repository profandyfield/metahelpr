#' metahelpr: Tutorials and helper functions for meta-analysis using metafor
#'
#' @description
#'
#' The goal of `metahelpr` is to offer tutorials and helper functions for conducting meta-analysis using the `metafor` package (Viechtbauer, 2010).
#'
#'
#' @section Interactive tutorials:
#'
#' **Getting started**:
#'
#' I recommend working through [this playlist of tutorials](https://youtube.com/playlist?list=PLEzw67WWDg83weG3idsgy4wuOIJAashA2&si=PiI-sDvqc1DkaWOq) on how to install, set up and work within R, RStudio and Quarto before starting the interactive tutorials.
#'
#' **Running a tutorial**:
#'
#' To run each tutorial execute
#'
#' `learnr::run_tutorial("name_of_tutorial", package = "metahelpr")`
#'
#' Replacing `name_of_tutorial` with the name in bold below. For example, to load the tutorial `meta_d` execute:
#'
#' `learnr::run_tutorial("meta_d", package = "metahelpr")`
#'
#' -   **meta_d**: A tutorial on using the [`metafor` package](https://wviechtb.github.io/metafor/) to conduct a meta-analysis using Cohen's \eqn{\hat{d}} as the effect size measure.
#' -   **meta_r**: A tutorial on using the [`metafor` package](https://wviechtb.github.io/metafor/) to conduct a meta-analysis using Pearson's $r$ as the effect size measure.
#'
#' @section Datasets:
#'
#' * [brewin_2024]: Raw data from the meta-analysis by Brewin & Field (2024).
#' * [brewin_es]: Reduced dataset from the meta-analysis by Brewin & Field (2024) that includes effect sizes.
#' * [pearce_2016]: A selection of variables from the dataset from the meta-analysis by Pearce & Field (2016).
#'
#' @section Helper functions:
#'
#' The package contains some (probably badly written and easily breakable) helper functions that I often use when I conduct meta-analysis.
#'
#' - [d_from_r]: estimate Cohen's *d* based on a biserial correlation coefficient, or when a biserial correlation isn't available the conversion uses the method described by  Mathur & VanderWeele (2019).
#' - [d_to_g]: converts Cohen's *d* to the unbiased Hedges' *g*.
#' - [forest_add_het_stats]: Add heterogeneity statistics to a forest plot
#' - [forest_subgroups]: plots a forest plot where effect sizes are split by subcategories of a factor.
#' - [get_mas]: Fit individual meta-analyses for each level of a categorical predictor and (optionally) collate the results into a tabulated form for printing in a quarto document.
#' - [get_mod_mas]: does the same thing as [get_mas] but allows the user to specify a moderator/predictor within the individual meta-analyses. So, it fits individual meta-analyses models with a moderator specified using the `rma.mv()` function from the `metafor` package within each level of a categorical variable. Optionally, the results can be collated into a tabulated form for printing in a quarto document.
#' - [get_pbm]: fits and collates publication bias models across categories of a predictor variable. It is assumed that you will supply two vectors of values one representing moderate publication bias and the other representing severe.
#' - [plot_bubble]: create a bubble plot (using `ggplot2`) based on a `regtest()` object.
#' - [pooled_var]: computes a full sample variance estimate based on means and variances from two subgroups.
#' - [regtest_tbl]: puts the results of `metafor::regtest()` into a tibble for reporting.
#' - [report_het]: collates information from heterogeneity tests and outputs text that summarizes the results in a format that will render nicely in quarto.
#' - [report_mod]: outputs text that reports (and renders nicely in quarto/Rmarkdown) the omnibus statistical tests from a moderation model.
#' - [report_par_tbl]: outputs a tibble of the table of coefficients of a meta-analyses object (`rma` or `rma.mv`) created using the `metafor` package. This function will mostly be useful for models containing predictors of effect sizes (so-called meta-regression).
#' - [report_pars]: outputs text that reports (and renders nicely in quarto/Rmarkdown) the individual effects from a meta-analysis model.
#' - [var_d_from_r]: Estimate the sampling variance of Cohen's *d* based on a correlation coefficient

#' @section References:
#' - Brewin, C. R., & Field, A. P. (2024). Meta-analysis shows trauma memories in PTSD lack coherence: A response to Taylor et al. (2022). *Clinical Psychological Science*.
#' - Mathur, M. B., & VanderWeele, T. J. (2020). A simple, interpretable conversion from Pearson's correlation to Cohen's *d* for continuous exposures. *Epidemiology*, 31(2), e16–e17. [doi.org/10.1097/EDE.0000000000001111](https://doi.org/10.1097/EDE.0000000000001111)
#' - Pearce, L. J., & Field, A. P. (2016). The impact of 'scary' TV and film on children’s internalizing emotions: A meta-analysis. *Human Communication Research*, 42, 98–121. [https://doi.org/doi:10.1111/hcre.12069](https://doi.org/doi:10.1111/hcre.12069)
#' - Viechtbauer, W. (2010). Conducting Meta-Analyses in R with the metafor Package. *Journal of Statistical Software*, 36, 1–48.
#'
#' @import learnr
#' @name metahelpr
#'
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL
