
<!-- README.md is generated from README.Rmd. Please edit that file -->

# metahelpr

<!-- badges: start -->
<!-- badges: end -->

The goal of `metahelpr` is to offer tutorials and helper functions for
conducting meta-analysis using the `metafor` package (Viechtbauer 2010).
The package contains some helper functions that I often use when I
conduct meta-analysis, but I offer no guarantees that they’ll actually
work for you because I’m not a particularly experienced programmer. Feel
free to contribute and make them better and more robust.

## Installation

You can install the development version of `metahelpr` from
[GitHub](https://github.com/) with:

``` r
# install.packages("remotes")
remotes::install_github("profandyfield/metahelpr")
```

## Interactive tutorials

I recommend working through [this playlist of
tutorials](https://youtube.com/playlist?list=PLEzw67WWDg83weG3idsgy4wuOIJAashA2&si=PiI-sDvqc1DkaWOq)
on how to install, set up and work within
<img src="./data-raw/images/r_logo.png" width="18"/> and
<img src="./data-raw/images/rstudio_logo.png" width="48"/> before
starting the interactive tutorials.

- **meta_d**: A tutorial on using the [`metafor`
  package](https://wviechtb.github.io/metafor/) to conduct a
  meta-analysis using Cohen’s $\hat{d}$ as the effect size measure.
- **meta_r**: A tutorial on using the [`metafor`
  package](https://wviechtb.github.io/metafor/) to conduct a
  meta-analysis using Pearson’s $r$ as the effect size measure.

## Data files

- `brewin_2024`: Raw data from the meta-analysis by Brewin and Field
  (2024).
- `brewin_es`: Reduced dataset from the meta-analysis by Brewin and
  Field (2024) that includes effect sizes.
- `pearce_2016`: A selection of variables from the dataset from the
  meta-analysis by Pearce and Field (2016)

## Helper functions

The package contains some (probably badly written and easily breakable)
helper functions that I often use when I conduct meta-analysis.

- `d_from_r`: estimate Cohen’s *d* based on a biserial correlation
  coefficient, or when a biserial correlation isn’t available the
  conversion uses the method described by Mathur and VanderWeele (2020).
- `d_to_g`: converts Cohen’s *d* to the unbiased Hedges’ *g*.
- `forest_add_het_stats`: Add heterogeneity statistics to a forest plot
- `forest_subgroups`: plots a forest plot where effect sizes are split
  by subcategories of a factor.
- `get_mas`: Fit individual meta-analyses for each level of a
  categorical predictor and (optionally) collate the results into a
  tabulated form for printing in a quarto document.
- `get_mod_mas`: does the same thing as `get_mas()` but allows the user
  to specify a moderator/predictor within the individual meta-analyses.
  So, it fits individual meta-analyses models with a moderator specified
  using the `rma.mv()` function from the `metafor` package within each
  level of a categorical variable. Optionally, the results can be
  collated into a tabulated form for printing in a quarto document.
- `get_pbm`: fits and collates publication bias models across categories
  of a predictor variable. It is assumed that you will supply two
  vectors of values one representing moderate publication bias and the
  other representing severe.
- `plot_bubble`: create a bubble plot (using `ggplot2`) based on a
  `regtest()` object.
- `pooled_var`: computes a full sample variance estimate based on means
  and variances from two subgroups.
- `regtest_tbl`: puts the results of `metafor::regtest()` into a tibble
  for reporting.
- `report_het`: collates information from heterogeneity tests and
  outputs text that summarizes the results in a format that will render
  nicely in quarto.
- `report_mod`: outputs text that reports (and renders nicely in
  quarto/Rmarkdown) the omnibus statistical tests from a moderation
  model.
- `report_par_tbl`: outputs a tibble of the table of coefficients of a
  meta-analyses object (`rma` or `rma.mv`) created using the `metafor`
  package. This function will mostly be useful for models containing
  predictors of effect sizes (so-called meta-regression).
- `report_pars`: outputs text that reports (and renders nicely in
  quarto/Rmarkdown) the individual effects from a meta-analysis model.
- `var_d_from_r`: Estimate the sampling variance of Cohen’s *d* based on
  a correlation coefficient

## References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-brewin_meta-analysis_2024" class="csl-entry">

Brewin, C. R., and Andy P. Field. 2024. “Meta-Analysis Shows Trauma
Memories in PTSD Lack Coherence: A Response to Taylor Et Al. (2022).”
*Clinical Psychological Science*.

</div>

<div id="ref-mathur2020" class="csl-entry">

Mathur, Maya B., and Tyler J. VanderWeele. 2020. “A Simple,
Interpretable Conversion from Pearson’s Correlation to Cohen’s
\<Em\>d\</Em\> for Continuous Exposures.” *Epidemiology* 31 (2): e16–17.
<https://doi.org/10.1097/EDE.0000000000001111>.

</div>

<div id="ref-pearce2016" class="csl-entry">

Pearce, Laura J., and Andy P. Field. 2016. “The Impact of ‘Scary’ TV and
Film on Children’s Internalizing Emotions: A Meta-Analysis.” *Human
Communication Research* 42 (January): 98–121.
<https://doi.org/doi:10.1111/hcre.12069>.

</div>

<div id="ref-metafor" class="csl-entry">

Viechtbauer, Wolfgang. 2010. “Conducting Meta-Analyses in {R} with the
{Metafor} Package.” <https://doi.org/10.18637/jss.v036.i03>.

</div>

</div>
