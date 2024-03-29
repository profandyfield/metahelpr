% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/meta_helpers.R
\name{get_mas}
\alias{get_mas}
\title{Fit individual meta-analyses for each level of a categorical predictor.}
\usage{
get_mas(
  tibble,
  predictor,
  es,
  var_es,
  es_id,
  study_id,
  digits = 2,
  p_digits = 3,
  summary = TRUE
)
}
\arguments{
\item{tibble}{A tibble containing the effect sizes, their variances,a study ID and an effect size ID}

\item{predictor}{Name of the categorical variable for which you want individual meta-analyses}

\item{es}{Name of the variable containing the effect sizes}

\item{var_es}{Name of the variable containing the variance estimate of the effect sizes}

\item{es_id}{Name of the variable that identifies unique effect sizes}

\item{study_id}{Name of the variable that identifies unique studies}

\item{digits}{number of decimal places to print in the output}

\item{p_digits}{number of decimal places for \emph{p}-values}

\item{summary}{if TRUE output a summary table, otherwise store the raw data and models}
}
\value{
If \code{summary = TRUE} the function returns a tibble containing the results of the meta-analysis for each category in a separate row.
Columns show
\itemize{
\item \emph{k}: number of effect sizes
\item \eqn{\sigma_b}: between study variance
\item \eqn{\sigma_w}: within-study variance
\item \emph{q}: Q statistic for heterogeneity,
\item \eqn{q_{df}}: degrees of freedom for the \emph{Q} statistic,
\item \eqn{q_p}: \emph{p}-value for the Q statistic
\item estimate: The pooled effect size estimate
\item ci: the 95\% confidence interval for the pooled effect size estimate
\item statistic: the test statistic for the pooled effect size estimate
\item p.value: the test statistic for the pooled effect size estimate
}

If \code{summary = FALSE} the function returns a tibble containing the models. Each row represents a category of the predictor, and columns contain
\itemize{
\item \code{name_of_predictor}: This column will be named after the predictor entered into the function and will list the categories (factor levels) in each row.
\item data: contains the raw data on which the model is based (that is, the data for the particular category)
\item model: each row contains the \code{rma.mv} object for the model fitted within each category
\item coefs: each row contains a table of coefficients for the model fitted within each category
}
}
\description{
\code{get_mas()} fits individual meta-analyses for each level of a categorical variable and (optionally) collates the results into a tabulated form for printing in a quarto document.
Models are fitted using using the \code{rma.mv()} function from the \code{metafor} package.
}
