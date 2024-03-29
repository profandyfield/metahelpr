% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/meta_helpers.R
\name{report_par_tbl}
\alias{report_par_tbl}
\title{Tabulate parameters from a meta-analysis object}
\usage{
report_par_tbl(metaobject, digits = 2, p_digits = 3, pred_labels)
}
\arguments{
\item{metaobject}{A meta-analyses object (\code{rma} or \code{rma.mv}) created using the \code{metafor} package}

\item{digits}{number of decimal places to print in the output}

\item{p_digits}{number of decimal places for \emph{p}-values}

\item{pred_labels}{A vector of text labels to use for each row of the table. For example, \code{pred_labels = c("Intercept", "FOA (Global) vs. FOA (Detailed)", "Not FOA (Disorganisation) vs. FOA (Detailed)", "Not FOA (Organisation) vs. FOA (Detailed)")}}
}
\value{
A tibble containing (for each predictor in the model) the name of the predictor (obtained from \code{pred_labels}), the effect size estimate, 95\% confidence interval, test statistic and \emph{p}-value.
}
\description{
\code{report_par_tbl()} is a helper function that outputs a tibble of the table of coefficients of a meta-analyses object (\code{rma} or \code{rma.mv})
created using the \code{metafor} package. This will mostly be useful for models containing predictors of effect sizes (so-called meta-regression).
}
