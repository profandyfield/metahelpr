% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/meta_helpers.R
\name{d_to_g}
\alias{d_to_g}
\title{Apply correction to Cohen's \emph{d} to make it the unbiased Hedges \eqn{g^\ast}}
\usage{
d_to_g(n, d)
}
\arguments{
\item{n}{the sample size on which \emph{d} is based}

\item{d}{the value of \emph{d}}
}
\value{
Hedges \eqn{g^\ast}
}
\description{
Although when applying \code{escalc} from \code{metafor} \emph{d} is converted to \eqn{g^\ast}, when the raw information is not available to compute \emph{d}
it may be necessary to find it indirectly and convert to Hedges statistic manually. \code{d_to_g()} does this adjustment.
}
\references{
Hedges, L. (1981). Distribution Theory for Glass’s Estimator of Effect Size and Related Estimators. \emph{Journal of Educational Statistics}, \emph{6}, 107–28.
}
