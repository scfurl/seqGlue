% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/feature_selection.R
\name{select_genes}
\alias{select_genes}
\title{Select features based on mean expression and variance}
\usage{
select_genes(
  obj,
  fit_min = 1,
  fit_max = Inf,
  logmean_ul = Inf,
  logmean_ll = -Inf,
  top_n = NULL
)
}
\arguments{
\item{obj}{the output of the \code{calc_dispersion} function.}

\item{fit_min}{the minimum multiple of the dispersion fit calculation; default = 1}

\item{fit_max}{the maximum multiple of the dispersion fit calculation; default = Inf}

\item{logmean_ul}{the maximum multiple of the dispersion fit calculation; default = Inf}

\item{logmean_ll}{the maximum multiple of the dispersion fit calculation; default = Inf}

\item{top}{top_n if specified, will override the fit_min and fit_max to select the top n most 
variant features.  logmena_ul and logmean_ll can still be used.}
}
\value{
an updated object (list) that records the selected features
}
\description{
Many genes typically co-vary with one another, and so the dimensionality of the
data can be reduced with a wide variety of different algorithms. After calculating dispersion for
a using the \code{calc_dispersion} function, the 
\code{select_genes} function allows the user to identify a set of genes
that will be used in downstream dimensionality reduction methods.
}
