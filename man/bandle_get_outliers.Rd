% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bandle-plots.R
\name{bandle_get_outliers}
\alias{bandle_get_outliers}
\title{Number of outliers at each iteration of MCMC}
\usage{
bandle_get_outliers(params)
}
\arguments{
\item{x}{Object of class \code{bandleParams}}
}
\value{
A \code{list} of length \code{length(x)}.
}
\description{
Helper function to get the number of outliers at each MCMC
iteration for each chain
}
