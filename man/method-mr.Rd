% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mr-method.R
\name{robustMahalanobis}
\alias{robustMahalanobis}
\alias{reprodScore}
\alias{mrMethod}
\title{robust Mahalanobis distance}
\usage{
robustMahalanobis(delta)

reprodScore(x, y, method = c("pearson"))

mrMethod(objectCond1, objectCond2, plot = TRUE)
}
\arguments{
\item{delta}{The difference profile to compute the squared mahalanobis distance}

\item{x}{Numeric vector to compute reproducibility score}

\item{y}{Numeric vector to compute reprodducibility score}

\item{method}{Correlation method. Default is Pearson}

\item{objectCond1}{A list of \code{\link[MSnbase:MSnSet-class]{MSnbase::MSnSet}}s where each is an experimental
replicate for the first condition, usually a control}

\item{objectCond2}{A list of \code{\link[MSnbase:MSnSet-class]{MSnbase::MSnSet}}s where each is an experimental
replicate for the second condition, usually a treatment}

\item{plot}{Whether to generate an MR plot as a side effect.}
}
\value{
The squared Mahalanobis distance

The R score

The MR score of the Ithzak et al. 2016/2017
}
\description{
These function implement the MR method of Itzhak et al
}
