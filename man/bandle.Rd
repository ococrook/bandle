% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bandle-function.R, R/bandle-sampler.R,
%   R/bandle-utils.R
\name{bandle}
\alias{bandle}
\alias{bandlePredict}
\alias{bandleProcess}
\alias{diffLoc}
\alias{diffLocalisationProb}
\alias{bootstrapdiffLocprob}
\alias{binomialDiffLocProb}
\alias{meanOrganelle}
\alias{kldirpg}
\alias{kldir}
\alias{prior_pred_dir}
\alias{prior_pred_pg}
\title{Differential localisation experiments using the bandle method}
\usage{
bandle(
  objectCond1,
  objectCond2,
  fcol = "markers",
  hyperLearn = "LBFGS",
  numIter = 1000,
  burnin = 100L,
  thin = 5L,
  u = 2,
  v = 10,
  lambda = 1,
  gpParams = NULL,
  hyperIter = 20,
  hyperMean = c(0, 0, 0),
  hyperSd = c(1, 1, 1),
  seed = NULL,
  pg = FALSE,
  pgPrior = NULL,
  tau = 0.2,
  dirPrior = NULL,
  maternCov = TRUE,
  PC = TRUE,
  pcPrior = c(0.5, 3, 100),
  nu = 2,
  propSd = c(0.3, 0.1, 0.05),
  numChains = 4L,
  BPPARAM = BiocParallel::bpparam()
)

bandlePredict(objectCond1, objectCond2, params, fcol = "markers")

bandleProcess(params)

diffLoc(
  objectCond1,
  objectCond2,
  fcol = "markers",
  hyperLearn = "MH",
  numIter = 1000,
  burnin = 100L,
  thin = 5L,
  u = 2,
  v = 10,
  lambda = 1,
  gpParams = NULL,
  hyperIter = 20,
  hyperMean = c(0, 0, 0),
  hyperSd = c(1, 1, 1),
  seed = NULL,
  pg = TRUE,
  pgPrior = NULL,
  tau = 0.2,
  dirPrior = NULL,
  maternCov = TRUE,
  PC = TRUE,
  pcPrior = c(0.5, 3, 100),
  nu = 2,
  propSd = c(0.3, 0.1, 0.05)
)

diffLocalisationProb(params)

bootstrapdiffLocprob(params, top = 20, Bootsample = 5000, decreasing = TRUE)

binomialDiffLocProb(params, top = 20, nsample = 5000, decreasing = TRUE)

meanOrganelle(object, fcol = "markers")

kldirpg(sigma, mu, alpha)

kldir(alpha, beta)

prior_pred_dir(object, iter = 5000, dirPrior = NULL, q = 15)

prior_pred_pg(
  objectCond1,
  objectCond2,
  tau = 0.2,
  lambda = 0.01,
  mu_prior = NULL,
  iter = 10000,
  q = 15
)
}
\arguments{
\item{objectCond1}{An instance of class \code{MSnSet}, usually the control dataset}

\item{objectCond2}{An instance of class \code{MSnSet}, usually the treatment dataset}

\item{fcol}{a feature column indicating which feature define the markers}

\item{hyperLearn}{Algorithm to learn posterior hyperparameters of the Gaussian
processes. Default is \code{LBFGS} and \code{MH} for metropolis-hastings is also implemented.}

\item{numIter}{The number of iterations of the MCMC
algorithm. Default is 1000. Though usually much larger numbers are used}

\item{burnin}{The number of samples to be discarded from the
begining of the chain. Default is 100.}

\item{thin}{The thinning frequency to be applied to the MCMC
chain.  Default is 5.}

\item{u}{The prior shape parameter for Beta(u, v). Default is 2}

\item{v}{The prior shape parameter for Beta(u, v). Default is 10.}

\item{lambda}{The \code{lambda} ridge parameter used for numerical stability.
Default is 0.01}

\item{gpParams}{Object of class gpparams. parameters from prior fitting of GPs
to each niche to accelerate inference. Default is NULL.}

\item{hyperIter}{The frequency of MCMC interation to update the hyper-parameters
default is 20}

\item{hyperMean}{The prior mean of the log normal prior of the GP parameters.
Default is 0 for each. Order is length-scale, amplitude and noise variance}

\item{hyperSd}{The prior standard deviation of the log normal prior fo the GP
parameters. Default is 1 for each. Order is length-scale, ampliture and noise
variance.}

\item{seed}{The random number seed.}

\item{pg}{\code{logical} indicating whether to use polya-gamma prior. Default is
\code{FALSE}.}

\item{pgPrior}{A matrix generated by pgPrior function. If param pg is TRUE
but pgPrior is NULL then a pgPrior is generated on the fly.}

\item{tau}{The \code{tau} parameter of the Polya-Gamma prior. Default is 0.2.}

\item{dirPrior}{The Dirichlet prior used. If NULL (default) will generate a
a default Dirichlet prior}

\item{maternCov}{\code{logical} indicated whether to use a matern or gaussian
covariance. Default is True and matern covariance is used}

\item{PC}{\code{logical} indicating whether to use a penalised complexity prior.
Default is TRUE.}

\item{pcPrior}{\code{numeric} of length 3 indicating the lambda paramters for the
penalised complexity prior. Default is \code{c(0.5, 3, 100)} and the order is
length-scale, amplitude and variance.}

\item{nu}{\code{integter} indicating the smoothness of the matern prior. Default
is 2.}

\item{propSd}{If MH is used to learn posterior hyperparameters then the proposal
standard deviations. A Gaussian random-walk proposal is used.}

\item{numChains}{\code{integer} indicating the number of parallel chains to run.
Defaults to 4.}

\item{BPPARAM}{BiocParallel parameter. Defaults to machine default backend
using bpparam()}

\item{params}{An instance of \code{bandleParams}}

\item{top}{The number of proteins for which to sample from the binomial distribution}

\item{Bootsample}{Number of Bootstramp samples. Default is 5000}

\item{decreasing}{Starting at protein most likely to be differentially localization}

\item{nsample}{how many samples to return from the binomial distribution}

\item{object}{An instance of class \code{MSnSet}}

\item{sigma}{the sigma parameter of the Polya-Gamma prior}

\item{mu}{the mu parameter of the Polya-Gamma prior}

\item{alpha}{The concentration parameter of the first Dirichlet distribution}

\item{beta}{The concentration parameter of the second Dirichlet distribution}

\item{iter}{Number of sample to use from prior predictive distribution.
Default is 10000}

\item{q}{The upper tail value. That is the prior probability of having more
than q differential localisations. Default is 15.}

\item{mu_prior}{The mean of the Polya-Gamma prior. Default is NULL which generates
a default Polya-Gamma prior.}
}
\value{
\code{bandle} returns an instance of class \code{bandleParams}

\code{bandlePredict} returns an instance of class
\code{\link[MSnbase:MSnSet-class]{MSnbase::MSnSet}} containing the localisation predictions as
a new \code{bandle.allocation} feature variable. The allocation
probability is encoded as \code{bandle.probability}
(corresponding to the mean of the distribution
probability). In addition the upper and lower quantiles of
the allocation probability distribution are available as
\code{bandle.probability.lowerquantile} and
\code{bandle.probability.upperquantile} feature variables. The
Shannon entropy is available in the \code{bandle.mean.shannon}
feature variable, measuring the uncertainty in the allocations
(a high value representing high uncertainty; the highest value
is the natural logarithm of the number of classes). An additional variable
indicating the differential localization probability is also added as
\code{bandle.differential.localisation}

\code{bandleProcess} returns an instance of class
\code{bandleParams} with its summary slot populated.

\code{bandle} returns an instance of class \code{bandleParams}

returns a named vector of differential localisation probabilities

returns a matrix of size Bootsample * top containing bootstrap

returns a list containing empirical binomial samples

returns a list of means and variances for each

returns a numeric indicating the KL divergence

a numeric indicating the KL divergence

A list contain the prior predictive distribution of
differential localisations, the mean number of differential localised proteins
and the probability than more than q are differentially localised

A list contain the prior predictive distribution of
differential localisations, the mean number of differential localised proteins
and the probability than more than q are differentially localised
}
\description{
These function implement the bandle model for dynamic mass spectrometry based
spatial proteomics datasets using MCMC for inference

These functions implement the bandle model for dynamic mass spectrometry based
spatial proteomics datasets using MCMC for inference, this is an internal sampling
function

These functions implement helper functions for the bandle method
}
\details{
The \code{bandle} function generate the sample from the posterior distributions
(object or class \code{bandleParams}) based on an annotated quantitative spatial
proteomics datasets (object of class \code{\link[MSnbase:MSnSet-class]{MSnbase::MSnSet}}). Both are then
passed to the \code{bandlePredict} function to predict the sub-cellular localisation
and compute the differential localisation probability of proteins. See the
vignette for examples

The \code{diffloc} function generate the sample from the posterior distributions
(object or class \code{bandleParam}) based on an annotated quantitative spatial
proteomics datasets (object of class \code{\link[MSnbase:MSnSet-class]{MSnbase::MSnSet}}). Both are then
passed to the \code{bandlePredict} function to predict the sub-cellular localisation
and compute the differential localisation probability of proteins. See the
vignette for examples
}
