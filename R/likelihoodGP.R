##' Internal R function to pass R to C++, not for external use.
##' 
##' 
##' @title Compute GP likelihood squared exponential kernel
##' @param Xk The data
##' @param tau The indexing parameters
##' @param h GP hyperparameters
##' @param nk Number of observations
##' @param D number of samples 
##' @return Returns gp negative log likelihood
##' @md
##' 
##' @rdname bandle-gp
likelihoodGP <- function(Xk, tau, h, nk, D){
  
  negloglike <- likelihoodGPcpp(Xk, tau, h, nk, D)
  
  return(negloglike)
}
##' @title Compute GP likelihood matern covariance
##' @param Xk The data
##' @param tau The indexing parameters
##' @param h GP hyperparameters
##' @param nk Number of observations
##' @param D number of samples 
##' @param materncov `logical` indicating whether matern covariance is used
##' @param nu matern smoothness parameter.
##' @return Returns gp negative log likelihood
##' @md
##' 
##' @rdname bandle-gp
likelihoodGPmatern <- function(Xk, tau, h, nk, D, materncov, nu){
  
  negloglike <- likelihoodGPcpp(Xk, tau, h, nk, D, materncov, nu)
  
  return(negloglike)
}
##' @title Compute GP posterior matern covariance
##' @param Xk The data
##' @param tau The indexing parameters
##' @param h GP hyperparameters
##' @param nk Number of observations
##' @param D number of samples 
##' @param materncov `logical` indicating whether matern covariance is used
##' @param nu matern smoothness parameter.
##' @param hyppar prior hyperparameters of the penalised complexity prior.
##' @return Returns the negative log posterior of the GP
##' @md
##' 
##' @rdname bandle-gp
posteriorGPmatern <- function(Xk, tau, h, nk, D, materncov, nu, hyppar){
  
  negpost <- likelihoodGPcpp(Xk, tau, h, nk, D, materncov, nu) - PCrhomvar(rho = exp(h[1]), a = exp(h[2]), lambda1 = hyppar[1],
                                                                              lambda2 = hyppar[2], log = TRUE) - Gumbel(1/exp(2*h[3]), lambda = hyppar[3], log = TRUE)
  
  return(negpost)
}
##' @title Type-2 Gumbel distribution
##' @param x observation
##' @param lambda scale parameter of the type-2 Gumbel distribution
##' @param log `logical` indicating whether to return `log`. Default is `TRUE`
##' @return Returns the likelihood of the type-2 GUmbel distribution
##' @md
##' 
##' @rdname bandle-gp
Gumbel <- function(x,
                   lambda,
                   log = TRUE){
  
  if (log == FALSE) {
    gumbel <- (lambda/2) * x^{-3/2} * exp(-lambda * x^{-1/2})
  } else {
    gumbel <- log(lambda/2) - 3*log(x)/2 - lambda * x^{-1/2}
  }
  
  return(gumbel)
}
##' @title Bivariate penalized complexity prior for length-scale and amplitude
##' @param rho length-scale parameter
##' @param a amplitude
##' @param lambda1 first parameter of distribution
##' @param lamdba2 second parameter of distribution
##' @param log `logical` indicating whether to return `log`. Default is `TRUE`
##' @return Returns the likelihood of the bivariate penalised complexity prior
##' @md
##' 
##' @rdname bandle-gp
PCrhomvar <- function(rho,
                      a,
                      lambda1,
                      lambda2,
                      log = TRUE){
  
  if (log == FALSE) {
    res <- (lambda1 * lambda2/2) * rho^{- 3/2} * exp(-lambda1 * rho^{-1/2} - lambda2 * a)
  } else {
    res <- log(lambda1 * lambda2/2) - 3 * log(rho)/2 - lambda1 * rho^{-1/2} - lambda2 * a
  }
  
  return(res) 
  
}