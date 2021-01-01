##' Internal R function to pass R to C++, not for external use.
##' 
##' 
##' @title Compute GP gradient
##' @param Xk The data
##' @param tau The indexing parameters
##' @param h GP hyperparameters
##' @param nk Number of observations
##' @param D number of samples 
##' @return Returns gp gradient
##' @md
##' 
##' @rdname bandle-gp
gradientGP <- function(Xk, tau, h, nk, D){
  
  grad <- gradientGPcpp(Xk, tau, h, nk, D)
  
  return(unlist(grad))
  
}
##' @title Compute GP gradient matern covariance
##' @param Xk The data
##' @param tau The indexing parameters
##' @param h GP hyperparameters
##' @param nk Number of observations
##' @param D number of samples 
##' @param materncov `logical` indicating whether matern covariance is used
##' @param nu matern smoothness parameter.
##' @return Returns gp gradient
##' @md
##' 
##' @rdname bandle-gp
gradientGPmatern <- function(Xk, tau, h, nk, D, materncov, nu){
  
  grad <- gradientGPcppmatern(Xk, tau, h, nk, D, nu)
  
  return(unlist(grad))
  
}
##' @title Compute GP posterior matern gradient
##' @param Xk The data
##' @param tau The indexing parameters
##' @param h GP hyperparameters
##' @param nk Number of observations
##' @param D number of samples 
##' @param materncov `logical` indicating whether matern covariance is used
##' @param nu matern smoothness parameter.
##' @param hyppar prior hyperparameters of the penalised complexity prior.
##' @return Returns the gradient of the posterior
##' @md
##' 
##' @rdname bandle-gp
posteriorgradientGPmatern <- function(Xk, tau, h, nk, D, materncov, nu, hyppar){
  
  grad <- gradientGPcppmatern(Xk, tau, h, nk, D, nu) + gradientlogprior(h, hyppar)
  
  return(unlist(grad))
  
}
##' @title  Compute the gradient of the log prior
##' @param  h numeric vector indicating value to evaluate
##' @param hyppar hyperaparameters of the prior
##' @return return the gradient of the log prior, length-scale, aamplitude and
##'  noise
##' @md
##' @rdname bandle-gp
gradientlogprior <- function(h, hyppar){
  
  gradrhoprior <- 3/2 - hyppar[1]*exp(- h[1]/2)/2 
  gradaprior <- hyppar[2]*exp(h[2])
  gradsigmaprior <- -3 + hyppar[3]*exp(h[3])
    
  gradprior <- as.matrix(c(gradrhoprior, gradaprior, gradsigmaprior))
  
  return(gradprior)
}
