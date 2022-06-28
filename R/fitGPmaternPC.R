##' The \code{fitGPmaternPC} function is a helper function to fit matern GPs to
##' data with penalised complexity priors on the hyperparameters.
##'
##' This set of functions allow users to fit GPs to their data. The
##' \code{fitGPmaternPC} function allows users to pass a vector of penalised
##' complexity hyperparameters using the \code{hyppar} argument. You must
##' provide a matrix with 3 columns and 1 row. The order of these 3 columns
##' represent the hyperparameters length-scale, amplitude, variance. We have
##' found that the `matrix(c(10, 60, 250), nrow = 1)` worked well for the
##' spatial proteomics datasets tested in Crook et al (2021). This was visually
##' assessed by passing these values and visualising the GP fit using the
##' \code{plotGPmatern} function (please see vignette for an example of the
##' output). Generally, (1) increasing the lengthscale parameter (the first
##' column of the \code{hyppar} matrix)  increases the spread of the covariance
##' i.e. the similarity between points, (2) increasing the amplitude parameter
##' (the second column of the \code{hyppar} matrix) increases the maximum value
##' of the covariance and lastly (3) decreasing the variance (third column of
##' the \code{hyppar} matrix) reduces the smoothness of the function to allow
##' for local variations. We strongly recommend users start with the recommended
##' parameters and change and assess them as necessary for their dataset by
##' visually evaluating the fit of the GPs using the \code{plotGPmatern}
##' function. Please see the vignettes for more details and examples.
##' 
##' @title Fit and plot matern GPs to spatial proteomics data.
##' @param object A instance of class `MSnSet`
##' @param fcol feature column to indicate markers. Default is "markers".
##' @param materncov `logical` indicating whether matern covariance is used,
##' else Gaussian covariance is used.
##' @param nu matern smoothness parameter. Default is 2.
##' @param hyppar The vector of penalised complexity hyperparameters, you must 
##' provide a matrix with 3 columns and 1 row. The order is hyperparameters
##' on length-scale, amplitude, variance.
##' @return Returns an object of class `gpParams` which stores the posterior
##' predictive means, standard deviations, variances and also the MAP 
##' hyperparamters for the GP. 
##' @md
##' @examples
##'## ====== fitGPmaternPC =====
##' library(pRolocdata)
##' data("tan2009r1")
##' set.seed(1)
##' tansim <- sim_dynamic(object = tan2009r1, 
##'                     numRep = 6L,
##'                    numDyn = 100L)
##' ## Please note that \code{hyppar} should be chosen carefully and tested
##' ## by checking the GP fit with the \code{\link{plotGPmatern}} function
##' ## (please see details above)
##' gpParams <- lapply(tansim$lopitrep, 
##' function(x) fitGPmaternPC(x, hyppar = matrix(c(10, 60, 100), nrow = 1)))
##' 
##' @rdname bandle-gpfit
fitGPmaternPC <- function(object = object,
                          fcol = "markers",
                          materncov = TRUE,
                          nu = 2,
                          hyppar = matrix(c(10, 60, 250), nrow = 1)) {
  
  stopifnot("object is not an instance of class MSnSet"=is(object, "MSnSet"))
  stopifnot("matercov must be a logical"=is(materncov, "logical"))
  stopifnot("hyppar must be a matrix"=is(hyppar, "matrix"))
  stopifnot("You must provide a matrix with 3 columns for hyperparmeters"
            =ncol(hyppar) == 3)
  if (!is.null(fcol) && !fcol %in% fvarLabels(object))
    stop("'", fcol, "' not found in feature variables.")
  
  ## storage
  componenthypers <- vector(mode = "list", 
                            length(getMarkerClasses(object, fcol = fcol)))
  
  ## size needed
  D <- ncol(object)
  K <- length(getMarkerClasses(object, fcol = fcol))
  
  # random grid sampling for starting values
  initialvalues <- seq(-5, 0, 0.5)
  init <- matrix(0, length(initialvalues), 3)
  for(i in seq_along(initialvalues)){
    init[i,] <- initialvalues[sample.int(length(initialvalues), 
                                         size = 3, replace = TRUE)]
  }
  
  # indexing sets
  idx <- seq.int(D)
  tau <- seq.int(D)
  
  # LBFGS routine to get hypers
  for (j in seq.int(K)) {
    
    exprs <- t(exprs(object[fData(object)[, fcol] == 
                              getMarkerClasses(object, fcol = fcol)[j], idx]))
    
    # optimisation step 
    res <- apply(init, 1, function(z){lbfgs(posteriorGPmatern,
                                            posteriorgradientGPmatern,
                                            vars = z,
                                            invisible = 1,
                                            epsilon = 1e-6,
                                            Xk = exprs,
                                            tau =  seq.int(D),
                                            nk = length(exprs)/D,
                                            D = D,
                                            materncov = materncov,
                                            nu = nu,
                                            hyppar = hyppar)})
    componenthypers[[j]] <- res[[which.min(lapply(res, function(x){max(x$value)}))]]$par
    
  }
  
  # put hypers here
  .hypers <- matrix(unlist(componenthypers), ncol = 3, byrow = TRUE)
  
  # extract important quantities
  rhomaternk <- exp(.hypers[,1])
  amaternk <- exp(.hypers[,2])
  sigma <- exp(2 * .hypers[,3])
  M <- vector(mode = "list", K)
  V <- vector(mode = "list", K)
  Var <- vector(mode = "list", K)
  
  # plotting routines
  for(j in seq.int(K)){
    Orgdata <- t(exprs(object[fData(object)$markers == 
                                getMarkerClasses(object, fcol = fcol)[j],idx]))
    nk <- table(fData(object)$markers)[getMarkerClasses(object, fcol = fcol)][j]
    S <- matrix(rep(seq.int(length(tau)), length(tau)), nrow = length(tau))
    params <- .hypers
    sigmak <- sigma[j]
    amatern <- amaternk[j]
    rhomatern <- rhomaternk[j]
    
    # trench computations needed
    covA <- matern(nu = nu, a = amatern, rho = rhomatern, tau = seq.int(D), D = D)
    R <- diag(1, D) + (nk * covA)/sigmak;
    trenchres <- trenchDetcpp(R[1,])
    Z <- trenchInvcpp(trenchres$v)
    invcov <- diag(1, nk*D)/sigmak - kronecker(matrix(1, nk, nk), Z %*% covA)/sigmak^2
    Kstar <- do.call(cbind, replicate(nk, covA, simplify = FALSE))
    Kstarstar <- rep(amatern^2 + sigmak, length(tau))
    M[[j]] <- Kstar %*% invcov %*% as.vector(Orgdata)
    V[[j]] <- as.matrix(sqrt(diag(diag(Kstarstar, length(tau)) - Kstar %*% invcov %*% t(Kstar))))
    Var[[j]] <- diag(rep(amatern^2, length(tau))) - Kstar %*% invcov %*% t(Kstar)
  }
  
  .res <- .gpParams(method = "fitGPmaternPC",
                    M = M, 
                    V = V, 
                    sigma = sigma, 
                    params = params)
  
  return(.res)
  
}

##' The \code{fitGPmatern} function fits matern GPs to data.
##' 
##' @title Fit matern GP to spatial proteomics data.
##' @param materncov `logical` indicating whether matern covariance is used.
##' @md
##' @examples 
##' ## ====== fitGPmatern =====
##' library(pRolocdata)
##' data("tan2009r1")
##' set.seed(1)
##' tansim <- sim_dynamic(object = tan2009r1, 
##'                     numRep = 6L,
##'                    numDyn = 100L)
##' gpParams <- lapply(tansim$lopitrep, function(x) fitGPmaternPC(x))
##' 
##' @rdname bandle-gpfit
fitGPmatern <- function(object = object,
                        fcol = "markers",
                        materncov = TRUE,
                        nu = 2) {
  
  stopifnot("object is not an instance of class MSnSet"=is(object, "MSnSet"))
  stopifnot("matercov must be a logical"=is(materncov, "logical"))
  if (!is.null(fcol) && !fcol %in% fvarLabels(object))
    stop("'", fcol, "' not found in feature variables.")
  
  ## storage
  componenthypers <- vector(mode = "list", 
                            length(getMarkerClasses(object, fcol = fcol)))
  
  ## dimensions needed
  D <- ncol(object)
  K <- length(getMarkerClasses(object, fcol = fcol))
  
  # random grid sampling for starting values
  initialvalues <- seq(-5, 0, 0.5)
  init <- matrix(0, length(initialvalues), 3)
  for(i in seq_along(initialvalues)){
    init[i, ] <- initialvalues[sample.int(length(initialvalues), 
                                          size = 3, replace = TRUE)]
  }
  
  # indexing sets
  idx <- seq.int(D)
  tau <- seq.int(D)
  
  # LBFGS routine to get hypers
  for (j in seq.int(K)) {
    
    exprs <- t(exprs(object[fData(object)[, fcol] == 
                              getMarkerClasses(object)[j], idx]))
    
    res <- apply(init, 1,function(z){lbfgs(likelihoodGPmatern,
                                           gradientGPmatern,
                                           vars = z,
                                           invisible = 1,
                                           epsilon = 1e-6,
                                           Xk = exprs,
                                           tau =  seq.int(D),
                                           nk = length(exprs)/D,
                                           D = D,
                                           materncov = materncov,
                                           nu = nu)})
    componenthypers[[j]] <- res[[which.min(lapply(res, function(x){max(x$value)}))]]$par
    
  }
  
  # get hypers
  .hypers <- matrix(unlist(componenthypers), ncol = 3, byrow = TRUE)
  rhomaternk <- exp(.hypers[,1])
  amaternk <- exp(.hypers[,2])
  sigma <- exp(2 * .hypers[,3])
  M <- vector(mode = "list", K)
  V <- vector(mode = "list", K)
  Var <- vector(mode = "list", K)
  
  # plotting
  for (j in seq.int(K)) {
    Orgdata <- t(exprs(object[fData(object)$markers == 
                                getMarkerClasses(object)[j],idx]))
    
    # require statistics
    nk <- table(fData(object)$markers)[getMarkerClasses(object)][j]
    S <- matrix(rep(seq.int(length(tau)), length(tau)), nrow = length(tau))
    params <- .hypers
    sigmak <- sigma[j]
    amatern <- amaternk[j]
    rhomatern <- rhomaternk[j]
    
    # trench compuations
    covA <- matern(nu = nu, a = amatern, rho = rhomatern, tau = seq.int(D), D = D)
    R <- diag(1, D) + (nk * covA)/sigmak;
    trenchres <- trenchDetcpp(R[1,])
    Z <- trenchInvcpp(trenchres$v)
    invcov <- diag(1, nk*D)/sigmak - kronecker(matrix(1, nk, nk), Z %*% covA)/sigmak^2
    Kstar <- do.call(cbind, replicate(nk, covA, simplify=FALSE))
    Kstarstar <- rep(amatern^2 + sigmak, length(tau))
    M[[j]] <- Kstar %*% invcov %*% as.vector(Orgdata)
    V[[j]] <- sqrt(diag(diag(Kstarstar, length(tau)) - Kstar %*% invcov %*% t(Kstar)))
    Var[[j]] <- diag(rep(amatern^2, length(tau))) - Kstar %*% invcov %*% t(Kstar)
  }
  
  ## output
  .res <- .gpParams(method = "fitGPmatern",
                    M = M, 
                    V = V,
                    sigma = sigma, 
                    params = params)
  
  return(.res)
  
}


##' The \code{plotGPmatern} function plots matern GPs
##'
##' @title Plot matern GPs to spatial proteomics data.
##' @param params The output of running `fitGPmatern`, `fitGPmaternPC` 
##' or `fitGP` which is of class `gpParams`
##' @param fcol feature column to indicate markers. Default is `"markers"`.
##' @md
##' @examples
##' ## ====== plotGPmatern =====
##' ## generate example data
##' library(pRolocdata)
##' data("tan2009r1")
##' set.seed(1)
##' tansim <- sim_dynamic(object = tan2009r1, 
##'                     numRep = 6L,
##'                    numDyn = 100L)
##' ## fit a GP
##' gpParams <- lapply(tansim$lopitrep, function(x) fitGP(x))
##' 
##' ## Overlay posterior predictives onto profiles
##' ## Dataset1 1
##' par(mfrow = c(2, 3))
##' plotGPmatern(tansim$lopitrep[[1]], gpParams[[1]])
##' 
##' ## Dataset 2, etc.
##' par(mfrow = c(2, 3))
##' plotGPmatern(tansim$lopitrep[[2]], gpParams[[2]])
##' @rdname bandle-gpfit
##' @return The functions `plotGPmatern` plot the posterior
##' predictives overlayed with the markers for each subcellular class.
plotGPmatern <- function(object = object,
                         params = params,
                         fcol = "markers") {  
  
  stopifnot("object is not an instance of class MSnSet"=is(object, "MSnSet"))
  stopifnot("params is not an instance of class gpParams"=is(params, "gpParams"))
  if (!is.null(fcol) && !fcol %in% fvarLabels(object))
    stop("'", fcol, "' not found in feature variables.")
  
  # ## size needed
  K <- length(getMarkerClasses(object, fcol = fcol))
  M <- params@M
  V <- params@V
  D <- ncol(object)
  
  # indexing sets
  idx <- seq.int(D)
  tau <- seq.int(D)
  
  # LBFGS routine to get hypers
  for (j in seq.int(K)) {
    
    exprs <- t(exprs(object[fData(object)[, fcol] == 
                              getMarkerClasses(object, fcol = fcol)[j], idx]))
  }
  
  # plotting routines
  for(j in seq.int(K)){
    Orgdata <- t(exprs(object[fData(object)$markers == 
                                getMarkerClasses(object, fcol = fcol)[j],idx]))
    matplot(x = idx, Orgdata, col = getStockcol()[j],
            pch = 19, type = "b", lty = 1, lwd = 1.5,
            main = paste(getMarkerClasses(object, fcol = fcol)[j]),
            xlab = "Fraction", ylab = "Normalised Abundance", cex.main = 2,
            ylim = c(min(Orgdata) - 0.05, max(Orgdata) + 0.05),
            cex.axis = 1.5, cex.main = 1.5,
            xaxt = "n", axes = FALSE)
    axis(2)
    axis(1, at = idx, labels = idx)
    points(seq_along(tau), M[[j]], col = "black", pch = 19, cex = 1.3,
           type = "b", lwd = 5, lty = 1)
    arrows(seq_along(tau),
           M[[j]]-1.96*V[[j]], seq_along(tau),
           M[[j]]+1.96*V[[j]], length=0.1,
           angle=90, code=3,
           col = "black", lwd = 3)
  }
  
}
