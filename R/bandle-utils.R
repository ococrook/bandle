##' These functions implement helper functions for the bandle method
##' 
##' 
##' @title Compute differential localisation probabilities from ms-based
##' experiments using the bandle method
##' @param params An instance of class `bandleParams`
##' @return  returns a named vector of differential localisation probabilities
##' @md
##' @examples 
##' library(pRolocdata)
##' data("tan2009r1")
##' set.seed(1)
##' tansim <- sim_dynamic(object = tan2009r1, 
##'                     numRep = 6L,
##'                    numDyn = 100L)
##' gpParams <- lapply(tansim$lopitrep, function(x) 
##' fitGPmaternPC(x, hyppar = matrix(c(0.5, 1, 100), nrow = 1)))
##' d1 <- tansim$lopitrep
##' control1 <- d1[1:3]
##' treatment1 <- d1[4:6]
##' mcmc1 <- bandle(objectCond1 = control1, objectCond2 = treatment1, gpParams = gpParams,
##'                                      fcol = "markers", numIter = 10L, burnin = 1L, thin = 2L,
##'                                      numChains = 2, BPPARAM = SerialParam(RNGseed = 1))
##' mcmc1 <- bandleProcess(mcmc1)
##' dp <- diffLocalisationProb(mcmc1)
##' 
##' @rdname bandle
diffLocalisationProb <- function(params) {
    
    # Must be a valid bandleParams object
    stopifnot(is(params, "bandleParams"))
    
    res <- rowSums(1 * (params@chains@chains[[1]]@niche[[1]] - params@chains@chains[[1]]@niche[[2]]) != 0)
    res <- res/ncol(params@chains@chains[[1]]@niche[[1]])
    
    return(res)
}

##' @title Obtain bootstrap uncertainty estimates for differential localisations
##' probabilities
##' @param params An instance of class `bandleParams`
##' @param top The number of proteins for which to compute bootstrap distributions
##'  default is 20.
##' @param Bootsample Number of Bootstramp samples. Default is 5000
##' @param decreasing Starting at proteins most likely to be differentially localised.
##' 
##' @return  returns a matrix of size Bootsample * top containing bootstrap 
##' @md
##' @examples 
##' library(pRolocdata)
##' data("tan2009r1")
##' set.seed(1)
##' tansim <- sim_dynamic(object = tan2009r1, 
##'                     numRep = 6L,
##'                    numDyn = 100L)
##' gpParams <- lapply(tansim$lopitrep, 
##' function(x) fitGPmaternPC(x, hyppar = matrix(c(0.5, 1, 100), nrow = 1)))
##' d1 <- tansim$lopitrep
##' control1 <- d1[1:3]
##' treatment1 <- d1[4:6]
##' mcmc1 <- bandle(objectCond1 = control1, objectCond2 = treatment1, gpParams = gpParams,
##'                                      fcol = "markers", numIter = 10L, burnin = 1L, thin = 2L,
##'                                      numChains = 2, BPPARAM = SerialParam(RNGseed = 1))
##' mcmc1 <- bandleProcess(mcmc1)
##' bdp <- bootstrapdiffLocprob(mcmc1)
##' @rdname bandle
bootstrapdiffLocprob <- function(params,
                                 top = 20,
                                 Bootsample = 5000,
                                 decreasing = TRUE) {
    
    res <- matrix(NA, ncol = Bootsample, nrow = top)
    probs <- diffLocalisationProb(params = params)
    
    prBootnames <- names(probs[order(probs, decreasing = decreasing)][seq.int(top)])
    rownames(res) <- prBootnames
    
    for (t in seq.int(Bootsample)){
        bootidx <- sample.int(ncol(params@chains@chains[[1]]@niche[[1]]),
                              replace = TRUE)
        res[,t] <- rowSums(1 * (params@chains@chains[[1]]@niche[[1]][prBootnames, bootidx] - 
                                    params@chains@chains[[1]]@niche[[2]][prBootnames, bootidx] != 0))/length(bootidx)
    }
    
    return(res)
}
##' @title Obtain uncertainty estimates on differential localisation directly from binomial distributions,
##' using the Jeffries interval
##' @param params An instance of `bandleParams`
##' @param top The number of proteins for which to sample from the binomial distribution
##' @param nsample how many samples to return from the binomial distribution
##' @param decreasing Starting at protein most likely to be differentially localization
##' 
##' @return returns a list containing empirical binomial samples
##' @md
##' @examples 
##' library(pRolocdata)
##' data("tan2009r1")
##' set.seed(1)
##' tansim <- sim_dynamic(object = tan2009r1, 
##'                     numRep = 6L,
##'                    numDyn = 100L)
##' gpParams <- lapply(tansim$lopitrep, 
##' function(x) fitGPmaternPC(x, hyppar = matrix(c(0.5, 1, 100), nrow = 1)))
##' d1 <- tansim$lopitrep
##' control1 <- d1[1:3]
##' treatment1 <- d1[4:6]
##' mcmc1 <- bandle(objectCond1 = control1, objectCond2 = treatment1, gpParams = gpParams,
##'                                      fcol = "markers", numIter = 10L, burnin = 1L, thin = 2L,
##'                                      numChains = 2, BPPARAM = SerialParam(RNGseed = 1))
##' mcmc1 <- bandleProcess(mcmc1)
##' dp <- binomialDiffLocProb(mcmc1)
##' @rdname bandle
binomialDiffLocProb <- function(params,
                                top = 20,
                                nsample = 5000,
                                decreasing = TRUE){
    
    # Must be a valid bandleParams object
    stopifnot(is(params, "bandleParams"))
    
    res <- matrix(NA, ncol = nsample, nrow = top)
    probs <- diffLocalisationProb(params = params)
    
    prnames <- names(probs[order(probs, decreasing = decreasing)][seq.int(top)])
    rownames(res) <- prnames
    
    diff <- rowSums(1 * (params@chains@chains[[1]]@niche[[1]] - params@chains@chains[[1]]@niche[[2]]) != 0)
    nT <- ncol(params@chains@chains[[1]]@niche[[1]])
    
    #Jeffrey's samples
    res <- t(sapply(diff[prnames], function(x) rbeta(n = nsample, shape1 = x + 1/2, shape2 = nT - x + 1/2)))
    
    rownames(res) <- prnames
    
    return(res)
}

##' The EFDR for a given threshold is equal to the sum over all proteins
##' that exceed that threshold of one minus the posterior probability of
##' differential localisations, divides by the total number of proteins
##' with probabilities of differential localisation greater than that
##' threshold.
##' 
##' @title Compute the expected False Discovery Rate 
##' @param prob A numeric indicating probabilities of differential localisation
##' @param threshold A numeric indicating the probability threshold. The default
##' is 0.90.
##' @return The expected false discovery rate for a given threshold
##' @md
##' @examples 
##' library(pRolocdata)
##' data("tan2009r1")
##' set.seed(1)
##' tansim <- sim_dynamic(object = tan2009r1, 
##'                     numRep = 6L,
##'                    numDyn = 100L)
##' gpParams <- lapply(tansim$lopitrep, function(x) 
##' fitGPmaternPC(x, hyppar = matrix(c(0.5, 1, 100), nrow = 1)))
##' d1 <- tansim$lopitrep
##' control1 <- d1[1:3]
##' treatment1 <- d1[4:6]
##' mcmc1 <- bandle(objectCond1 = control1, objectCond2 = treatment1, gpParams = gpParams,
##'                                      fcol = "markers", numIter = 10L, burnin = 1L, thin = 2L,
##'                                      numChains = 2, BPPARAM = SerialParam(RNGseed = 1))
##' mcmc1 <- bandleProcess(mcmc1)
##' dp <- diffLocalisationProb(mcmc1)
##' EFDR(dp, threshold = 0.5)
##' 
##' @rdname bandle
EFDR <- function(prob, threshold = 0.90) {
    
    stopifnot("prob must be numeric"=is(prob, "numeric"))
    stopifnot("prob must be probabilities"=max(prob) <= 1)
    stopifnot("prob must be probabilities"=min(prob) >= 0)
    stopifnot("threshold must be numeric"=is(threshold, "numeric"))
    stopifnot("threshold must be a single values"=length(threshold) == 1)
    
    .out <- sum((1 - prob) * I(prob > threshold)) / sum(I(prob > threshold))
    return(.out)
}


##' @title Computes Organelle means and variances using markers
##' @param object a instance of class `MSnset`
##' @param fcol a feature column indicating which feature define the markers
##' @return returns a list of means and variances for each 
##' @md
##' @examples 
##' library(pRolocdata)
##' data("tan2009r1")
##' meanOrganelle(object = tan2009r1)
##' 
##' @rdname bandle
meanOrganelle <- function(object, fcol = "markers"){
    
    stopifnot(is(object, "MSnSet"))
    
    M <- V <- matrix(NA, nrow = length(getMarkerClasses(object, fcol = fcol)), ncol = ncol(object))
    rownames(M) <- rownames(V) <-  getMarkerClasses(object, fcol = fcol)
    for (j in getMarkerClasses(object, fcol = fcol)) {
        M[j, ] <- colMeans(Biobase::exprs(object)[fData(object)[, fcol] == j,])
        V[j, ] <- apply(Biobase::exprs(object)[fData(object)[, fcol] == j,], 2, var)
    }
    
    return(list(M = M, V = V))
}

##' @title Computes the Kullback-Leiber divergence between Polya-Gamma and 
##' Dirichlet priors
##' @param sigma the sigma parameter of the Polya-Gamma prior. A positive-definite
##' symmetric matrix.
##' @param mu the mu parameter of the Polya-Gamma prior. A vector of means
##' @param alpha the alpha (concentration) parameter of the Dirichlet prior
##' @return returns a numeric indicating the KL divergence
##' @md
##' @examples 
##' kldirpg(sigma = diag(c(1,1,1)), mu = c(0,0,0), alpha = 1)
##' 
##' @rdname bandle
kldirpg <- function(sigma = diag(1,1,1),
                    mu = c(0,0,0),
                    alpha = c(1)) {
    
    stopifnot("sigma must be matrix"=is(sigma, "matrix"))
    stopifnot("mu must be numeric"=is(mu, "numeric"))
    stopifnot("alpha must be numeric"=is(alpha, "numeric"))
    stopifnot("dimensions of sigma must match length of mu"=length(mu)==ncol(sigma))
    
    D <- length(mu)
    entpg <- -log((2*pi*exp(1))^D)/2 - sum(log(eigen(sigma)$values))/2 # compute determinant from eigenvalues
    entdir <- sum(lgamma(alpha)) - lgamma(sum(alpha))
    entdirpg_1 <- c(alpha) * (log(1 + exp(-mu)) + exp(mu) * diag(sigma)/(1 + exp(mu))^2)
    entdirpg_2 <- c(alpha) * (log(1 + exp(mu)) - exp(mu) * diag(sigma)/(1 + exp(mu))^2)
    entdirpgsum_1 <- sum(entdirpg_1)[-D]
    entdirpgsum_2 <- sum(seq.int(1, D - 1) * entdirpg_2[2:D])
    
    res <- entpg + entdir + entdirpgsum_1 + entdirpgsum_2
    
    return(res)
}

##' @title  Compute the KL divergence between two Dirichlet distributions
##' 
##' @param alpha The concentration parameter of the first Dirichlet distribution
##' @param beta The concentration parameter of the second Dirichlet distribution
##' @return a numeric indicating the KL divergence
##' @md
##' @examples 
##' kldir(c(1,1), c(3,1))
##' 
##' @rdname bandle
kldir <- function(alpha, beta) {
    
    stopifnot("alpha must be numeric"=is(alpha, "numeric"))
    stopifnot("beta must be numeric"=is(beta, "numeric"))
    
    res <- lgamma(sum(alpha)) - sum(lgamma(alpha)) - lgamma(sum(beta)) + sum(lgamma(beta)) +
        sum((alpha - beta) * (digamma(alpha) - digamma(sum(alpha))))
    
    
    return(res)
    
}

##' @title A function to compute the prior predictive distribution of the 
##' Dirichet prior.
##' 
##' @param object An instance of class `MSnSet`
##' @param fcol Feature column indicating the markers. Default is "markers"
##' @param iter Number of sample to use from prior predictive distribution.
##'  Default is 5000
##' @param dirPrior The Dirichlet prior used. If NULL (default) will generate a 
##'  a default Dirichlet prior. This should be a matrix with the same dimensions
##'  as the number of subcellular niches. The diagonal terms correspond
##'  to the prior probability of not differentially localising. The (i,j)
##'  term corresponds to prior probabilty of differntially localising between 
##'  niche i and j. 
##' @param q The upper tail value. That is the prior probability of having more 
##'  than q differential localisations. Default is 15. 
##' @return A list contain the prior predictive distribution of
##'   differential localisations, the mean number of differential localised proteins
##'   and the probability than more than q are differentially localised 
##' @md
##' 
##' @examples 
##' library(pRolocdata)
##' data("tan2009r1")
##' 
##' out <- prior_pred_dir(object = tan2009r1)
##' 
##' @rdname bandle     
prior_pred_dir <- function(object,
                           fcol = "markers",
                           iter = 5000,
                           dirPrior = NULL,
                           q = 15) {
    
    stopifnot("object must be of class MSnSet"=is(object, "MSnSet"))
    K <- length(getMarkerClasses(object, fcol = fcol))
    
    stopifnot("dirPrior must have dimensions equal to the number of
                  niches"=dim(dirPrior)==c(K, K))
    
    if (is.null(dirPrior)) {
        dirPrior <- diag(rep(1, K)) + matrix(0.05, nrow = K, ncol = K)
    }
    
    priornotAlloc <- vector(length = iter)
    nkknown <- table(getMarkers(object, verbose = FALSE, fcol = fcol))[getMarkerClasses(object, fcol = fcol)]
    for (i in seq.int(iter)) {
        
        concentration <- diag(nkknown) + dirPrior
        currentweights <- t(sampleDirichlet(K^2, c(concentration)))
        priornotAlloc[i] <- sum(currentweights[1, -c((K + 1) * seq.int(K) - K)])
    }
    
    # average number of differential localisations
    meannotAlloc <- mean(priornotAlloc) * nrow(unknownMSnSet(object, fcol = fcol))
    
    # probability of having greater than q
    tailnotAlloc <- sum((priornotAlloc * nrow(unknownMSnSet(object, fcol = fcol))) > q)/iter
    
    return(list(priornotAlloc = priornotAlloc,
                meannotAlloc = meannotAlloc,
                tailnotAlloc = tailnotAlloc))
}

##' @title A function to compute the prior predictive distribution of the 
##' Polya-Gamma prior.
##' 
##' @param objectCond1 An instance of class `MSnSet`, usually the control dataset
##' @param objectCond2 An instance of class `MSnSet`, usually the treatment dataset
##' @param fcol The feature column indiating the markers. Default is "markers"
##' @param tau The `tau` parameter of the Polya-Gamma prior. Default is 0.2.
##' @param lambda The `lambda` ridge parameter used for numerical stability. 
##'  Default is 0.01
##' @param mu_prior The mean of the Polya-Gamma prior. Default is NULL which generates
##'  a default Polya-Gamma prior.  
##' @param iter Number of sample to use from prior predictive distribution.
##'  Default is 10000
##' @param q The upper tail value. That is the prior probability of having more 
##'  than q differential localisations. Default is 15. 
##' @return A list contain the prior predictive distribution of
##'   differential localisations, the mean number of differential localised proteins
##'   and the probability than more than q are differentially localised 
##' @md
##' @examples 
##' library(pRolocdata)
##' data("tan2009r1")
##' set.seed(1)
##' tansim <- sim_dynamic(object = tan2009r1, 
##'                     numRep = 6L,
##'                    numDyn = 100L)
##' d1 <- tansim$lopitrep
##' control1 <- d1[1:3]
##' treatment1 <- d1[4:6]
##' out <- prior_pred_pg(objectCond1 = control1[[1]],
##' objectCond2 = treatment1[[1]])
##' 
##' 
##' @rdname bandle     
prior_pred_pg <- function(objectCond1,
                          objectCond2,
                          fcol = "markers",
                          tau = 0.2,
                          lambda = 0.01,
                          mu_prior = NULL,
                          iter = 10000,
                          q = 15) {
    
    stopifnot(is(objectCond1, "MSnSet"))
    stopifnot(is(objectCond2, "MSnSet"))
    
    # Expressions from data and important summaries
    mydata <- Biobase::exprs(objectCond1)
    M <- colMeans(mydata)
    V <- cov(mydata)
    nkknown <- table(getMarkers(objectCond1, verbose = FALSE, fcol = fcol))[getMarkerClasses(objectCond1, fcol = fcol)]
    K <- length(getMarkerClasses(objectCond1, fcol = fcol))
    
    sigma1 <- covOrganelle(object = objectCond1, fcol = fcol)
    sigma2 <- covOrganelle(object = objectCond2, fcol = fcol)
    w <- rep(1, K^2)
    
    n_vec <- c(diag(nkknown))
    kappa <- n_vec - 1/2
    sigma_post <- solve(diag(w) + tau * kronecker(solve(t(sigma1) + diag(lambda, K)), solve(t(sigma2) + diag(lambda, K))))
    
    if (is.null(mu_prior)) {
        mu_prior <- rep(-9, K^2)
        mu_prior[c((K + 1) * seq.int(K) - K)] <- mu_prior[c((K + 1) * seq.int(K) - K)] + 1
    }
    
    
    mu_post <- sigma_post %*% (kappa + tau * kronecker(solve(t(sigma1) + diag(lambda, K)),
                                                       solve(t(sigma2) + diag(lambda, K))) %*% mu_prior)
    
    priornotAllocpg <- vector(length = iter)
    for (i in seq.int(iter)) {
        phi <- mu_post + chol(sigma_post) %*% rnorm(n = K^2, mean = 0, sd = rep(1, K^2))
        
        currentweights <- rep(0, length(mu_post))
        
        stick <- 1
        for (j in seq.int(K^2 - 1)) {
            currentweights[j] <- (1/(1 + exp(-phi[j]))) * stick # Stick-breaking construction
            stick <- stick - currentweights[j]
        }
        currentweights[K^2] <- stick
        
        w <- rcpp_pgdraw(rep(1, K^2), phi)
        priornotAllocpg[i] <- sum(currentweights[-c((K + 1) * seq.int(K) - K)])
    }
    
    # average number of differential localisations
    meannotAlloc <- mean(priornotAllocpg) * nrow(unknownMSnSet(objectCond1, fcol = fcol))
    varnotAlloc <- var(priornotAllocpg)
    
    # probability of having greater than 15 
    tailnotAlloc <- sum((priornotAllocpg * nrow(unknownMSnSet(objectCond1, fcol = fcol))) > q)/iter
    
    
    return(list(priornotAllocpg = priornotAllocpg,
                meannotAlloc = meannotAlloc,
                varnotAlloc = varnotAlloc,
                tailnotAlloc = tailnotAlloc))
}
