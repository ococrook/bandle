##' Internal sampling function, not for outside use documented for completness
##' 
##' 
##' @title sample allocations, probabilities and compute loglikilihoods
##' @param loglikelihoods the log likelihoods
##' @param currentweights the current allocations weights
##' @param alloctemp the current protein allocations
##' @param cond the control = 1, treatment = 2
##' @return returns samples for protein allocations, log likelihoods and probabilities
##' @md
##' 
##' @rdname bandle-internal
proteinAllocation <- function(loglikelihoods,
                              currentweights,
                              alloctemp,
                              cond) {
    
    K <- ncol(currentweights)
    loglikelihoods_comb <- Reduce("+", loglikelihoods)
    if (cond == 1) {
        logconditionalprior <- t(log(currentweights[, alloctemp[[2]]])) # columns given condition 2
    } else {
        logconditionalprior <- log(currentweights[alloctemp[[1]], ]) # rows given condition 1
    }
    conditionalAlloc <- loglikelihoods_comb + logconditionalprior 
    cc <- apply(conditionalAlloc, 1, max)
    conditionalAlloc <- conditionalAlloc - cc # correct for underflow
    allocprobtemp <- exp(conditionalAlloc)/rowSums(exp(conditionalAlloc))
    #alloctemp <- apply(allocprobtemp, 1, function(x) sample.int(n = K, size = 1, replace = F, prob = x))
    alloctemp <- sampleAlloccpp(allocprobtemp)
    
    return(list(alloctemp = alloctemp,
                loglikelihoods_comb = loglikelihoods_comb,
                allocprobtemp = allocprobtemp))
}

##' @title computer outlier allocations probabilties
##' @param outlierlikelihood the outlier log likelihoods 
##' @param loglikelihoods the log likelihoods
##' @param epsilon the outlier component weight
##' @param alloctemp the current protein allocations
##' @param cond the control = 1, treatment = 2
##' @return returns outlier probabilities
##' @md
##' 
##' @rdname bandle-internal
outlierAllocationProbs <- function(outlierlikelihood,
                                   loglikelihoods,
                                   epsilon,
                                   alloctemp,
                                   cond) {
    
    subset <- cbind(seq.int(length(alloctemp[[cond]])), alloctemp[[cond]])
    outlierlikelihood_comb <- Reduce("+", outlierlikelihood)
    allocnotOutprob <- log(1 - epsilon[cond]) + loglikelihoods[subset]
    allocOutprob <- log(epsilon[cond]) + outlierlikelihood_comb
    
    return(list(allocnotOutprob = allocnotOutprob, allocOutprob = allocOutprob))
    
}
##' @title sample outlier probabilities
##' @param allocoutlierprob the outlier probabilities
##' @return returns outlier allocations
##' @md
##' 
##' @rdname bandle-internal
sampleOutlier <- function(allocoutlierprob){
    
    c <- apply(allocoutlierprob , 1, max) 
    allocoutlierprob <- allocoutlierprob - c # correct for underflow
    allocoutlierprob <- exp(allocoutlierprob)/rowSums(exp(allocoutlierprob))
    
    outlier <- sampleOutliercpp(allocoutlierprob) # reversed sample so 2nd entry is prob of 0
    # outlier <- apply(allocoutlierprob, 1, function(z){
    #   sample(x = c(1, 0), 1, prob = z)}
    # ) # reversed sample so 2nd entry is prob of 0
    
    return(outlier)
}
##' @title compute organelle covariances
##' @param object An instance of class `MSnSet`
##' @param fcol feature column indicating marker data. Default is "markers"
##' @return returns covariance of organelles using marker proteins
##' @md
##' 
##' @examples 
##' library(pRolocdata)
##' data("tan2009r1")
##' covOrganelle(object = tan2009r1)
##' 
##' 
##' @rdname bandle-internal
covOrganelle <- function(object, fcol = "markers"){
    
    stopifnot("object must be a class of MSnSet"=is(object, "MSnSet"))
    
    M <- matrix(NA, nrow = length(getMarkerClasses(object, fcol = fcol)), ncol = ncol(object))
    rownames(M) <- getMarkerClasses(object, fcol = fcol)
    for (j in getMarkerClasses(object, fcol = fcol)) {
        M[j, ] <- colMeans(Biobase::exprs(object)[fData(object)[, fcol] == j, ])
    }
    sigma <- cov(t(M))
    
    return(sigma)
}
##' @title Compute empirical Bayes Polya-Gamma prior
##' @param object_cond1 A list of instance of class `MSnSets` usually control
##' @param object_cond2 A list of instance of class `MSnSets` usually treatment
##' @param K the number of organelle classes
##' @param pgPrior The Polya-Gamma if user provided. Default is NULL to obtain value
##'  empirically
##' @param fcol The feature column containing the markers.
##' @return returns the Polya-Gamma prior
##' @md
##' 
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
##' out <- pg_prior(object_cond1 = control1,
##'  object_cond2 = treatment1, K = 11) 
##' 
##' 
##' @rdname bandle-internal
pg_prior <- function(object_cond1, 
                     object_cond2, 
                     K, 
                     pgPrior = NULL, 
                     fcol = "markers") {
    
    if (is.null(pgPrior)) {
        mu_prior <- rep(-7, K^2)
        mu_prior[c(K+1 * seq.int(K) - K)] <- mu_prior[c(K+1 * seq.int(K) - K)] + 1
        sigma1 <- covOrganelle(object = object_cond1[[1]], fcol = fcol)
        sigma2 <- covOrganelle(object = object_cond2[[1]], fcol = fcol)
    }
    pgPrior <- list(mu_prior = mu_prior, sigma1 =  sigma1, sigma2 = sigma2)
    return(pgPrior)
}
##' @title Sample mixture weights given the polya-gamma prior
##' @param nk_mat The summary matrix of allocations
##' @param pgPrior The Polya-Gamma prior
##' @param w The Polya-Gamma auxiliary variable 
##' @param K The number of organelle classes
##' @param tau The empirical bayes parameter for the Polya-Gamma variable. 
##'  Defaults to 0.2. 
##' @return returns A sample of the weights using Polya-Gamma priors. 
##' @md
##' 
##' @rdname bandle-internal
sample_weights_pg <- function(nk_mat,
                              pgPrior,
                              w,
                              K,
                              tau = 0.2) {
    
    # Polya-Gamma prior
    mu_prior <- pgPrior$mu_prior
    sigma1 <- pgPrior$sigma1
    sigma2 <- pgPrior$sigma2
    
    # Polya-Gamma Sampler
    n_vec <- c(nk_mat)
    kappa <- n_vec - 1/2
    sigma_post <- solve(diag(w) + tau * kronecker(solve(t(sigma1) + diag(0.01, K)), solve(t(sigma2) + diag(0.01, K))))
    mu_post <- sigma_post %*% (kappa + tau * kronecker(solve(t(sigma1) + diag(0.01, K)), solve(t(sigma2) + diag(0.01, K))) %*% mu_prior)
    phi <- mu_post + chol(sigma_post) %*% rnorm(n = K^2, mean = 0, sd = rep(1, K^2))
    
    #stick breaking construction
    currentweights <- rep(0, length(mu_post))
    stick <- 1
    for (j in seq.int(K^2 - 1)) {
        currentweights[j] <- (1/(1 + exp(-phi[j]))) * stick # Stick-breaking construction
        stick <- stick - currentweights[j]
    }
    currentweights[K^2] <- stick
    
    # Sample Polya-Gamma Variable
    w <- rcpp_pgdraw(rep(1, K^2), phi)
    
    return(list(currentweights = currentweights, w = w))
}
##' @title Sample mixture weights given the Dirichlet prior
##' @param nk_mat The summary matrix of allocations
##' @param dirPrior The Dirichlet prior
##' @return returns A sample of the weights using Dirichlet prior. 
##' @md
##' 
##' @rdname bandle-internal
sample_weights_dir <- function(nk_mat, dirPrior){
    
    #sample weights from dirichlet by normalising gammas
    K <- ncol(nk_mat)
    concentration <- dirPrior + nk_mat
    currentweights <- t(sampleDirichlet(K^2, c(concentration)))
    return(currentweights)
}