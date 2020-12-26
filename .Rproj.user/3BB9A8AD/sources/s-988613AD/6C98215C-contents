
proteinAllocation <- function(loglikelihoods, currentweights, alloctemp, cond) {
    
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
    
    return(list(alloctemp = alloctemp, loglikelihoods_comb = loglikelihoods_comb, allocprobtemp = allocprobtemp))
}

outlierAllocationProbs <- function(outlierlikelihood, loglikelihoods, epsilon, alloctemp, cond) {
    
    subset <- cbind(1:length(alloctemp[[cond]]), alloctemp[[cond]])
    outlierlikelihood_comb <- Reduce("+", outlierlikelihood)
    allocnotOutprob <- log(1 - epsilon[cond]) + loglikelihoods[subset]
    allocOutprob <- log(epsilon[cond]) + outlierlikelihood_comb
    
    return(list(allocnotOutprob = allocnotOutprob, allocOutprob = allocOutprob))
    
}

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


covOrganelle <- function(mydata, fcol = "markers"){
    
    M <- matrix(NA, nrow = length(getMarkerClasses(mydata, fcol = fcol)), ncol = ncol(mydata))
    rownames(M) <- getMarkerClasses(mydata, fcol = fcol)
    for (j in getMarkerClasses(mydata, fcol = fcol)) {
        M[j, ] <- colMeans(exprs(mydata)[fData(mydata)[, fcol] == j,])
    }
    sigma <- cov(t(M))
    
    return(sigma)
}

pg_prior <- function(object_cond1, object_cond2, K, pgPrior = NULL) {
    
    if (is.null(pgPrior)) {
        mu_prior <- rep(-7, K^2)
        mu_prior[c(K+1 * seq.int(1:K) - K)] <- mu_prior[c(K+1 * seq.int(1:K) - K)] + 1
        sigma1 <- covOrganelle(mydata = object_cond1[[1]])
        sigma2 <- covOrganelle(mydata = object_cond2[[1]])
    }
    pgPrior <- list(mu_prior = mu_prior, sigma1 =  sigma1, sigma2 = sigma2)
    return(pgPrior)
}


sample_weights_pg <- function(nk_mat, pgPrior, w, K, tau = tau) {
    
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

sample_weights_dir <- function(nk_mat, dir_prior){
    
    #sample weights from dirichlet by normalising gammas
    K <- ncol(nk_mat)
    concentration <- dir_prior + nk_mat
    currentweights <- t(sampleDirichlet(K^2, c(concentration)))
    return(currentweights)
}