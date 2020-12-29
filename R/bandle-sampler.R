diffLoc <- function(objectCond1,
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
                    propSd = c(0.3, 0.1, 0.05)){
    
        # Setting seed manually
    if (is.null(seed)) {
        seed <- sample(.Machine$integer.max, 1)
    }
    .seed <- as.integer(seed)  
    set.seed(.seed)
    
    
    # Elementary checks
    stopifnot(length(hyperMean)==3)
    stopifnot(length(hyperSd)==3)
    stopifnot(sum(hyperSd > 0)==3)
    
    ## Samples to be retained as a result of burnin and thinning
    toRetain <- seq.int(burnin + 1L, numIter, thin)
    numRetained <- length(toRetain)
    
    # Normalising data
    normCond1 <- lapply(objectCond1,  function(x) normalise(x, method = "center.mean"))
    normCond2 <- lapply(objectCond2,  function(x) normalise(x, method = "center.mean"))
    
    # Bring data together
    object_cmb <- c(cond1 = normCond1, cond2 = normCond2)
    
    # numCond should be 2
    numCond <- 2
    numRepl <- length(objectCond1)
    
    # expressions from data and important summaries
    exprs_cmb <- lapply(object_cmb, exprs)
    M <- lapply(exprs_cmb, colMeans)
    V <- lapply(exprs_cmb, function(x) lambda * diag(diag(cov(x))) + diag(rep(10^{-6}, ncol(x))))
    
    # Getting data dimensions
    D <- ncol(object_cmb[[1]])
    K <- length(getMarkerClasses(object_cmb[[1]]))
    
    # construct empirical Bayes Polya-Gamma prior
    if (is.null(pgPrior)) {
        pgPrior <- pg_prior(normCond1, normCond2, K = K, pgPrior = NULL)
    }
    
    # Fit GPs to markers
    componenthypers <- lapply(object_cmb, function(x) vector(mode = "list", K))
    if (maternCov == FALSE) {
        res <- lapply(object_cmb, function(x) fitGP(x, fcol = fcol))
    } else if ((maternCov == TRUE) & (PC = FALSE) & (is.null(gpParams))) {
        res <- lapply(object_cmb, function(x) fitGPmatern(x, fcol = fcol, nu = nu))
    } else if ((maternCov ==TRUE) & (PC = TRUE) & (is.null(gpParams))) {
        res <- lapply(object_cmb, function(x) fitGPmaternPC(x, fcol = fcol, nu = nu, hyppar = pc_prior))  
    } else {
        res <- gpParams
    }
    
    # separate data into known and unknown
    unknown_cmb <- lapply(object_cmb, unknownMSnSet)
    exprsUnknown_cmb <- lapply(unknown_cmb, exprs)
    exprsKnown <- lapply(object_cmb, function(x) exprs(markerMSnSet(x, fcol = fcol)))
    
    # separate conditions allocations
    allocKnown <- lapply(object_cmb[c(1, numRepl + 1)], function(x) seq.int(K)[fData(markerMSnSet(x, fcol = fcol))[, fcol]])
    numProtein <- lapply(unknown_cmb[c(1, numRepl + 1)], nrow)
    
    # Some storage
    alloc <- lapply(numProtein, function(x) matrix(0, nrow = x, ncol = numRetained))
    allocOut <- lapply(numProtein, function(x) matrix(0, nrow = x, ncol = numRetained))
    outlierprob <- lapply(numProtein, function(x) matrix(0, nrow = x, ncol = numRetained))
    weights <- array(0, c(K, K, numRetained))
    epsilons <- matrix(0, nrow = 2, ncol = numRetained)
    allocprob <- lapply(numProtein, function(x) array(0, c(x, numRetained, K)))
    loglikelihoods <- lapply(rep(numProtein, numRepl), function(x) matrix(0, nrow = x, ncol = K))
    
    #random allocation of unknown Proteins, allocations are condition specific
    alloctemp <- lapply(numProtein, function(x) sample.int(n = K, size = x, replace = TRUE))
    for (i in seq.int(numCond)) {
        
        object_cmb[[i]] <- knnClassification(object_cmb[[i]], k = 10)
        alloc[[i]][, 1] <- fData(object_cmb[[i]][rownames(unknown_cmb[[i]])])$knn
    }
    
    # number of proteins allocated to each component
    nkknown <- lapply(object_cmb[c(1, numRepl + 1)], function(x)
        table(getMarkers(x, verbose = FALSE))[getMarkerClasses(x)])
    
    #number initial allocated to outlier component
    outlier <- vector(mode = "list", length = 2)
    for (i in seq.int(numCond)) {
        outlier[[i]] <- allocOut[[i]][,1]
    }
    
    sampleGPMean <- lapply(object_cmb, function(x) vector(mode = "list", length = K))
    centereddata <- lapply(object_cmb, function(x) vector(mode = "list", length = K))
    hypers <- lapply(object_cmb, function(x) vector(mode = "list", length = numRetained))
    .hypers <- vector(mode = "list", length = length(object_cmb))
    
    for (i in seq.int(object_cmb)) {
        hypers[[i]][[1]] <- res[[i]]$params
        .hypers[[i]] <- hypers[[i]][[1]]
    }
    # intialise polya-gamma auxiliary variables
    w <- rep(1, K^2)
    
    for(t in 2:numIter){
        
        # Between data allocation tally
        nk_mat <- diag(nkknown[[1]]) + table(factor(alloctemp[[1]], levels = 1:K),
                                             factor(alloctemp[[2]], levels = 1:K))
        
        # Within data allocation tally
        nk <- list(cond1 = rowSums(nk_mat), cond2 = colSums(nk_mat))
        
        if((t %% 1) ==0){
            print(t)
        }
        
        gpnk_cond1 <- gpnk_cond2 <- vector(mode = "numeric", length = K)
        gpnk <- list(gpnk_cond1 = gpnk_cond1, gpnk_cond2 = gpnk_cond2)
        
        for (i in seq.int(object_cmb)) {
            
            j <- ceiling(i/numRepl) # allocs don't change across replicates
            for (l in seq.int(K)) {
                gpnk[[j]][l] <- sum(alloctemp[[j]]*outlier[[j]] == l) + nkknown[[1]][l]
            }
            if (maternCov == FALSE) {
                centereddata[[i]] <- centeredData(Xknown = exprsKnown[[i]],
                                                  BX = allocKnown[[j]],
                                                  Xunknown = exprsUnknown_cmb[[i]], 
                                                  BXun = alloctemp[[j]]*outlier[[j]],
                                                  hypers = .hypers[[i]], 
                                                  nk = gpnk[[j]], tau = seq.int(D), D = D, K)
            } else if (maternCov == TRUE) {  
                centereddata[[i]] <- centeredDatamatern(Xknown = exprsKnown[[i]],
                                                        BX = allocKnown[[j]],
                                                        Xunknown = exprsUnknown_cmb[[i]], 
                                                        BXun = alloctemp[[j]]*outlier[[j]],
                                                        hypers = .hypers[[i]], 
                                                        nk = gpnk[[j]], tau = seq.int(D), D = D, K, nu = nu)
            }
        }
        
        # Polya-Gamma Sampler or Dirichlet weights
        if (pg == TRUE) {
            pgres <- sample_weights_pg(nk_mat = nk_mat, pgPrior = pgPrior, K = K, w = w, tau = tau)
            currentweights <- pgres$currentweights
            w <- pgres$w
        } else {
            if (is.null(dirPrior)){
                dir_prior <- diag(rep(1, K)) + matrix(0.05, nrow = K, ncol = K)
            }
            currentweights <- sample_weights_dir(nk_mat = nk_mat, dir_prior = dir_prior)
        }
        
        currentweights <- matrix(currentweights, ncol = K, nrow = K)

        #extract noise component from hypers
        sigmak <- lapply(hypers, function(x) exp(2 * x[[1]][, 3])) # fixed for the moment
        sigmasqrt <- lapply(sigmak, sqrt)
        loglikelihoods <- comploglikelist(centereddata, sigmasqrt)
        
        resAllocCond1 <- proteinAllocation(loglikelihoods = loglikelihoods[1:numRepl],
                                           currentweights = currentweights, alloctemp = alloctemp, cond = 1)
        resAllocCond2 <- proteinAllocation(loglikelihoods = loglikelihoods[(1+numRepl):(numRepl*numCond)],
                                           currentweights = currentweights, alloctemp = alloctemp, cond = 2)
        
        alloctemp <- list(resAllocCond1$alloctemp, resAllocCond2$alloctemp)
        loglikelihoods_cond1 <- resAllocCond1$loglikelihoods_comb
        loglikelihoods_cond2 <- resAllocCond2$loglikelihoods_comb
        
        # sample epsilon
        tau1 <- lapply(outlier, function(x) sum(x == 1) + nrow(exprsKnown[[1]]))
        tau2 <- lapply(outlier, function(x) sum(x == 0))
        epsilon <- rbeta(n = 2, shape1 = u + unlist(tau2), shape2 = v + unlist(tau1))
        
        
        allocnotOutprob <- vector(mode = "list")
        allocOutprob <- vector(mode = "list")
        # Sample outlier allocations first condition 1 then condition 2
        # outlier likelihood
        outlierlikelihood <- vector(mode = "list")
        for (i in seq.int(exprsUnknown_cmb)) {
            outlierlikelihood[[i]] <- dmvtCpp(exprsUnknown_cmb[[i]], mu_ = M[[i]], sigma_ = V[[i]], df_ = 4, log_ = TRUE, isChol_ = F)
        }
        
        resOutCond1 <- outlierAllocationProbs(outlierlikelihood = outlierlikelihood[1:numRepl],
                                              loglikelihoods = loglikelihoods_cond1, epsilon = epsilon, alloctemp = alloctemp, cond = 1)
        resOutCond2 <- outlierAllocationProbs(outlierlikelihood = outlierlikelihood[(1+numRepl):(numRepl*numCond)],
                                              loglikelihoods = loglikelihoods_cond2, epsilon = epsilon, alloctemp = alloctemp, cond = 2)
        
        allocnotOutprob[[1]] <- resOutCond1[[1]]
        allocOutprob[[1]] <- resOutCond1[[2]]
        
        allocnotOutprob[[2]] <- resOutCond2[[1]]
        allocOutprob[[2]] <- resOutCond2[[2]]
        
        
        # Sample outliers
        allocoutlierprob <- cbind(allocnotOutprob[[1]], allocOutprob[[1]])
        outlier[[1]] <- sampleOutlier(allocoutlierprob = allocoutlierprob)
        
        allocoutlierprob <- cbind(allocnotOutprob[[2]], allocOutprob[[2]])
        outlier[[2]] <- sampleOutlier(allocoutlierprob = allocoutlierprob)
        
        #update hypers
        if((t %% hyperIter) == 0){
            if (maternCov == FALSE) {  
                if(hyperLearn == "LBFGS"){
                }else if(hyperLearn == "MH"){
                    # Between data allocation tally
                    nk_mat <- diag(nkknown[[1]]) + table(factor(alloctemp[[1]], levels = 1:K),
                                                         factor(alloctemp[[2]], levels = 1:K))
                    
                    # Within data allocation tally
                    nk <- list(cond1 = rowSums(nk_mat), cond2 = colSums(nk_mat)) 
                    
                    for (i in seq.int(object_cmb)) {
                        for(j in seq.int(K)){
                            
                            l <- ceiling(i/numRepl) # allocs don't change across replicates
                            Y <- makeComponent(X = exprsKnown[[i]],
                                               BX = allocKnown[[l]],
                                               Y = exprsUnknown_cmb[[i]],
                                               BY = alloctemp[[l]] * outlier[[l]],
                                               j = j)
                            componenthypers[[i]][[j]] <- metropolisGP(inith = .hypers[[i]][j,],
                                                                      X = t(Y),
                                                                      tau = seq.int(D), nk = nk[[l]][j],
                                                                      D = D,
                                                                      niter = 1,
                                                                      hypMean = hypMean,
                                                                      hypSd = hypSd)$h[ ,2]
                            
                            
                        }
                        # stores current hyperparameters invisibily
                        .hypers[[i]] <- matrix(unlist(componenthypers[[i]]), ncol = 3, byrow = TRUE)
                    }  
                }
            } else if ((maternCov == TRUE) & (hyperLearn == "MH")) {
                # Between data allocation tally
                nk_mat <- diag(nkknown[[1]]) + table(factor(alloctemp[[1]], levels = 1:K),
                                                     factor(alloctemp[[2]], levels = 1:K))
                
                for (i in seq.int(object_cmb)) {
                    for(j in seq.int(K)) {
                        l <- ceiling(i/numRepl) # allocs don't change across replicates
                        Y <- makeComponent(X = exprsKnown[[i]],
                                           BX = allocKnown[[l]],
                                           Y = exprsUnknown_cmb[[i]],
                                           BY = - alloctemp[[l]] * 0,
                                           j = j)
                        
                        hyperstoupdate <- c(exp(.hypers[[i]][j,1:2]), exp(2 * .hypers[[i]][j,3]))
                        
                        .pc_prior <- pc_prior[j, ]
                        newhypers <- metropolisGPmatern(inith = hyperstoupdate,
                                                        X = t(Y),
                                                        tau = seq.int(D),
                                                        nk = nk[[l]][j],
                                                        D = D,
                                                        niter = 1,
                                                        nu = nu,
                                                        hyppar = .pc_prior,
                                                        propsd = propsd)$h[ ,2]
                        
                        componenthypers[[i]][[j]] <- c(log(newhypers[1:2]), log(newhypers[3])/2)
                    }
                    # stores current hyperparameters invisibily
                    .hypers[[i]] <- matrix(unlist(componenthypers[[i]]), ncol = 3, byrow = TRUE)
                    
                }
            }
        }  
        
        ## Only store iterations that are going to be retained
        if(t %in% toRetain) {
            
            s <- which(t == toRetain) # index of variable to save
            for(i in seq.int(object_cmb)) {
                hypers[[i]][[s]] <- .hypers[[i]]
            }  
            weights[, , s] <- currentweights
            for (i in seq.int(numCond)) {
                alloc[[i]][, s] <- alloctemp[[i]]
                allocOut[[i]][, s] <- outlier[[i]]
                if (i == 1) {
                    allocprob[[i]][, s, ] <- resAllocCond1$allocprobtemp
                } else {
                    allocprob[[i]][, s, ] <- resAllocCond2$allocprobtemp
                }
            }  
            epsilons[, s] <- epsilon
        }  
    }
    
    nicheParam <- list()
    .niche <- alloc
    .nicheProb <- allocprob
    .outlier <- allocOut
    for (i in seq.int(numCond)){

        for (j in seq.int(numRepl)){
            if (i == 1){
                nicheParam[[j]] <- .nicheParam(dataset = "control",
                                               replicate = as.integer(j),
                                               K = K,
                                               D = D, 
                                               method = "bandle",
                                               params = hypers[[j]]) 
            } else{
                nicheParam[[numRepl + j]] <- .nicheParam(dataset = "treatment",
                                                          replicate = j,
                                                          K = K,
                                                          D = D, 
                                                          method = "bandle",
                                                          params = hypers[[numRepl + j]])     
            }    
        }
    }
    
    # set correct dimension names
    dimnames(weights)[[1]] <- dimnames(weights)[[2]] <- getMarkerClasses(object = objectCond1[[1]],
                                                                         fcol = fcol)
    rownames(epsilons) <- c("Dataset 1", "Dataset 2")
    .niche <-  lapply(.niche, function(x){ rownames(x) <- rownames(objectCond1[[1]]); x})
    .nicheProb <- lapply(.nicheProb, function(x) {
                                        dimnames(x)[[1]] <- rownames(objectCond1[[1]])
                                        dimnames(x)[[3]] <- getMarkerClasses(objectCond1[[1]]); x})
    .outlier <- lapply(.outlier, function(x){ rownames(x) <- rownames(objectCond1[[1]]); x})
    
    ## construct bandleChains object
    .out <- .bandleChain(dataset = "bandleExperiment",
                         replicates = length(object_cmb),
                         n = length(toRetain),
                         K = K,
                         N = numProtein[[1]], 
                         weights = weights,
                         epsilons = epsilons,
                         niche = .niche,
                         nicheProb = .nicheProb,
                         outlier = .outlier,
                         nicheParams = .nicheParams(params = nicheParam))
    
    return(.out)
    
}
