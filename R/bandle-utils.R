## diff loc helpers


probsameorganelle <- function(params) {
    
    res <- vector(mode = "numeric", length = nrow(unknownMSnSet(object)))
    for (i in 1:nrow(unknownMSnSet(object))) {
        res[i] <- sum(params$alloc[[1]][i, ] == params$alloc[[2]][i, ])
    }
    names(res) <- rownames(unknownMSnSet(object))
    
    return(res)
}

## Bootstrap Probabilities

bootstrapdiffloc <- function(object, params, top = 20, Bootsample = 5000, decreasing = FALSE) {
    
    res <- matrix(NA, ncol = Bootsample, nrow = top)
    probs <- probsameorganelle(object = object, params = params)
    
    prBootnames <- names(probs[order(probs, decreasing = decreasing)][seq.int(top)])
    rownames(params$alloc[[1]]) <- rownames(params$alloc[[2]]) <- rownames(unknownMSnSet(object))
    
    rownames(res) <- prBootnames
    for (t in seq.int(Bootsample)){
        bootidx <- sample.int(ncol(params$alloc[[1]]), replace = TRUE)
        for (i in prBootnames) {
            res[i,t] <- sum(params$alloc[[1]][i, bootidx] == params$alloc[[2]][i, bootidx])
        }
    }
    
    return(res)
}

## Compute means of each organelle just using markers
meanOrganelle <- function(mydata, fcol = "markers"){
    
    M <- V <- matrix(NA, nrow = length(getMarkerClasses(mydata, fcol = fcol)), ncol = ncol(mydata))
    rownames(M) <- rownames(V) <-  getMarkerClasses(mydata, fcol = fcol)
    for (j in getMarkerClasses(mydata, fcol = fcol)) {
        M[j, ] <- colMeans(exprs(mydata)[fData(mydata)[, fcol] == j,])
        V[j, ] <- apply(exprs(mydata)[fData(mydata)[, fcol] == j,], 2, var)
    }
    
    return(list(M = M, V = V))
}

kldirpg <- function(sigma, mu, alpha) {
    
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

kldir <- function(alpha, beta) {
    
    res <- lgamma(sum(alpha)) - sum(lgamma(alpha)) - lgamma(sum(beta)) + sum(lgamma(beta)) +
        sum((alpha - beta) * (digamma(alpha) - digamma(sum(alpha))))
    
    
    return(res)
    
}



# prior predictive checks for dirichlet prior

prior_pred_dir <- function(object,
                           iter = 5000,
                           dir_prior = NULL,
                           q = 15) {
    
    K <- length(getMarkerClasses(object))
    if (is.null(dir_prior)) {
        dir_prior <- diag(rep(1, K)) + matrix(0.05, nrow = K, ncol = K)
    }
    
    priornotAlloc <- vector(length = iter)
    nkknown <- table(getMarkers(object, verbose = FALSE))[getMarkerClasses(object)]
    for (i in seq.int(iter)) {
        
        concentration <- diag(nkknown) + dir_prior
        currentweights <- t(sampleDirichlet(K^2, c(concentration)))
        priornotAlloc[i] <- sum(currentweights[1, -c((K + 1) * seq.int(1:(K)) - K)])
    }
    
    # average number of differential localisations
    meannotAlloc <- mean(priornotAlloc) * nrow(unknownMSnSet(object))
    
    # probability of having greater than q
    tailnotAlloc <- sum((priornotAlloc * nrow(unknownMSnSet(object))) > q)/iter
    
    return(list(priornotAlloc = priornotAlloc,
                meannotAlloc = meannotAlloc,
                tailnotAlloc = tailnotAlloc))
}

# prior predictive checks for Polya-Gamma


prior_pred_pg <- function(object,
                          tau = 0.2,
                          lambda = 0.01,
                          mu_prior = NULL,
                          iter = 10000,
                          q = 15) {
    
    # expressions from data and important summaries
    mydata <- exprs(object[[1]])
    M <- colMeans(mydata)
    V <- cov(mydata)
    nkknown <- table(getMarkers(object[[1]], verbose = FALSE))[getMarkerClasses(object[[1]])]
    K <- length(getMarkerClasses(object[[1]]))
    
    sigma1 <- covOrganelle(mydata = object[[1]])
    sigma2 <- covOrganelle(mydata = object[[4]])
    w <- rep(1, K^2)
    
    n_vec <- c(diag(nkknown))
    kappa <- n_vec - 1/2
    sigma_post <- solve(diag(w) + tau * kronecker(solve(t(sigma1) + diag(lambda, K)), solve(t(sigma2) + diag(lambda, K))))
    
    if (is.null(mu_prior)) {
        mu_prior <- rep(-9, K^2)
        mu_prior[c((K + 1) * seq.int(1:(K)) - K)] <- mu_prior[c((K + 1) * seq.int(1:(K)) - K)] + 1
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
        priornotAllocpg[i] <- sum(currentweights[-c((K + 1) * seq.int(1:(K)) - K)])
    }
    
    # average number of differential localisations
    meannotAlloc <- mean(priornotAllocpg) * nrow(unknownMSnSet(object[[1]]))
    varnotAlloc <- var(priornotAllocpg)
    
    # probability of having greater than 15 
    tailnotAlloc <- sum((priornotAllocpg * nrow(unknownMSnSet(object[[1]]))) > q)/iter
    
    
    return(list(priornotAllocpg = priornotAllocpg,
                meannotAlloc = meannotAlloc,
                varnotAlloc = varnotAlloc,
                tailnotAlloc = tailnotAlloc))
}


mcmc_plot_probs <- function(param, fname, n = 1, bw = 0.05, scale = "width", trim = TRUE) {
    Organelle <- Probability <- NULL
    stopifnot(length(fname) == 1)
    dfr <- as.data.frame(param[fname, , ])
    colnames(dfr) <- dimnames(param)[[3]]
    dfr_long <- data.frame(Organelle = rep(names(dfr), each = nrow(dfr)),
                           Probability = unlist(dfr, use.names = FALSE),
                           row.names = NULL,
                           stringsAsFactors = FALSE)
    gg2 <- ggplot(dfr_long,
                  aes(Organelle, Probability,
                      width = (Probability))) +
        geom_violin(aes(fill = Organelle), scale = scale, bw = bw)
    gg2 <- gg2 + theme_bw() +
        scale_fill_manual(values = pRoloc::getStockcol()[seq_len(nrow(dfr))]) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1),
              axis.title.x = element_blank())
    gg2 <- gg2 +
        ylab("Membership Probability") +
        ggtitle(paste0("Distribution of Subcellular Membership for Protein ", fname ))
    gg2 <- gg2 +
        theme(legend.position = "none")
    return(gg2)
}


require(akima)
require(fields)


spatial2D <- function(object,
                      dims = c(1, 2),
                      cov.function = wendland.cov,
                      theta = 2,
                      derivative = 2,
                      k = 1,
                      breaks = c(0.99, 0.95, 0.9, 0.85, 0.8, 0.75, 0.7),
                      aspect = 0.5) {
    
    # generate pca plot and create data from with probabilties
    .pca <- plot2D(object, dims = dims, plot = FALSE)
    probs <- data.frame(x = .pca[, 1], y = .pca[, 2], mcmc.prob = t(allocproborg)[rownames(object),])
    colnames(probs) <- c(c("x", "y"), getMarkerClasses(object))
    eigs <- colnames(.pca)
    
    # put data in appropriate long format
    probs.lst <- list()
    for(j in getMarkerClasses(object)) {
        probs.lst[[j]] <- probs[, c("x", "y", j)]
        colnames(probs.lst[[j]]) <- c("x", "y", "probability")
    }
    probs.lst.df <- plyr::ldply(probs.lst, .fun = function(x) x, .id = "organelle")
    
    # Create storage
    coords <- list()
    locations <- list()
    df <- list()
    # Create appropriate spatial grid
    for (j in getMarkerClasses(object)) {
        idxOrg <- c(probs.lst.df$organelle == j)
        coords[[j]] <- akima::interp(x = probs.lst.df$x[idxOrg],
                                     y = probs.lst.df$y[idxOrg],
                                     z = probs.lst.df$probability[idxOrg],
                                     extrap=FALSE, linear = TRUE, duplicate = TRUE) # interpolate onto appropriate grid
        coords[[j]]$z[is.na(coords[[j]]$z)] <- 0 # NaNs beyond data set to 0
        locations[[j]] <- cbind(rep(coords[[j]]$x, 40), rep(coords[[j]]$y, each = 40)) # get grid
        smoothedprobs <- fields::smooth.2d(coords[[j]]$z, x = locations[[j]],
                                           cov.function = cov.function,
                                           theta = theta,
                                           derivative = derivative, k = k) # spatial smoothing of probabilities
        # normalisation and formatting
        zmat <- matrix(smoothedprobs$z, ncol = 1)
        zmat <- zmat/max(zmat)
        df[[j]] <- data.frame(x = rep(smoothedprobs$x, 64), y = rep(smoothedprobs$y, each = 64), z = zmat)
    }
    # format data
    df.lst <- plyr::ldply(df, .fun = function(x) x, .id = "organelle") 
    df.lst <- df.lst %>%
        mutate(organelle = factor(organelle)) 
    K <- length(getMarkerClasses(object))
    cols <- getStockcol()[1:K] # get appropriate colours
    
    
    gg <- ggplot(
        data = df.lst,
        aes(x = x, y = y, z = z, color = organelle)) +
        coord_fixed() + 
        geom_contour(breaks = breaks, size = 1.2, aes(alpha = stat(level))) + 
        geom_point(alpha = 0) + 
        xlab(paste0(eigs[1])) + 
        ylab(paste0(eigs[2])) +
        scale_alpha(guide = "none") + 
        theme(legend.position = "right", 
              text = element_text(size = 12)) +
        scale_color_manual(values = cols) +
        scale_fill_manual(values = cols) +
        theme_minimal() + 
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              aspect.ratio = aspect,
              panel.border = element_rect(colour = "black", fill = NA, size = 1),
              plot.title = element_text(hjust = 0.5, size = 20),
              legend.text=element_text(size = 14)) +
        ggtitle(label = "Spatial variation of localisation probabilities") 
    gg
    return(gg)
}

