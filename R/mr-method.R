##' These function implement the MR method of Itzhak et al
##' 
##' @title robust Mahalanobis distance
##' @param delta The difference profile to compute the squared mahalanobis distance
##' @return The squared Mahalanobis distance 
##' @md
##' @rdname method-mr
robustMahalanobis <- function(delta) {
    
    # Compute robust covariance
    rCov <- robustbase::covMcd(delta)
    
    sqMan <- matrix(NA, nrow = nrow(delta), ncol = 1)
    for (i in seq.int(nrow(sqMan))) {
        sqMan[i] <- (delta[i,] - colMeans(delta)) %*% solve(rCov$cov) %*% (delta[i,] - colMeans(delta))
    }
    
    return(sqMan)
} 
##' @title Compute the reproducibility score 
##' @param x Numeric vector to compute reproducibility score
##' @param y Numeric vector to compute reprodducibility score
##' @param method Correlation method. Default is Pearson
##' @return The R score
##' @md
##' @rdname method-mr
reprodScore <- function(x, y, method = c("pearson")) {
    
    Rscore <- matrix(NA, nrow = nrow(x), ncol = 1)
    for (i in seq.int(nrow(x))){
        Rscore[i] <- cor(x[i,], y[i,], method = method)
    }
    
    return(Rscore)
}
##' @title Apply the MR method to spatial proteomics data
##' @param objectCond1 A list of [`MSnbase::MSnSet`]s where each is an experimental
##' replicate for the first condition, usually a control
##' @param objectCond2 A list of [`MSnbase::MSnSet`]s where each is an experimental
##' replicate for the second condition, usually a treatment
##' @param plot Whether to generate an MR plot as a side effect. 
##' @return The MR score of the Ithzak et al. 2016/2017
##' @md
##' @rdname method-mr
mrMethod <- function(objectCond1,
                     objectCond2,
                     plot = TRUE) {
    
    stopifnot(class(objectCond1[[1]]) == "MSnSet")
    stopifnot(class(objectCond2[[1]]) == "MSnSet")
    stopifnot(length(objectCond1) >= 2)
    stopifnot(length(objectCond2) >= 2)
    
    exprsObj1 <- lapply(objectCond1, exprs)
    exprsObj2 <- lapply(objectCond2, exprs)
    
    
    # Compute delta matrices
    delta <- lapply(1:length(exprsObj1),
                    function(n) exprsObj1[[n]] - exprsObj2[[n]])
    
    # Compute robust mahalanobis distance
    sqMan <- lapply(delta, function(x) robustMahalanobis(x))
    
    # Under null distance is chi squared distributed (disputed)
    df <- ncol(delta[[1]])
    pchi <- lapply(sqMan, function(x) pchisq(x, df = df, lower.tail = FALSE))
    
    # As stated take the maximum p-value from each test
    summarisedPvalue <- apply(do.call(cbind, pchi), 1, max)
    
    # Adjust for multiple testing of the cubed pValue (not recommended but
    # is part of the mr method)
    summarisedPvalue <- p.adjust(summarisedPvalue^3, method = "BH")
    
    # add Names
    names(summarisedPvalue) <- rownames(objectCond1[[1]])
    
    # Compute M-score
    Mscore <- -log10(summarisedPvalue)
    
    # compute reproducibility sore
    Rscore1 <- lapply(1:length(exprsObj1),
                     function(x) reprodScore(x = delta[[x]], y = delta[[(x)%%length(objectCond1) + 1]]))
    
    # take smallest score
    Rscore <- apply(do.call(cbind, Rscore1), 1, function(x) min(x, na.rm = TRUE))
    
    # add Names
    names(Rscore) <- names(Mscore) <- rownames(objectCond1)
    
    if (plot == TRUE) {
        plot(Mscore, Rscore, pch = 19)
    }
    
    ## Compute ROC, with normalized scores
    MRscore <- ((Rscore - mean(Rscore))/sd(Rscore)) * ((Mscore - mean(Mscore[!is.infinite(Mscore)]))/sd(Mscore[!is.infinite(Mscore)]))
    names(MRscore) <- rownames(objectCond1[[1]])
    MRscore[is.infinite(MRscore)] <- max(MRscore[!is.infinite(MRscore)])
    
    return(list(Rscore = Rscore, Mscore = Mscore, MRscore = MRscore, sqMan = sqMan))
}
