##' These function implement the MR method of Itzhak et al
##' 
##' @title robust Mahalanobis distance
##' @params delta The difference profile to compute the squared mahalanobis distance
##' @return The squared Mahalanobis distance 
##' @md
##' @rdname method-mr
robustMahalanobis <- function(delta) {
    
    # Compute robust covariance
    rCov <- covMcd(delta)
    
    sqMan <- matrix(NA, nrow = nrow(delta), ncol = 1)
    for (i in seq.int(nrow(sqMan))) {
        sqMan[i] <- (delta[i,] - colMeans(delta)) %*% solve(rCov$cov) %*% (delta[i,] - colMeans(delta))
    }
    
    return(sqMan)
} 
##' @title Compute the reproducibility score 
##' @params x Numeric vector to compute reproducibility score
##' @params y Numeric vector to compute reprodducibility score
##' @params method Correlation method. Default is Pearson
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
##' @params plot Whether to generate an MR plot as a side effect. 
##' @return The MR score of the Ithzak et al. 2016/2017
##' @md
##' @rdname method-mr
mrMethod <- function(objectCond1,
                     objectCond2,
                     plot = TRUE) {
    
    stopifnot(class(objectCond1)[[1]] == "MSnSet")
    stopifnot(class(objectCond2)[[1]] == "MSnSet")
    stopifnot(length(objectCond1) == 3)
    stopifnot(length(objectCond2) == 3)
    
    test1 <- objectCond1[[1]]
    test1b <- objectCond1[[2]]
    test1c <- objectCond1[[3]]
    
    test2 <- objectCond2[[1]]
    test2b <- objectCond2[[2]]
    test2c <- objectCond2[[3]]
    
    # Compute delta matrices
    delta1 <- exprs(test1) - exprs(test2)
    delta2 <- exprs(test1b) - exprs(test2b)
    delta3 <- exprs(test1c) - exprs(test2c)
    
    # Compute robust mahalanobis distance
    sqMan1 <- robustMahalanobis(delta = delta1)
    sqMan2 <- robustMahalanobis(delta = delta2)
    sqMan3 <- robustMahalanobis(delta = delta3)
    
    
    # Under null distance is chi squared distributed
    df <- ncol(delta1)
    pchi1 <- pchisq(sqMan1, df = df, lower.tail = FALSE)
    pchi2 <- pchisq(sqMan2, df = df, lower.tail = FALSE)
    pchi3 <- pchisq(sqMan3, df = df, lower.tail = FALSE)
    
    
    # As stated take the maximum p-value from each test
    summarisedPvalue <- apply(cbind(pchi1, pchi2, pchi3), 1, max)
    
    # Adjust for multiple testing of the cubed pValue (not recommended but
    # is part of the mr method)
    summarisedPvalue <- p.adjust(summarisedPvalue^3, method = "BH")
    
    # add Names
    names(summarisedPvalue) <- rownames(test1)
    
    # Compute M-score
    Mscore <- -log10(summarisedPvalue)
    
    # compute reproducibility sore
    Rscore1 <- reprodScore(x = delta1, y = delta2)
    Rscore2 <- reprodScore(x = delta2, y = delta3)
    Rscore3 <- reprodScore(x = delta3, y = delta1)
    
    # take smallest score
    Rscore <- apply(cbind(Rscore1, Rscore2, Rscore3), 1, function(x) min(x, na.rm = TRUE))
    
    # add Names
    names(Rscore) <- names(Mscore) <- rownames(test1)
    
    if (plot == TRUE) {
        plot(Mscore, Rscore, pch = 19)
    }
    
    ## Compute ROC, with normalized scores
    MRscore <- ((Rscore - mean(Rscore))/sd(Rscore)) * ((Mscore - mean(Mscore[!is.infinite(Mscore)]))/sd(Mscore[!is.infinite(Mscore)]))
    names(MRscore) <- rownames(test1)
    MRscore[is.infinite(MRscore)] <- max(MRscore[!is.infinite(MRscore)])
    
    return(list(Rscore = Rscore, Mscore = Mscore, MRscore = MRscore, sqMan1 = sqMan1, sqMan2 = sqMan2, sqMan3 = sqMan3))
}