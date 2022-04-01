##' These function implement the MR method of Itzhak et al
##' 
##' @title robust Mahalanobis distance
##' @param delta The difference profile to compute the squared mahalanobis distance
##' @return The squared Mahalanobis distance 
##' @md
##' @examples 
##' ## Generate some example data
##' library("pRolocdata")
##' data("tan2009r1")
##' set.seed(1)
##' tansim <- sim_dynamic(object = tan2009r1, 
##'                       numRep = 4L,
##'                       numDyn = 100L)
##' data <- tansim$lopitrep
##' control <- data[1:2]
##' treatment <- data[3:4]
##' 
##' ## compute delta matrix
##' deltaMatrix <- exprs(control[[1]]) - exprs(treatment[[1]])
##' res <- bandle:::robustMahalanobis(deltaMatrix)
##' @rdname method-mr
robustMahalanobis <- function(delta) {
    
    # Compute robust covariance
    rCov <- robustbase::covMcd(delta)
    
    sqMan <- matrix(NA, nrow = nrow(delta), ncol = 1)
    for (i in seq.int(nrow(sqMan))) {
        sqMan[i] <- (delta[i,] - colMeans(delta)) %*% chol2inv(chol(rCov$cov + 10^{-8})) %*% (delta[i,] - colMeans(delta))
    }
    
    return(sqMan)
} 

##' @title Compute the reproducibility score 
##' @param x Numeric vector to compute reproducibility score
##' @param y Numeric vector to compute reproducibility score
##' @param method Correlation method. Default is Pearson
##' @return The R score
##' @md
##' @examples
##' ##' @examples 
##' ## Generate some example data
##' library("pRolocdata")
##' data("tan2009r1")
##' set.seed(1)
##' tansim <- sim_dynamic(object = tan2009r1, 
##'                       numRep = 4L,
##'                       numDyn = 100L)
##' data <- tansim$lopitrep
##' control <- data[1:2]
##' treatment <- data[3:4]
##' 
##' ## compute delta matrix
##' deltaMatrix1 <- exprs(control[[1]]) - exprs(treatment[[1]])
##' deltaMatrix2 <- exprs(control[[2]]) - exprs(treatment[[2]])
##' mr_score <- bandle:::reprodScore(deltaMatrix1, deltaMatrix2)
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
##' @return The MR score of the Ithzak et al. 2016/2017
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
##' mr1 <- mrMethod(objectCond1 = control1, objectCond2 = treatment1)
##' plot(mr1$Mscore, mr1$Rscore, pch = 21, 
##'      xlab = "MScore", ylab = "RScore")
##' @rdname method-mr
mrMethod <- function(objectCond1,
                     objectCond2,
                     method = "2017") {
    
    stopifnot(is(objectCond1[[1]], "MSnSet"))
    stopifnot(is(objectCond2[[1]], "MSnSet"))
    stopifnot(length(objectCond1) >= 2)
    stopifnot(length(objectCond2) >= 2)
    
    exprsObj1 <- lapply(objectCond1, exprs)
    exprsObj2 <- lapply(objectCond2, exprs)
    
    
    # Compute delta matrices
    delta <- lapply(seq.int(length(exprsObj1)),
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
    if (method == "2017") {
    summarisedPvalue <- p.adjust(summarisedPvalue^3, method = "BH")
    } else if (method == "2016") {
    summarisedPvalue <- p.adjust(summarisedPvalue, method = "BH")    
    }
    # add Names
    names(summarisedPvalue) <- rownames(objectCond1[[1]])
    
    # Compute M-score
    Mscore <- -log10(summarisedPvalue)
    
    # compute reproducibility sore
    Rscore1 <- lapply(seq.int(length(exprsObj1)),
                     function(x) reprodScore(x = delta[[x]], y = delta[[(x)%%length(objectCond1) + 1]]))
    
    # take smallest score
    Rscore <- apply(do.call(cbind, Rscore1), 1, function(x) min(x, na.rm = TRUE))
    
    # add Names
    names(Rscore) <- names(Mscore) <- rownames(objectCond1)
    
    # if (plot == TRUE) {
    #     plot(Mscore, Rscore, pch = 19)
    # }
    
    ## Compute ROC, with normalized scores
    MRscore <- ((Rscore - mean(Rscore[!is.infinite(Rscore)]))/sd(Rscore[!is.infinite(Rscore)])) * ((Mscore - mean(Mscore[!is.infinite(Mscore)]))/sd(Mscore[!is.infinite(Mscore)]))
    names(MRscore) <- rownames(objectCond1[[1]])
    MRscore[is.infinite(MRscore)] <- max(MRscore[!is.infinite(MRscore)])
    
    return(list(Rscore = Rscore, Mscore = Mscore, MRscore = MRscore, sqMan = sqMan))
}