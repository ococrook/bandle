
library(pRolocdata)
data("tan2009r1")
set.seed(1)
tansim <- sim_dynamic(object = tan2009r1, 
                     numRep = 6L,
                    numDyn = 100L)
gpParams <- lapply(tansim$lopitrep, function(x) 
fitGPmaternPC(x, hyppar = matrix(c(0.5, 1, 100), nrow = 1)))
d1 <- tansim$lopitrep
control1 <- d1[1:3]
treatment1 <- d1[4:6]

test_that("differential localisation computation", {

    set.seed(1)
    mcmc1 <- bandle(objectCond1 = control1,
                    objectCond2 = treatment1,
                    gpParams = gpParams,
                    fcol = "markers",
                    numIter = 10L, burnin = 1L, thin = 2L, numChains = 2,
                    BPPARAM = SerialParam(RNGseed = 1))
    mcmc1 <- bandleProcess(mcmc1)
    dp1 <- diffLocalisationProb(mcmc1)
    expect_length(dp1, length(rownames(unknownMSnSet(object = tan2009r1, fcol = "markers"))))
    expect_true(all(dp1 <= 1))
    expect_true(all(dp1 >= 0))

})

test_that("differential localisation computation 2", {
    
    set.seed(1)
    mcmc1 <- bandle(objectCond1 = control1,
                    objectCond2 = treatment1,
                    gpParams = gpParams,
                    fcol = "markers",
                    numIter = 10L, burnin = 1L, thin = 2L, numChains = 2,
                    BPPARAM = SerialParam(RNGseed = 1))
    mcmc1 <- bandleProcess(mcmc1)
    .top <- 20
    .bootsample <- 100
    bdp <- bootstrapdiffLocprob(mcmc1, top = .top, Bootsample = .bootsample,
                                decreasing = TRUE)
    expect_length(bdp[1,], .bootsample)
    expect_length(bdp[,1], .top)
    expect_equal(order(bdp[,1], decreasing = TRUE), seq.int(.top))
    
})

test_that("differential localisation computation 3", {
    
    set.seed(1)
    mcmc1 <- bandle(objectCond1 = control1,
                    objectCond2 = treatment1,
                    gpParams = gpParams,
                    fcol = "markers",
                    numIter = 10L, burnin = 1L, thin = 2L, numChains = 2,
                    BPPARAM = SerialParam(RNGseed = 1))
    mcmc1 <- bandleProcess(mcmc1)
    .top <- 20
    .nsample <- 100
    dp <- binomialDiffLocProb(mcmc1, top = .top, nsample = .nsample,
                              decreasing = TRUE)
    probs <- diffLocalisationProb(params = mcmc1)
    probs <- probs[order(probs, decreasing = TRUE)]
    expect_length(dp[1,], .nsample)
    expect_length(dp[,1], .top)
    expect_equal(names(probs)[seq.int(.top)], rownames(dp))
    
})


