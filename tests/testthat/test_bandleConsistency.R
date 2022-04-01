context("bandle function test")

library("pRolocdata")
data("tan2009r1")
set.seed(1)
tansim <- sim_dynamic(object = tan2009r1, 
                      numRep = 4L,
                      numDyn = 100L)
d2 <- d1 <- tansim$lopitrep
control2 <- control1 <- d1[1:2]
treatment2 <- treatment1 <- d1[3:4]
i <- which(fvarLabels(d1[[1]]) == "markers")
stopifnot(length(i) == 1)
fvarLabels(control2[[1]])[i] <- "xx"
fvarLabels(control2[[2]])[i] <- "xx"
# fvarLabels(control2[[3]])[i] <- "xx"
fvarLabels(treatment2[[1]])[i] <- "xx"
fvarLabels(treatment2[[2]])[i] <- "xx"
# fvarLabels(treatment2[[3]])[i] <- "xx"

.times <- 2
.seed <- 1

test_that("bandle consistency", {
    .numIter <- 5
    set.seed(1)
    gpParams <- lapply(tansim$lopitrep, 
                       function(x) fitGPmaternPC(x, hyppar = matrix(c(0.5, 1, 100), nrow = 1)))
    mcmc1 <- bandle(objectCond1 = control1, objectCond2 = treatment1, gpParams = gpParams,
                    fcol = "markers", numIter = .numIter, burnin = 1L, thin = 2L,
                    numChains = 2, BPPARAM = SerialParam(RNGseed = 1))
    set.seed(1)
    mcmc2 <- bandle(objectCond1 = control2, objectCond2 = treatment2, gpParams = gpParams,
                    fcol = "xx", numIter = .numIter, burnin = 1L, thin = 2L,
                    numChains = 2, BPPARAM = SerialParam(RNGseed = 1))
    mcmc1 <- bandleProcess(mcmc1)
    mcmc2 <- bandleProcess(mcmc2)
    ans1 <- bandlePredict(objectCond1 = control1, objectCond2 = treatment1, params = mcmc1, fcol = "markers")
    ans2 <- bandlePredict(objectCond1 = control2, objectCond2 = treatment2, params = mcmc2, fcol = "xx")
    expect_equal(ans2, ans1, check.attributes = FALSE)
})

test_that("bandle consistency 2", {
    .numIter <- 5
    set.seed(1)
    gpParams <- lapply(tansim$lopitrep, 
                       function(x) fitGPmatern(x))
    mcmc1 <- bandle(objectCond1 = control1, objectCond2 = treatment1, gpParams = gpParams,
                    fcol = "markers", numIter = .numIter, burnin = 1L, thin = 2L,
                    numChains = 2, BPPARAM = SerialParam(RNGseed = 1))
    set.seed(1)
    mcmc2 <- bandle(objectCond1 = control2, objectCond2 = treatment2, gpParams = gpParams,
                    fcol = "xx", numIter = .numIter, burnin = 1L, thin = 2L,
                    numChains = 2, BPPARAM = SerialParam(RNGseed = 1))
    mcmc1 <- bandleProcess(mcmc1)
    mcmc2 <- bandleProcess(mcmc2)
    ans1 <- bandlePredict(objectCond1 = control1, objectCond2 = treatment1, params = mcmc1, fcol = "markers")
    ans2 <- bandlePredict(objectCond1 = control2, objectCond2 = treatment2, params = mcmc2, fcol = "xx")
    expect_equal(ans2, ans1, check.attributes = FALSE)
})

test_that("bandle consistency 3", {
    .numIter <- 5
    set.seed(1)
    gpParams <- lapply(tansim$lopitrep, 
                       function(x) fitGP(x))
    mcmc1 <- bandle(objectCond1 = control1, objectCond2 = treatment1, gpParams = gpParams,
                    fcol = "markers", numIter = .numIter, burnin = 1L, thin = 2L,
                    numChains = 2, BPPARAM = SerialParam(RNGseed = 1))
    set.seed(1)
    mcmc2 <- bandle(objectCond1 = control2, objectCond2 = treatment2, gpParams = gpParams,
                    fcol = "xx", numIter = .numIter, burnin = 1L, thin = 2L,
                    numChains = 2, BPPARAM = SerialParam(RNGseed = 1))
    mcmc1 <- bandleProcess(mcmc1)
    mcmc2 <- bandleProcess(mcmc2)
    ans1 <- bandlePredict(objectCond1 = control1, objectCond2 = treatment1, params = mcmc1, fcol = "markers")
    ans2 <- bandlePredict(objectCond1 = control2, objectCond2 = treatment2, params = mcmc2, fcol = "xx")
    expect_equal(ans2, ans1, check.attributes = FALSE)
})

test_that("bandle consistency 4", {
    .numIter <- 5
    set.seed(1)
    gpParams <- lapply(tansim$lopitrep, 
                       function(x) fitGPmaternPC(x, hyppar = matrix(c(0.5, 1, 100), nrow = 1)))
    mcmc1 <- bandle(objectCond1 = control1, objectCond2 = treatment1, gpParams = gpParams,
                    fcol = "markers", numIter = .numIter, burnin = 1L, thin = 2L,
                    numChains = 2, BPPARAM = SerialParam(RNGseed = 1), pg = TRUE)
    set.seed(1)
    mcmc2 <- bandle(objectCond1 = control2, objectCond2 = treatment2, gpParams = gpParams,
                    fcol = "xx", numIter = .numIter, burnin = 1L, thin = 2L,
                    numChains = 2, BPPARAM = SerialParam(RNGseed = 1), pg = TRUE)
    mcmc1 <- bandleProcess(mcmc1)
    mcmc2 <- bandleProcess(mcmc2)
    ans1 <- bandlePredict(objectCond1 = control1, objectCond2 = treatment1, params = mcmc1, fcol = "markers")
    ans2 <- bandlePredict(objectCond1 = control2, objectCond2 = treatment2, params = mcmc2, fcol = "xx")
    expect_equal(ans2, ans1, check.attributes = FALSE)
})

test_that("bandle consistency 5", {
    .numIter <- 20
    set.seed(1)
    gpParams <- lapply(tansim$lopitrep, 
                       function(x) fitGPmaternPC(x, hyppar = matrix(c(0.5, 1, 100), nrow = 1)))
    mcmc1 <- bandle(objectCond1 = control1, objectCond2 = treatment1, gpParams = gpParams,
                    fcol = "markers", numIter = .numIter, burnin = 1L, thin = 2L,
                    numChains = 2, BPPARAM = SerialParam(RNGseed = 1), hyperLearn = "MH", hyperIter = 5,
                    pcPrior = matrix(c(0.5, 1, 100), nrow = 11, ncol = 3, byrow = TRUE))
    set.seed(1)
    mcmc2 <- bandle(objectCond1 = control2, objectCond2 = treatment2, gpParams = gpParams,
                    fcol = "xx", numIter = .numIter, burnin = 1L, thin = 2L,
                    numChains = 2, BPPARAM = SerialParam(RNGseed = 1), hyperLearn = "MH", hyperIter = 1,
                    pcPrior = matrix(c(0.5, 1, 100), nrow = 11, ncol = 3, byrow = TRUE))
    mcmc1 <- bandleProcess(mcmc1)
    mcmc2 <- bandleProcess(mcmc2)
    ans1 <- bandlePredict(objectCond1 = control1, objectCond2 = treatment1, params = mcmc1, fcol = "markers")
    ans2 <- bandlePredict(objectCond1 = control2, objectCond2 = treatment2, params = mcmc2, fcol = "xx")
    expect_equal(ans2, ans1, check.attributes = FALSE)
})
