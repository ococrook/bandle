context("bandle function test")

data("tan2009r1")
set.seed(1)
tansim <- sim_dynamic(object = tan2009r1, 
                      numRep = 6L,
                      numDyn = 100L)
d2 <- d1 <- tansim$lopitrep
control2 <- control1 <- d1[1:3]
treatment2 <- treatment1 <- d1[4:6]
i <- which(fvarLabels(d1[[1]]) == "markers")
stopifnot(length(i) == 1)
fvarLabels(control2[[1]])[i] <- "xx"
fvarLabels(control2[[2]])[i] <- "xx"
fvarLabels(control2[[3]])[i] <- "xx"
fvarLabels(treatment2[[1]])[i] <- "xx"
fvarLabels(treatment2[[2]])[i] <- "xx"
fvarLabels(treatment2[[3]])[i] <- "xx"

.times <- 2
.seed <- 1

test_that("bandle consistency", {
    .numIter <- 5
    set.seed(1)
    gpParams <- lapply(tansim$lopitrep, function(x) fitGPmaternPC(x, hyppar = c(0.5, 1, 100)))
    mcmc1 <- bandle(objectCond1 = control1, objectCond2 = treatment1, gpParams = gpParams,
                    fcol = "markers", numIter = .numIter, burnin = 1L, thin = 2L,
                    numChains = 2, BPPARAM = SerialParam(RNGseed = 1))
    set.seed(1)
    mcmc2 <- bandle(objectCond1 = control2, objectCond2 = treatment2, gpParams = gpParams,
                    fcol = "xx", numIter = .numIter, burnin = 1L, thin = 2L,
                    numChains = 2, BPPARAM = SerialParam(RNGseed = 1))
    ans1 <- bandlePredict(objectCond1 = control1, objectCond2 = treatment1, params = mcmc1, fcol = "markers")
    ans2 <- bandlePredict(objectCond1 = control2, objectCond2 = treatment2, params = mcmc2, fcol = "xx")
    expect_equal(ans1[[1]], ans2[[1]], check.attributes = FALSE)
})
