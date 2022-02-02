library(pRolocdata)
data("tan2009r1")
set.seed(1)
.numRep <- 6L
.numDyn <- 100L

tansim <- sim_dynamic(object = tan2009r1,
                      numRep = .numRep, numDyn = .numDyn, )


test_that("sim dynamic", {
   
    expect_length(tansim, 2)
    expect_length(tansim$lopitrep, .numRep)
    expect_equal(featureNames(tansim$lopitrep[[1]]),
                 featureNames(tansim$lopitrep[[2]]))
    expect_equal(featureNames(tansim$lopitrep[[1]]),
                 featureNames(tan2009r1))
    expect_length(tansim$perm1_names, .numDyn)
    expect_true(all(tansim$perm1_names %in% featureNames(tan2009r1)))

})
