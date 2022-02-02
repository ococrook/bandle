context("bandle function test")

library("pRolocdata")
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

test_that("prior consistency", {
    
    set.seed(1)
    ans1 <- prior_pred_dir(object = control1[[1]])
    set.seed(1)
    ans2 <- prior_pred_dir(object = control2[[1]], fcol = "xx")
    expect_equal(ans2, ans1, check.attributes = FALSE)
    
    set.seed(1)
    ans3 <- prior_pred_pg(objectCond1 = control1[[1]],
                          objectCond2 = treatment1[[1]], fcol = "markers")
    set.seed(1)
    ans4 <- prior_pred_pg(objectCond1 = control2[[1]],
                          objectCond2 = treatment2[[1]], fcol = "xx")
    expect_equal(ans3, ans4, check.attributes = FALSE)    
})
