library(pRolocdata)
data("tan2009r1")
data("tan2009r2")
set.seed(1)
tansim <- sim_dynamic(object = tan2009r1, 
                      numRep = 3L,
                      numDyn = 100L)
tansim2 <- sim_dynamic(object = tan2009r1, 
                      numRep = 3L,
                      numDyn = 100L)
d1 <- c(tansim$lopitrep, tansim2$lopitrep)

test_that("mr Method", {
    
    control1 <- d1[1:3]
    treatment1 <- d1[4:6]
    mr1 <- mrMethod(objectCond1 = control1, objectCond2 = treatment1)
    
    expect_length(mr1, length(control1) + 1)
    expect_length(mr1[[length(control1) + 1]], length(control1))
    
    expect_length(mr1[[1]], length(rownames(tan2009r1)))
    expect_true(all(mr1[[1]] >= -1))
    expect_true(all(mr1[[1]] <= 1))
    expect_true(all(mr1[[2]] >= 0 ))
    
})
