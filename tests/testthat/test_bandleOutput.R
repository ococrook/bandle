context("bandle output tests")

library("pRolocdata")
data("tan2009r1")
set.seed(1)
tansim <- sim_dynamic(object = tan2009r1, 
                      numRep = 4L,
                      numDyn = 100L)

control <- tansim[[1]][1:2]
treatment <- tansim[[1]][3:4]

pData(control[[1]])$Condition <- pData(control[[2]])$Condition <- "Control"
pData(treatment[[1]])$Condition <- pData(treatment[[2]])$Condition <- "Treatment"
pData(control[[1]])$Replicate <- pData(treatment[[1]])$Replicate <- 1
pData(control[[2]])$Replicate <- pData(treatment[[2]])$Replicate <- 2


test_that("bandlePredict pData output", {
  .numIter <- 5
  set.seed(1)
  
  ## before running bandle
  pdata_con1 <- pData(control[[1]])
  pdata_con2 <- pData(control[[2]])
  pdata_tr1 <- pData(treatment[[1]])
  pdata_tr2 <- pData(treatment[[2]])
  
  ## run
  gpParams <- lapply(tansim$lopitrep, 
                     function(x) fitGPmaternPC(x, hyppar = matrix(c(0.5, 1, 100), nrow = 1)))
  res <- bandle(objectCond1 = control, objectCond2 = treatment, gpParams = gpParams,
                  fcol = "markers", numIter = .numIter, burnin = 1L, thin = 2L,
                  numChains = 2, BPPARAM = SerialParam(RNGseed = 1))
  
  res <- bandleProcess(res)
  
  ## get results and compare
  ans <- bandlePredict(objectCond1 = control, 
                        objectCond2 = treatment, 
                        params = res, fcol = "markers")
  res_pdata_con1 <- pData(ans[[1]][[1]])
  res_pdata_con2 <- pData(ans[[1]][[2]])
  res_pdata_tr1 <- pData(ans[[2]][[1]])
  res_pdata_tr2 <- pData(ans[[2]][[2]])
  
  expect_equal(pdata_con1, res_pdata_con1, check.attributes = FALSE)
  expect_equal(pdata_con2, res_pdata_con2, check.attributes = FALSE)
  expect_equal(pdata_tr1, res_pdata_tr1, check.attributes = FALSE)
  expect_equal(pdata_tr2, res_pdata_tr2, check.attributes = FALSE)
})

test_that("bandlePredict fData output", {
  .numIter <- 5
  set.seed(1)
  
  ## before running bandle get fData
  fdata_con1 <- fvarLabels(control[[1]])
  fdata_con2 <- fvarLabels(control[[2]])
  fdata_tr1 <- fvarLabels(treatment[[1]])
  fdata_tr2 <- fvarLabels(treatment[[2]])
  
  ## run
  gpParams <- lapply(tansim$lopitrep, 
                     function(x) fitGPmaternPC(x, hyppar = matrix(c(0.5, 1, 100), nrow = 1)))
  res <- bandle(objectCond1 = control, objectCond2 = treatment, gpParams = gpParams,
                fcol = "markers", numIter = .numIter, burnin = 1L, thin = 2L,
                numChains = 2, BPPARAM = SerialParam(RNGseed = 1))
  
  res <- bandleProcess(res)
  
  ## get results and compare
  ans <- bandlePredict(objectCond1 = control, 
                       objectCond2 = treatment, 
                       params = res, fcol = "markers")
  res_fdata_con1 <- fvarLabels(ans[[1]][[1]])
  res_fdata_con2 <- fvarLabels(ans[[1]][[2]])
  res_fdata_tr1 <- fvarLabels(ans[[2]][[1]])
  res_fdata_tr2 <- fvarLabels(ans[[2]][[2]])
  
  bandle_cols <- c("bandle.allocation", "bandle.probability", 
                   "bandle.probability.lowerquantile",
                   "bandle.probability.upperquantile", "bandle.mean.shannon", 
                   "bandle.differential.localisation", "bandle.outlier",
                   "bandle.joint")
  
  ## check bandle results are only appended to the FIRST MSnSet of each condition 
  expect_equal(res_fdata_con1, c(fdata_con1, bandle_cols), check.attributes = FALSE)
  expect_equal(res_fdata_tr1, c(fdata_tr1, bandle_cols), check.attributes = FALSE)
  
  ## check the other fvarLabels are the same 
  expect_equal(res_fdata_con2, fdata_con2, check.attributes = FALSE)
  expect_equal(res_fdata_tr2, fdata_tr2, check.attributes = FALSE)
})
