simres <- sim_dynamic(object = tan2009r1)

tan2009rep <- simres$lopitrep
perm1_names <- simres$perm1_names

# computational methods
mrres <- mrMethod(c(tan2009rep))
par(mfrow = c(3,4))
pc_prior <- matrix(NA, ncol = 3, 11)
pc_prior[seq.int(1:11),] <- matrix(rep(c(0.05,60,100), each = 11), ncol = 3)

resdynpg <- bandle(objectCond1 = c(tan2009rep[1:3]),
                   objectCond2 = tan2009rep[4:6],
                   numIter = 100,
                   burnin = 10,
                    pg = FALSE,
                   hyperLearn = "LBFGS", maternCov = TRUE, nu = 2, pcPrior = pc_prior, numChains = 1)
