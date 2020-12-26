## Simulating dynamic experiment using bootstrap method

sim_dynamic <- function(object,
                        subnum = NULL,
                        knn_par = 10,
                        fcol = "markers",
                        numRep = 6,
                        method = "wild",
                        batch = FALSE,
                        frac_perm = FALSE,
                        nu = 2,
                        numDyn = 20) {

lopitmarkerset <- markerMSnSet(object)
lopitunknownset <- unknownMSnSet(object)

Nunk <- nrow(lopitunknownset)

if (!is.null(subnum)){
  subs <- sample.int(Nunk, size = subnum)
  object <- combine(lopitmarkerset, lopitunknownset[subs, ])
}

## Perform K-NN with K classification to get clustering
object <- knnClassification(object = object, k = knn_par, fcol = fcol)

orgM <- meanOrganelle(mydata = object, fcol = fcol)$M

## Use means as fitted function for regression and compute residuals
residual <- matrix(NA, nrow = nrow(object), ncol = ncol(object))
residual <- exprs(object) - orgM[fData(object)[, "knn"], ]

## Bootstrap residuals
numRep <- numRep # number of additional datasets

lopitrep <- list()
K <- length(getMarkerClasses(object, fcol = fcol))

# Methods for creating replicates
if (method == "const") {
  nu <- nu
} else if (method == "rand") {
  nu <- runif(nrow(object), 1, nu)
} else if (method == "wild") {
  nu <- runif(K, 1, nu)[fData(object)$knn] # cluster specific residual inflated noise
}

for (i in seq.int(numRep)) {
  
  myobject <- object
  idxBoot <- matrix(sample.int(ncol(myobject),
                               size = ncol(myobject) * nrow(myobject),
                               replace = TRUE),
                               nrow = nrow(myobject), ncol(myobject))
  bootres <- matrix(NA,  nrow = nrow(myobject), ncol(myobject))
  for (k in seq.int(nrow(myobject))) {
    bootres[k, ] <- residual[k, idxBoot[k, ]]
  }
  lopitrep[[i]] <- myobject
  newExprs <- orgM[fData(myobject)$knn, ] + nu * bootres ## add boostraps
  colnames(newExprs) <- colnames(myobject)
  rownames(newExprs) <- rownames(myobject)
  exprs(lopitrep[[i]]) <- newExprs
  
}

if (batch == "rand") {
  for (i in seq.int(numRep)){
    batch_effect <- rgamma(1, 1, 1) # sample batch effect
    g <- sample(c(1:ncol(myobject)), size = 1) # sample which fraction to add the effect to
    exprs(lopitrep[[i]])[, g] <- exprs(lopitrep[[i]])[, g] + batch_effect
  }
} else if (batch == "systematic") {
  g <- sample(c(1:ncol(myobject)), size = 1) # sample which fraction to add the effect to
  batch_effect <- rgamma(1, 1, 1) # sample batch effect
  for (i in seq.int(numRep)) {
    exprs(lopitrep[[i]])[, g] <- exprs(lopitrep[[i]])[, g] + batch_effect
  }
}

if (frac_perm == TRUE) {
  for (i in seq.int(numRep)){
    permute <- sample.int(ncol(myobject))
    exprs(lopitrep[[i]]) <- exprs(lopitrep[[i]])[, permute]
  }
}


## simulate translocation events
mlopit <- unknownMSnSet(object)
perm1 <- sample.int(nrow(mlopit), size = numDyn, replace = TRUE)
perm1_names <- rownames(mlopit)[perm1]

replicateStats <- lapply(lopitrep, function(x) meanOrganelle(x, fcol = fcol))

# Simulate dynamics # gahh need to use correct names

for (i in perm1_names) {
  K <- length(getMarkerClasses(object, fcol = fcol))
  r <- sample.int(K, size = 1)
  for (j in 4:6) {
    exprs(lopitrep[[j]])[i, ] <- replicateStats[[j]]$M[r, ] + 
      rnorm(ncol(object), mean = 0, sd = sqrt(replicateStats[[j]]$V[r, ]))
  }
}

return(list(lopitrep = lopitrep, perm1_names = perm1_names))
}
