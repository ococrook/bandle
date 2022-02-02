##' A function to simulate dynamic spatial proteomics data using a bootstrap method
##' 
##' 
##' @title Generate a dynamic spatial proteomics experiment
##' @param object A instance of class `MSnSet` from which to generate a spatial
##' proteomics dataset.
##' @param subsample how many proteins to subsample to speed up analysis. Default is NULL.
##' @param knn_par the number of nearest neighbours to use in KNN classification to simulate
##'  dataset.
##'  Default is 10
##' @param fcol feature column to indicate markers. Default is "markers". Proteins
##' with unknown localisations must be encoded as "unknown". 
##' @param numRep The total number of datasets to generate. Default is 6. An
##' integer must be provided
##' @param method The bootstrap method to use to simulate dataset. Default is "wild".
##'  refer to BANDLE paper for more details.
##' @param batch Whether or not to include batch effects. Default is FALSE.
##' @param frac_perm whether or not to permute the fractions. Default is FALSE
##' @param nu parameter to generate residual inflated noise. Default is 2. See BANDLE
##' paper for more details
##' @param numDyn An integer number of protein to simulate dynamic transitions. Default is 20 
##' @return returns simulate dynamic lopit datasets and the name of the relocalated protein.  
##' @md
##' 
##' @examples
##' library(pRolocdata)
##' data("tan2009r1")
##' set.seed(1)
##' tansim <- sim_dynamic(object = tan2009r1, numRep = 6L, numDyn = 100L)
##' 
##' 
##' @rdname bandle-sim
sim_dynamic <- function(object,
                        subsample = NULL,
                        knn_par = 10L,
                        fcol = "markers",
                        numRep = 6L,
                        method = "wild",
                        batch = FALSE,
                        frac_perm = FALSE,
                        nu = 2,
                        numDyn = 20L) {
  # checks
  stopifnot("object is not an instance of class MSnSet"=is(object, "MSnSet"))
  stopifnot("number of Replicates is not valid, provide an integer"=is(numRep, "integer"))
  stopifnot("number of Translocations is not valid, provide an integer"=is(numDyn, "integer"))
  stopifnot("Knn parameers must be an integer"=is(knn_par, "integer"))
  stopifnot("method is not valid"=method %in% c("wild", "rand", "const"))
  stopifnot("batch is not valid"=batch %in% c(FALSE, "rand", "systematic"))
    
  # generate marker and unknown MSnSets
  lopitmarkerset <- markerMSnSet(object, fcol = fcol)
  lopitunknownset <- unknownMSnSet(object, fcol = fcol)

  # number of proteins that are unannotated
  Nunk <- nrow(lopitunknownset)

  # if subsampling start to subsample
  if (!is.null(subsample)){
    subs <- sample.int(Nunk, size = subsample)
    object <- combine(lopitmarkerset, lopitunknownset[subs, ])
  }

  ## Perform K-NN with K classification to get clustering
  object <- knnClassification(object = object, k = knn_par, fcol = fcol)

  # Get Mean of each organelle
  orgM <- meanOrganelle(object = object, fcol = fcol)$M

  ## Use means as fitted function for regression and compute residuals
  residual <- matrix(NA, nrow = nrow(object), ncol = ncol(object))
  residual <- exprs(object) - orgM[fData(object)[, "knn"], ]

  ## Bootstrap residuals
  numRep <- numRep # number of additional datasets

  ## storage for datasets
  lopitrep <- list()
  K <- length(getMarkerClasses(object, fcol = fcol))

  
  # Methods for creating replicates, constant, random or wild boostrap
  if (method == "const") {
    nu <- nu
  } else if (method == "rand") {
    nu <- runif(nrow(object), 1, nu)
  } else if (method == "wild") {
    nu <- runif(K, 1, nu)[fData(object)$knn] # cluster specific residual inflated noise
  }

  # code to generate replicates
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

  # should we had batch effects
  if (batch == "rand") {
    for (i in seq.int(numRep)){
      batch_effect <- rgamma(1, 1, 1) # sample batch effect
      g <- sample(c(seq.int(ncol(myobject))), size = 1) # sample which fraction to add the effect to
      exprs(lopitrep[[i]])[, g] <- exprs(lopitrep[[i]])[, g] + batch_effect
    }
  } else if (batch == "systematic") {
    g <- sample(c(seq.int(ncol(myobject))), size = 1) # sample which fraction to add the effect to
    batch_effect <- rgamma(1, 1, 1) # sample batch effect
    for (i in seq.int(numRep)) {
      exprs(lopitrep[[i]])[, g] <- exprs(lopitrep[[i]])[, g] + batch_effect
    }
  }

  # should we permute the fractions
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

  ## compute statistics for simulated replicates
  replicateStats <- lapply(lopitrep, function(x) meanOrganelle(x, fcol = fcol))

  # Simulate dynamics, make sure to use correct 
  for (i in perm1_names) {
    K <- length(getMarkerClasses(object, fcol = fcol))
    r <- sample.int(K, size = 1)
    for (j in seq.int(numRep/2 + 1, numRep)) {
      exprs(lopitrep[[j]])[i, ] <- replicateStats[[j]]$M[r, ] + 
        rnorm(ncol(object), mean = 0, sd = sqrt(replicateStats[[j]]$V[r, ]))
    }
  }
  
  # return simulate replicates and the names of dynamic translocations
  return(list(lopitrep = lopitrep, perm1_names = perm1_names))
  }
