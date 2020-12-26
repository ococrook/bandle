##' @slot chains `list()` containing the individual full MCMC chain
##'     results in an `bandleChains` instance. Each element must be a
##'     valid `bandleChain` instance.
##' @md
##' @rdname bandleParams
.bandleChains <- setClass("bandleChains",
                        slots = c(chains = "list"),
                        validity = function(object) {
                            msg <- validMsg(NULL, NULL)
                            cls <- sapply(object@chains,
                                          function(x) inherits(x, "bandleChain"))
                            if (!all(cls))
                                msg <- validMsg(msg, "Not all items are bandleChains.")
                            if (is.null(msg)) TRUE
                            else msg
                        })

##' @slot posteriorEstimates A `DataFrame` documenting the posteriors
##'  in an `bandleSummary` instance
##' @slot diagnostics A `matrix` of dimensions 1 by 2 containing the
##'     `bandleSummary` diagnostics.
##' @slot bandle.joint A `matrix` of dimensions N by K storing the joint
##'     probability in an `bandleSummary` instance for each of the first condition
##' @md
##' @rdname bandleParams
.bandleSummary <- setClass("bandleSummary",
                         slots = c(posteriorEstimates = "DataFrame",
                                   diagnostics = "matrix",
                                   bandle.joint = "matrix"))


##' @slot chains `list()` containing the individual bandle Summaries for 
##' different conditions results in an `bandleSummaries` instance. Each element must be a
##'     valid `bandleSummary` instance.
.bandleSummaries <- setClass("bandleSummaries",
                            slots = c(summaries = "list"),
                            validity = function(object) {
                                msg <- validMsg(NULL, NULL)
                                cls <- sapply(object@summaries,
                                              function(x) inherits(x, "bandleSummary"))
                                if (!all(cls))
                                    msg <- validMsg(msg, "Not all items are bandleSummaries.")
                                if (is.null(msg)) TRUE
                                else msg
                            })

##' The `bandleParams` infrastructure is used to store and process MCMC results for
##' bandle model from Crook et al 2021
##'
##' Objects of the `bandleParams` class are created with the `bandle()` function
##' These objects store the *priors* for the model and the results of the MCMC
##' chains, which themselves are stored as an instance of class `bandleChains` and
##' can be accessed with the `chains()` function. A summary of the `bandleChains`
##' (or class `bandleSummary`) can be further computed with the `bandleProcess`
##' function.
##' 
##' see the *bandle* vignette for examples
##' @title Infrastructure to to store and process MCMC results
##' @slot method A `character()` storing the bandle method name
##' @slot priors A `list()` with the priors for the parameters
##' @slot seed An `integer()` with the random number generation seed.
##' @slot summary Object of class `bandleSummary` the summarised MCMC results
##' available in the `bandleParams` instance.
##' @slot chains Object of class `bandleChains` containing the full MCMC results
##' in the `bandleParams` instance
##' @md
##' @rdname `bandleParams`
.bandleParams <- setClass("bandleParams",
                          slots = c(method = "character",
                                    priors = "list",
                                    seed = "integer",
                                    summary = "bandleSummary",
                                    chains = "bandleChains"))

##'@param object An instance of appropriate class.
##'@rdname bandleParams
chains <- function(object) {
    stopifnot(inherits(object, "bandleParams"))
    object@chains
}

##' @md
##' @rdname bandleParams
setMethod("show", "bandleParams",
          function(object) {
              cat("Object of class \"", class(object), "\"\n", sep = "")
              cat("Method:", object@method, "\n")
              cat("Number of chains:", length(object@chains), "\n")
              if (nrow(object@summary@differential.localisation))
                  cat("Summary available\n")
              invisible(NULL)
          })

##' @slot datset `character` indicating which dataset i.e control or treatment
##' @slot replicate `integer` an integer indicating which replicate
##' @slot K `integer(1)` indicating the number of components.
##' @slot D `integer(1)` indicating the number of samples.
##' @slot method `character(1)` defining the method used. Currently
##'     `bandle`
##' @slot mk `matrix(K, D)`
##' @slot lambdak `numeric(K)`
##' @slot nuk `numeric(K)`
##' @slot sk `array(K, D, D)`
##' @md
##' @noRd
.nicheParam <- setClass("nicheParam",
                            slots = c(dataset = "character",
                                      relicate = "integer",
                                      K = "integer",
                                      D = "integer",
                                      method = "character"),
                            prototype = prototype(
                                method = "bandle"
                            ),
                            validity = function(object) {
                                msg <- validMsg(NULL, NULL)
                                K <- object@K
                                D <- object@D
                                if (object@method != "bandle")
                                    msg <- validMsg(msg, "Wrong method")
                                if (is.null(msg)) TRUE
                                else msg
                            })

##' @rdname MCMCParams
setMethod("show", "nicheParam",
          function(object) {
              cat("Object of class \"", class(object), "\"\n", sep = "")
              cat(" method:", object@method, "\n")
              cat(" Datasets included:", object@dataset)
              cat(" Number of replicates:", object@replicate)
              cat(" Number of components:", object@K, "\n")
              cat(" Number of samples:", object@D, "\n")
              invisible(NULL)
          })


##' @title Container for a single bandle chain results
##'
##' @slot dataset `character` indicating the dataset usaully control or treatment
##' @slot replicate `integer` indicating the number of dataset replicate
##' @slot n `integer(1)` indicating the number of MCMC interactions.
##'     Stored in an `bandleChain` instance.
##' @slot K `integer(1)` indicating the number of components. Stored
##'     in an `bandleChain` instance.
##' @slot N `integer(1)` indicating the number of proteins. Stored in
##'     an `bandleChain` instance.
##' @slot niche `matrix(N, n)` component allocation results of an
##'     `bandleChain` instance.
##' @slot nicheProb `matrix(N, n, K)` component allocation
##'     probabilities of an `bandleChain` instance.
##' @slot outlier `matrix(N, n)` outlier allocation results.
##' @slot outlierProb `matrix(N, n, 2)` outlier allocation
##'     probabilities of an `bandleChain` instance.
##' @md
##' @rdname bandleParams
.bandleChain <- setClass("bandleChain",
                       slots = c(dataset = "character",
                                 replicate = "integer",
                                 n = "integer",
                                 K = "integer",
                                 N = "integer",
                                 niche = "matrix",
                                 nicheProb = "array",
                                 outlier = "matrix",
                                 nicheParam = "nicheParam"),
                       validity = function(object) {
                           msg <- validMsg(NULL, NULL)
                           N <- object@N
                           n <- object@n
                           K <- object@K
                           if (!identical(nrow(object@niche), N))
                               msg <- validMsg(msg, "Wrong number of proteins in niche")
                           if (!identical(nrow(object@outlier), N))
                               msg <- validMsg(msg, "Wrong number of proteins in outlier")
                           if (!identical(ncol(object@niche), n))
                               msg <- validMsg(msg, "Wrong number of iterations in niche")
                           if (!identical(ncol(object@outlier), n))
                               msg <- validMsg(msg, "Wrong number of iterations in outlier")
                           if (!identical(rownames(object@niche), rownames(object@nicheProb)))
                               msg <- validMsg(msg, "Component rownames don't match")
                           if (!identical(rownames(object@niche), rownames(object@nicheProb)))
                               msg <- validMsg(msg, "Outlier rownames don't match")
                           if (!identical(rownames(object@outlier), rownames(object@niche)))
                               msg <- validMsg(msg, "Proteins don't match between niche and outlier")
                           if (!identical(dim(object@nicheProb)[3], K))
                               msg <- validMsg(msg, "Wrong number of components in niche probability")
                           if (is.null(msg)) TRUE
                           else msg
                       })

##' @rdname bandleParams
setMethod("show", "bandleChain",
          function(object) {
              cat(" Object of class \"", class(object), "\"\n", sep = "")
              cat(" Number of components:", object@K, "\n")
              cat(" Number of proteins:", object@N, "\n")
              cat(" Number of iterations:", object@n, "\n")
              invisible(NULL)
          })

##' @rdname bandleParams
setMethod("length", "bandleChains",
          function(x) length(x@chains))

##' @rdname bandleParams
setMethod("length", "bandleParams",
          function(x) length(chains(x)))

##' @rdname bandleParams
setMethod("length", "bandleSummaries",
          function(x) length(Summaries(x)))

##'@param object An instance of appropriate class.
##'@rdname bandleParams
summaries <- function(object) {
    stopifnot(inherits(object, "bandleParams"))
    object@summaries
}

##' @param x Object to be subset.
##' @param i An `integer()`. Should be of length 1 for `[[`.
##' @param j Missing.
##' @param drop Missing.
##' @md
##' @rdname bandleParams
setMethod("[[", "bandleChains",
          function(x, i, j = "missing", drop = "missing") x@chains[[i]])

##' @rdname bandleParams
setMethod("[[", "bandleParams",
          function(x, i, j = "missing", drop = "missing") chains(x)[[i]])

##' @rdname bandleParams
setMethod("[", "bandleChains",
          function(x, i, j = "missing", drop = "missing") {
              if (any(i > length(x)))
                  stop("Index out of bounds. Only ", length(x), " chain(s) available.")
              x@chains <- x@chains[i]
              x
          })

##' @rdname bandleParams
setMethod("[", "bandleParams",
          function(x, i, j = "missing", drop = "missing") {
              if (any(i > length(x)))
                  stop("Index out of bounds. Only ", length(x), " chain(s) available.")
              x@chains <- chains(x)[i]
              x
          })

##' @rdname bandleParams
setMethod("show", "bandleChains",
          function(object) {
              cat(" Object of class \"", class(object), "\"\n", sep = "")
              cat(" Number of chains:", length(object), "\n")
              invisible(NULL)
          })

##' @rdname bandleParams
setMethod("show", "bandleSummaries",
          function(object) {
              cat(" Object of class \"", class(object), "\"\n", sep = "")
              cat(" Number of chains:", length(object), "\n")
              invisible(NULL)
          })

##' @param x Object to be subset.
##' @param i An `integer()`. Should be of length 1 for `[[`.
##' @param j Missing.
##' @param drop Missing.
##' @md
##' @rdname bandleParams
setMethod("[[", "bandleSummaries",
          function(x, i, j = "missing", drop = "missing") x@summaries[[i]])

##' @rdname bandleParams
setMethod("[[", "bandleSummaries",
          function(x, i, j = "missing", drop = "missing") summaries(x)[[i]])

##' @rdname bandleParams
setMethod("[", "bandleSummaries",
          function(x, i, j = "missing", drop = "missing") {
              if (any(i > length(x)))
                  stop("Index out of bounds. Only ", length(x), " summary(s) available.")
              x@summaries <- x@summaries[i]
              x
          })



