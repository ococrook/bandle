##' @slot chains `list()` containing the individual full MCMC chain
##'     results in an `bandleChains` instance. Each element must be a
##'     valid `bandleChain` instance.
##' @md
##' @rdname bandleParams
.bandleChains <- setClass("bandleChains",
                          slots = c(chains = "list"),
                          validity = function(object) {
                              msg <- validMsg(NULL, NULL)
                              cls <- vapply(object@chains,
                                            function(x) inherits(x, "bandleChain"),
                                            logical(1))
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
##' @md
##' @rdname bandleParams    
.bandleSummaries <- setClass("bandleSummaries",
                             slots = c(summaries = "list"),
                             validity = function(object) {
                                 msg <- validMsg(NULL, NULL)
                                 cls <- vapply(object@summaries,
                                               function(x) inherits(x, "bandleSummary"),
                                               logical(1))
                                 if (!all(cls))
                                     msg <- validMsg(msg, "Not all items are bandleSummary(s).")
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
##' 
##' @slot method A `character()` storing the bandle method name
##' @slot priors A `list()` with the priors for the parameters
##' @slot seed An `integer()` with the random number generation seed.
##' @slot summary Object of class `bandleSummary` the summarised MCMC results
##' available in the `bandleParams` instance.
##' @slot chains Object of class `bandleChains` containing the full MCMC results
##' in the `bandleParams` instance
##' @return An object of class `bandleParams` which stores the main results
##' for the analysis when using bandle
##' 
##' @md
##' @rdname bandleParams
.bandleParams <- setClass("bandleParams",
                          slots = c(method = "character",
                                    priors = "list",
                                    seed = "integer",
                                    summary = "bandleSummaries",
                                    chains = "bandleChains"))

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
##' @rdname bandleParams
.nicheParam <- setClass("nicheParam",
                        slots = c(dataset = "character",
                                  replicate = "integer",
                                  K = "integer",
                                  D = "integer",
                                  method = "character",
                                  params = "list"),
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

##' @slot params `list()` containing the individual `nicheParam` objects
##'     results in an `bandleParams` instance. Each element must be a
##'     valid `bandleParam` instance.
##' @md
##' @rdname bandleParams
.nicheParams <- setClass("nicheParams",
                         slots = c(params = "list"),
                         validity = function(object) {
                             msg <- validMsg(NULL, NULL)
                             cls <- vapply(object@params,
                                           function(x) inherits(x, "nicheParam"),
                                           logical(1))
                             if (!all(cls))
                                 msg <- validMsg(msg, "Not all items are nicheParams.")
                             if (is.null(msg)) TRUE
                             else msg
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
                                   replicates = "integer",
                                   n = "integer",
                                   K = "integer",
                                   N = "integer",
                                   weights = "array",
                                   epsilons = "matrix",
                                   niche = "list",
                                   nicheProb = "list",
                                   outlier = "list",
                                   nicheParams = "nicheParams"),
                         validity = function(object) {
                             msg <- validMsg(NULL, NULL)
                             N <- object@N
                             n <- object@n
                             K <- object@K
                             if (!identical(nrow(object@niche[[1]]), N))
                                 msg <- validMsg(msg, "Wrong number of proteins in niche")
                             if (!identical(nrow(object@outlier[[1]]), N))
                                 msg <- validMsg(msg, "Wrong number of proteins in outlier")
                             if (!identical(ncol(object@niche[[1]]), n))
                                 msg <- validMsg(msg, "Wrong number of iterations in niche")
                             if (!identical(ncol(object@outlier[[1]]), n))
                                 msg <- validMsg(msg, "Wrong number of iterations in outlier")
                             if (!identical(rownames(object@niche), rownames(object@nicheProb)))
                                 msg <- validMsg(msg, "Component rownames don't match")
                             if (!identical(rownames(object@niche), rownames(object@nicheProb)))
                                 msg <- validMsg(msg, "Outlier rownames don't match")
                             if (!identical(rownames(object@outlier), rownames(object@niche)))
                                 msg <- validMsg(msg, "Proteins don't match between niche and outlier")
                             if (!identical(dim(object@nicheProb[[1]])[3], K))
                                 msg <- validMsg(msg, "Wrong number of components in niche probability")
                             if (is.null(msg)) TRUE
                             else msg
                         })