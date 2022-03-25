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
              invisible(NULL)
          })

##' @rdname bandleParams
setMethod("show", "nicheParam",
          function(object) {
              cat("Object of class \"", class(object), "\"\n", sep = "")
              cat(" method:", object@method, "\n")
              cat(" Datasets included:", object@dataset, "\n")
              cat(" Replicate:", object@replicate, "\n")
              cat(" Number of components:", object@K, "\n")
              cat(" Number of samples:", object@D, "\n")
              invisible(NULL)
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
          function(x) length(summaries(x)))

##' @rdname bandleParams
setMethod("length", "nicheParams",
          function(x) length(x@params))

##' @rdname bandleParams
setMethod("length", "nicheParams",
          function(x) length(params(x)))


##'@param object An instance of appropriate class.
##'@rdname bandleParams
posteriorEstimates <- function(object) {
    stopifnot(inherits(object, "bandleSummary"))
    object@posteriorEstimates
}

##' @rdname bandleParams
setMethod("posteriorEstimates", "bandleSummary",
          function(object) object@posteriorEstimates)

##'@param object An instance of appropriate class.
##'@rdname bandleParams
summaries <- function(object) {
    object@summary@summaries
}

##'@param object An instance of appropriate class.
##'@rdname bandleParams
params <- function(object) {
    stopifnot(inherits(object, "nicheParams"))
    object@params
}

##'@param object An instance of appropriate class.
##'@rdname bandleParams
bandleJoint <- function(object) {
    stopifnot(inherits(object, "bandleSummary"))
    object@bandle.joint
}

##'@rdname bandleParams
setMethod("bandleJoint", "bandleSummary",
          function(object) object@bandle.joint)


##' @param x Object to be subset.
##' @param i An `integer()`. Should be of length 1 for `[[`.
##' @param j Missing.
##' @param drop Missing.
##' 
##' @md
##' @rdname bandleParams
setMethod("[[", "bandleChains",
          function(x, i, j = "missing", drop = "missing") x@chains[[i]])

##' @param x Object to be subset.
##' @param i An `integer()`. Should be of length 1 for `[[`.
##' @param j Missing.
##' @param drop Missing.
##' 
##' @md
##' @rdname bandleParams
setMethod("[[", "bandleParams",
          function(x, i, j = "missing", drop = "missing") params(x)[[i]])

##' @param x Object to be subset.
##' @param i An `integer()`. Should be of length 1 for `[[`.
##' @param j Missing.
##' @param drop Missing.
##' 
##' @md
##' @rdname bandleParams
setMethod("[", "bandleChains",
          function(x, i, j = "missing", drop = "missing") {
              if (any(i > length(x)))
                  stop("Index out of bounds. Only ", length(x), " chain(s) available.")
              x@chains <- x@chains[i]
              x
          })
##' @param x Object to be subset.
##' @param i An `integer()`. Should be of length 1 for `[[`.
##' @param j Missing.
##' @param drop Missing.
##' 
##' @md
##' @rdname bandleParams
setMethod("[", "bandleParams",
          function(x, i, j = "missing", drop = "missing") {
              if (any(i > length(x)))
                  stop("Index out of bounds. Only ", length(x), " chain(s) available.")
              x@chains <- chains(x)[i]
              x
          })
##' @param object of bandleChains
##' @rdname bandleParams
setMethod("show", "bandleChains",
          function(object) {
              cat(" Object of class \"", class(object), "\"\n", sep = "")
              cat(" Number of chains:", length(object), "\n")
              invisible(NULL)
          })
##' @param  object of class bandleSummaries
##' @rdname bandleParams
setMethod("show", "bandleSummaries",
          function(object) {
              cat(" Object of class \"", class(object), "\"\n", sep = "")
              cat(" Number of Summary(s):", length(object), "\n")
              invisible(NULL)
          })

##' @param x Object to be subset.
##' @param i An `integer()`. Should be of length 1 for `[[`.
##' @param j Missing.
##' @param drop Missing.
##' 
##' @md
##' @rdname bandleParams
setMethod("[[", "bandleSummaries",
          function(x, i, j = "missing", drop = "missing") x@summaries[[i]])

##' @param x Object to be subset.
##' @param i An `integer()`. Should be of length 1 for `[[`.
##' @param j Missing.
##' @param drop Missing.
##' 
##' @md
##' @rdname bandleParams
setMethod("[[", "bandleSummaries",
          function(x, i, j = "missing", drop = "missing") summaries(x)[[i]])

##' @param x Object to be subset.
##' @param i An `integer()`. Should be of length 1 for `[[`.
##' @param j Missing.
##' @param drop Missing.
##' 
##' @md
##' @rdname bandleParams
setMethod("[", "bandleSummaries",
          function(x, i, j = "missing", drop = "missing") {
              if (any(i > length(x)))
                  stop("Index out of bounds. Only ", length(x), " summary(s) available.")
              x@summaries <- x@summaries[i]
              x
          })
##' @param x Object to be subset.
##' @param i An `integer()`. Should be of length 1 for `[[`.
##' @param j Missing.
##' @param drop Missing.
##' 
##' @md
##' @rdname bandleParams
setMethod("[[", "nicheParams",
          function(x, i, j = "missing", drop = "missing") x@params[[i]])

##' @param x Object to be subset.
##' @param i An `integer()`. Should be of length 1 for `[[`.
##' @param j Missing.
##' @param drop Missing.
##' 
##' @md
##' @rdname bandleParams
setMethod("[[", "nicheParams",
          function(x, i, j = "missing", drop = "missing") params(x)[[i]])

##' @param x Object to be subset.
##' @param i An `integer()`. Should be of length 1 for `[[`.
##' @param j Missing.
##' @param drop Missing.
##' 
##' @md
##' @rdname bandleParams
setMethod("[", "nicheParams",
          function(x, i, j = "missing", drop = "missing") {
              if (any(i > length(x)))
                  stop("Index out of bounds. Only ", length(x), " param(s) available.")
              x@params <- x@params[i]
              x
          })

##' @param object object of class nicheParams.
##' @md
##' @rdname bandleParams
setMethod("show", "nicheParams",
          function(object) {
              cat(" Object of class \"", class(object), "\"\n", sep = "")
              cat(" Number of params:", length(object@params), "\n")
              invisible(NULL)
          })
