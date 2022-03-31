##' These functions implement plotting functions for bandle objects
##' 
##' 
##' @title Generate a violin plot showing the probabilitiy of protein
##'   localisation to different organelles
##' @param params An instance of class `bandleParams`
##' @param fname The name of the protein to plot
##' @param cond Which conditions do we want to plot. Must be `1` or `2`. Default
##'   is `1`
##' @param n The chain from which we plot the probability distribution. Default
##'   is 1.
##' @param bw The bandwidth use in probability distribution smoothing of
##'   geom_violin Default is 0.05.
##' @param scale Scaling of geom_violin. Defaults to width.
##' @param trim trim parameter of geom_violin. Defaults to true.
##' @return  returns a named vector of differential localisation probabilities
##' @md
##' 
##' @examples
##' library(pRolocdata)
##' data("tan2009r1")
##' set.seed(1)
##' tansim <- sim_dynamic(object = tan2009r1, 
##'                     numRep = 6L,
##'                    numDyn = 100L)
##' gpParams <- lapply(tansim$lopitrep, function(x) 
##' fitGPmaternPC(x, hyppar = matrix(c(0.5, 1, 100), nrow = 1)))
##' d1 <- tansim$lopitrep
##' control1 <- d1[1:3]
##' treatment1 <- d1[4:6]
##' mcmc1 <- bandle(objectCond1 = control1,
##'  objectCond2 = treatment1, gpParams = gpParams,
##'  fcol = "markers", numIter = 25L, burnin = 1L, thin = 2L,
##'  numChains = 2, BPPARAM = SerialParam(RNGseed = 1))
##' mcmc_plot_probs(params = mcmc1, fname = rownames(tan2009r1)[1])
##' 
##' @rdname bandle-plots-prob
mcmc_plot_probs <- function(params,
                            fname,
                            cond = 1,
                            n = 1,
                            bw = 0.05,
                            scale = "width",
                            trim = TRUE) {
    
    stopifnot(is(params, "bandleParams"))
    stopifnot(is(fname, "character"))
    Organelle <- Probability <- NULL
    
    ch <- chains(params)[[n]]
    dfr <- as.data.frame(ch@nicheProb[[cond]][fname, , ])
    dfr_long <- data.frame(Organelle = rep(names(dfr), each = nrow(dfr)),
                           Probability = unlist(dfr, use.names = FALSE),
                           row.names = NULL,
                           stringsAsFactors = FALSE)
    gg2 <- ggplot(dfr_long,
                  aes(Organelle, Probability),
                      width = Probability) +
        geom_violin(aes(fill = Organelle), scale = scale, bw = bw)
    gg2 <- gg2 + theme_bw() +
        scale_fill_manual(values = pRoloc::getStockcol()[seq_len(nrow(dfr))]) +
        theme(text = element_text(size=15),
              axis.text.x = element_text(angle = 45, hjust = 1),
              axis.title.x = element_blank())
    gg2 <- gg2 +
        ylab("Membership Probability") +
        ggtitle(paste0("Distribution of Subcellular Membership for Protein ", fname ))
    gg2 <- gg2 +
        theme(legend.position = "none")
    return(gg2)
}
##' @title Generate a PCA plot with smoothed probability contours
##' @param object An instance of class `MSnSet` to provide the pca coordinates
##' @param params An instance of class `bandleParams`
##' @param fcol Feature columns that defines the markers. Defaults to "markers".
##' @param dims The PCA dimensions to plot. Defaults to `c(1,2)`
##' @param cov.function The covariance function for the smoothing kernel.
##'   Defaults to wendland.cov
##' @param n The chain from which we plot the probability distribution. Default
##'   is 1.
##' @param theta The theta parameter of the wendland.cov. Defaults to 2.
##' @param derivative The derivative paramter of the wendland.cov. Defaults to
##'   2.
##' @param k The k parameter of the wendland.cov
##' @param breaks The levels at which to plot the contours. Defaults to c(0.99,
##'   0.95, 0.9, 0.85, 0.8, 0.75, 0.7)
##' @param aspect The aspect ratio of the pca plots. Defaults to 0.5.
##' @param cond Which conditions do we want to plot. Must be `1` or `2`. Default
##'   is `1`
##' @return  returns a named vector of differential localisation probabilities
##' @md
##'
##' @rdname bandle-plots-spatial

spatial2D <- function(object,
                      params, 
                      fcol = "markers",
                      dims = c(1, 2),
                      cov.function = NULL,
                      theta = 2,
                      derivative = 2,
                      k = 1,
                      cond = 1,
                      n = 1,
                      breaks = c(0.99, 0.95, 0.9, 0.85, 0.8, 0.75, 0.7),
                      aspect = 0.5) {
    
    stopifnot(is(object, "MSnSet"))
    stopifnot(is(params, "bandleParams"))
    organelle <- x <- y <- z <- level <- NULL
    
    if (is.null(cov.function)){
        cov.function <- fields::wendland.cov
    }
    
    ## Compute mean probabilities 
    ch <- chains(params)[[n]]
    probs <- apply(ch@nicheProb[[cond]], c(1, 3), mean)
    
    ## create allocation matrix for markers
    markerSet <- markerMSnSet(object, fcol = fcol)
    .probmat <- matrix(0, nrow = nrow(markerSet), 
                       ncol = length(getMarkerClasses(object)))
    .class <- fData(markerSet)[, fcol]
    for (j in seq_len(nrow(markerSet))) {
        ## give markers prob 1
        .probmat[j, as.numeric(factor(.class), seq(1, length(unique(.class))))[j]] <- 1
    }
    colnames(.probmat) <- getMarkerClasses(object)
    rownames(.probmat) <- rownames(markerSet)
    
    .probs <- rbind(probs, .probmat)
    
    # generate pca plot and create data from with probabilties
    .pca <- plot2D(object, dims = dims, plot = FALSE)
    probs <- data.frame(x = .pca[, 1],
                        y = .pca[, 2],
                        mcmc.prob = .probs[rownames(object),])
    colnames(probs) <- c(c("x", "y"), getMarkerClasses(object))
    eigs <- colnames(.pca)
    
    # put data in appropriate long format
    probs.lst <- list()
    for(j in getMarkerClasses(object)) {
        probs.lst[[j]] <- probs[, c("x", "y", j)]
        colnames(probs.lst[[j]]) <- c("x", "y", "probability")
    }
    probs.lst.df <- plyr::ldply(probs.lst, .fun = function(x) x, .id = "organelle")
    
    # Create storage
    coords <- list()
    locations <- list()
    df <- list()
    # Create appropriate spatial grid
    for (j in getMarkerClasses(object)) {
        idxOrg <- c(probs.lst.df$organelle == j)
        coords[[j]] <- interp::interp(x = probs.lst.df$x[idxOrg],
                                     y = probs.lst.df$y[idxOrg],
                                     z = probs.lst.df$probability[idxOrg],
                                     extrap=FALSE, linear = TRUE, 
                                     duplicate = "mean") 
                                     # interpolate onto appropriate grid
        coords[[j]]$z[is.na(coords[[j]]$z)] <- 0 # NaNs beyond data set to 0
        locations[[j]] <- cbind(rep(coords[[j]]$x, 40), 
                                rep(coords[[j]]$y, each = 40)) # get grid
        smoothedprobs <- fields::smooth.2d(coords[[j]]$z, x = locations[[j]],
                                           cov.function = cov.function,
                                           theta = theta,
                                           derivative = derivative, k = k) 
                                           # spatial smoothing of probabilities
        # normalisation and formatting
        zmat <- matrix(smoothedprobs$z, ncol = 1)
        zmat <- zmat/max(zmat)
        df[[j]] <- data.frame(x = rep(smoothedprobs$x, 64), 
                              y = rep(smoothedprobs$y, each = 64), 
                              z = zmat)
    }
    # format data
    df.lst <- plyr::ldply(df, .fun = function(x) x, .id = "organelle") 
    df.lst <- df.lst %>%
        mutate(organelle = factor(organelle)) 
    K <- length(getMarkerClasses(object))
    col <- getStockcol()[seq.int(K)] # get appropriate colours
    
    
    gg <- ggplot(
        data = df.lst,
        aes(x = .data$x, y = .data$y, z = .data$z, color = organelle)) +
        coord_fixed() + 
        geom_contour(breaks = breaks, size = 1.2, aes(alpha = stat(level))) + 
        geom_point(alpha = 0) + 
        xlab(paste0(eigs[1])) + 
        ylab(paste0(eigs[2])) +
        scale_alpha(guide = "none") + 
        theme(legend.position = "right", 
              text = element_text(size = 12)) +
        scale_color_manual(values = col) +
        scale_fill_manual(values = col) +
        theme_minimal() + 
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              aspect.ratio = aspect,
              panel.border = element_rect(colour = "black", fill = NA, size = 1),
              plot.title = element_text(hjust = 0.5, size = 20),
              legend.text=element_text(size = 14)) +
        ggtitle(label = "Spatial variation of localisation probabilities") 
    return(gg)
}

##' Produces a chord diagram (circos plot) or an alluvial plot (also known as a
##' Sankey diagram) to show changes in location between two conditions or
##' datasets.
##' @title Generates a chord diagram or alluvial plot for visualising changes in
##'   localisation between two conditions/datasets
##' @param params An instance of class \code{bandleParams} or an instance of
##'   class \code{MSnSetList} of length 2.
##' @param type A \code{character} specifying the type of visualisation to plot.
##'   One of \code{"alluvial"} (default) or \code{"chord"}.
##' @param all A logical specifying whether to count all proteins or only show
##'   those that have changed in location between conditions. Default is
##'   `FALSE`.
##' @param fcol If \code{params} is a \code{list} of \code{MSnSets}. Then
##'   \code{fcol} must be defined. This is a \code{character} vector of length 2
##'   to set different labels for each dataset. If only one label is specified,
##'   and the \code{character} is of length 1 then this single label will be
##'   used to identify the annotation column in both datasets.
##' @param col A list of colours to define the classes in the data. If not
##'   defined then the default \code{pRoloc} colours in \code{getStockCol()} are
##'   used.
##' @param labels Logical indicating whether to display class/organelle labels
##'   for the chord segments or alluvial stratum. Default is \code{TRUE}.
##' @param labels.par If \code{type} is \code{"alluvial"}. Label style can be
##'   specified as one of \code{"adj"}, \code{"repel"}. Default is \code{"adj"}.
##' @param spacer A `numeric`. Default is 4. Controls the white space around the
##'   circos plotting region.
##' @param cex Text size. Default is 1.
##' @param ... Additional arguments passed to the `chordDiagram` function.
##' @return Returns a directional circos/chord diagram showing the translocation
##'   of proteins between conditions. If \code{type = "alluvial"} ouput is a
##'   \code{ggplot} object.
##' @rdname bandle-plots-translocations

plotTranslocations <- function(params,
                               type = "alluvial",
                               all = FALSE,
                               fcol,
                               col,
                               labels = TRUE,
                               labels.par = "adj", 
                               cex = 1,
                               spacer = 4,
                               ...) {
    
    stopifnot(inherits(params, "bandleParams") | inherits(params, "list"))
    Condition <- x <- y <- z <- stratum <- value <- NULL

    # check method is one of chord or alluvial
    type <- match.arg(type, c("alluvial", "chord"))
    labels.par <- match.arg(labels.par, c("adj", "repel"))
    
    # get results from params
    if (inherits(params, "bandleParams")) {
        res1 <- summaries(params)[[1]]@posteriorEstimates$bandle.allocation
        res2 <- summaries(params)[[2]]@posteriorEstimates$bandle.allocation
        cl1 <- colnames(summaries(params)[[1]]@bandle.joint)
        cl2 <- colnames(summaries(params)[[2]]@bandle.joint)
    }
    if (inherits(params, "list")) {
        stopifnot(unlist(lapply(params, function(z) inherits(z, "MSnSet"))))
        params <- commonFeatureNames(params) ## keep only intersection between datasets
        params <- list(params[[1]], params[[2]])
        if (missing(fcol)) stop(message("Missing fcol, please specify feature columns"))
        if (length(fcol) == 1) {
            fcol <- rep(fcol, 2)
            message(
              c("------------------------------------------------",
                "\nIf length(fcol) == 1 it is assumed that the",
                "\nsame fcol is to be used for both datasets",
                "\nsetting fcol = c(", fcol[1], ",", fcol[2],")",
                "\n----------------------------------------------")
            )
        }
        for (i in seq(fcol)) {
            if (!is.null(fcol[i]) && !fcol[i] %in% fvarLabels(params[[i]]))
                stop("No fcol found in MSnSet, please specify a valid fcol ",
                     immediate. = TRUE)
        }
        res1 <- fData(params[[1]])[, fcol[1]]
        res2 <- fData(params[[2]])[, fcol[2]]
        cl1 <- names(table(res1))
        cl2 <- names(table(res2))
    }
    fct.lev <- union(cl1, cl2)
    res1_lev <- factor(res1, fct.lev)
    res2_lev <- factor(res2, fct.lev)
 
    # create data frame of translocations
    dat <- data.frame(x = res1_lev, y = res2_lev)
    dat$z <- 1
    datdf <- dat %>% group_by(x, y, .drop = FALSE) %>%
        dplyr:::summarise(count=sum(z), .groups = "keep")
    if (!all) {
        torm <- which(datdf$x == datdf$y)
        datdf <- datdf[-torm, ]
    }
    df <- as.data.frame(datdf)
    
    # add colour scheme if not provided
    if (missing(col)) {
        setStockcol(NULL)
        grid.col <- setNames(getStockcol()[seq(fct.lev)], fct.lev)
        if (length(fct.lev) > length(getStockcol()))
            grid.col <- setNames(rainbow(length(fct.lev)), fct.lev)
    } else {
        if (length(fct.lev) > length(col)) 
           stop(message("Not enough colours specified for subcellular classes"))
        grid.col <- col
        if (is.null(names(col))) {
          names(grid.col) <- fct.lev
          warning(message("Please specify a named character of colours\n to ensure matching of correct colours to subcellular classes "))
        } 
    }
    
    # ------------chord/circos plot
    if (type == "chord") {
        # clear any previous plot
        circos.clear()

        # create circos
        par(mar = c(1,1,1,1)*spacer, cex = 1, xpd = NA)
        circos.par(gap.degree = 4)
        chordDiagram(df, annotationTrack = "grid",
                     preAllocateTracks = 1,
                     grid.col = grid.col,
                     directional = 1,
                     direction.type = c("diffHeight", "arrows"),
                     link.arr.type = "big.arrow", ...)
        
        # annotate tracking regions and customise
        if (labels) {
            circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
                xlim = get.cell.meta.data("xlim")
                ylim = get.cell.meta.data("ylim")
                sector.name = get.cell.meta.data("sector.index")
                circos.text(CELL_META$xcenter, 
                            ylim[1] + cm_h(2), 
                            sector.name, 
                            facing = "clockwise",
                            niceFacing = TRUE, 
                            adj = c(0, 0.5),
                            cex = cex,
                            col=grid.col[sector.name],
                            font = 2)
                circos.axis(h = "bottom",
                            labels.cex = .6,
                            sector.index = sector.name
                )
            }, bg.border = NA)
        } else {
            circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
                xlim = get.cell.meta.data("xlim")
                ylim = get.cell.meta.data("ylim")
                sector.name = get.cell.meta.data("sector.index")
                circos.axis(h = "top",
                            labels.cex = .6,
                            major.tick.length = 1,
                            sector.index = sector.name,
                            track.index = 2)
            }, bg.border = NA)
        }
        circos.clear()
    }
    
    # -----------alluvial plot
    if (type == "alluvial") {
        
        # remove zero counts
        torm <- which(df$count == 0)
        df <- df[-torm, ]

        # set colours for alluvial plot (this is a little tricky as a specific 
        # ordering is required for ggalluvial)
        names(df) <- c("Condition1", "Condition2", "value")
        levs1 <- levels(df$Condition1) 
        levs2 <- levels(df$Condition2)
        # levs1 <- levs1[table(df$Condition1)!=0]
        # levs2 <- levs2[table(df$Condition2)!=0]
        res1 <- unique(df$Condition1)
        res2 <- unique(df$Condition2)
        cond1_col <- grid.col[levs1[levs1 %in% res1]]
        cond2_col <- grid.col[levs2[levs2 %in% res2]]
        columncol <- c(cond1_col, cond2_col)
        stratcol <- c(rev(cond1_col), rev(cond2_col))
        
        # convert to long format
        df_expanded <- df[rep(row.names(df), df$value), ]
        df_expanded <- df_expanded %>%
            mutate(id = row_number()) %>%
            pivot_longer(-c(value, id), names_to = "Condition", values_to = "label")
        
        # plot alluvial diagram
        q <- ggplot(df_expanded, aes(x = .data$Condition, stratum = .data$label,
                                     alluvium = id, fill = .data$label)) +
            geom_flow(width = 0) +
            scale_fill_manual(values = columncol) +
            scale_color_manual(values = stratcol) +
            geom_stratum(width = 1/8, color = "white") +
            scale_x_discrete(
                expand = c(.25, .25)
            ) +
            scale_y_continuous(breaks = NULL) +
            theme_minimal() +
            theme(
                axis.ticks.y = element_blank(),
                axis.text.y = element_blank(),
                axis.text.x = element_text(size=12),
                panel.grid.major.y = element_blank(),
                panel.grid.major.x = element_blank(),
                axis.title.x=element_blank()
            ) +
            theme(legend.position = "none") +
            ylab(NULL)
        
        # add labels to alluvial stratum
        if (labels == "TRUE") {
            if (labels.par == "adj") {
                q <- q +
                    geom_text(
                        aes(
                            label = after_stat(stratum),
                            hjust = ifelse(Condition == "Condition1", 1, 0),
                            x = as.numeric(factor(Condition)) + .075 * ifelse(Condition == "Condition1", -1, 1),
                            color = after_stat(stratum)
                        ),
                        stat = "stratum", fontface = "bold", size = 4
                    )
            }
            if (labels.par == "repel") {
                q <- q +
                    ggrepel::geom_text_repel(
                        aes(label = ifelse(after_stat(x) == 1, as.character(after_stat(stratum)), "")),
                        stat = "stratum", size = 4, direction = "y", nudge_x = -.6) +
                    ggrepel::geom_text_repel(
                        aes(label = ifelse(after_stat(x)  == 2, as.character(after_stat(stratum)), "")),
                        stat = "stratum", size = 4, direction = "y", nudge_x = .6)
            }
            labels.par <- match.arg(labels.par, c("repel", "adj"))
        } 
        return(q)
    }
    
}
##' Produces a rank plot to analyse convergence of MCMC algorithm
##' @title Generates a histogram of ranks (a rank plot) for convergence
##' @param params An instance of class \code{bandleParams}
##' @return Returns the ranks of the number of outliers in each chain. The
##' side effect returns rank plots. Number of rank plots is equal to the number
##' of chains
##' @md
##' 
##'  
##' @rdname bandle-plots-convergence
plotConvergence <- function(params){
    
    stopifnot("params must be an object of
              class bandleParams"=is(params, "bandleParams"))
    
    n <- as.numeric(seq(params@chains@chains[[1]]@n))
    out <- vapply(params@chains@chains, 
                  function(x) colSums(x@outlier$cond11), n)
    toplot <- apply(out, 1, rank)
    sapply(seq.int(nrow(toplot)), function(i) 
        hist(toplot[i,], main = "convergence rank plot", xlab = "rank", 
             col = alpha(getStockcol()[i], 0.7)))
    return(toplot)
}

##' Produces a table summarising differential localisation results
##' @title Generates a table for visualising changes in
##'   localisation between two conditions/datasets
##' @param params An instance of class \code{bandleParams} or an instance of
##'   class \code{MSnSetList} of length 2.
##' @param all A logical specifying whether to count all proteins or only show
##'   those that have changed in location between conditions. Default is
##'   `FALSE`.
##' @param fcol If \code{params} is a \code{list} of \code{MSnSets}. Then
##'   \code{fcol} must be defined. This is a \code{character} vector of length 2
##'   to set different labels for each dataset. If only one label is specified,
##'   and the \code{character} is of length 1 then this single label will be
##'   used to identify the annotation column in both datasets.
##' @return Returns a summary table of translocations of proteins between conditions. 
##' @rdname bandle-plots-translocations-table
plotTable <- function(params,
                      all = FALSE,
                      fcol) {
    
    stopifnot(inherits(params, "bandleParams") | inherits(params, "list"))
    Condition <- x <- y <- z <- stratum <- value <- grid.col <- NULL
  
    # get results from params
    if (inherits(params, "bandleParams")) {
        res1 <- summaries(params)[[1]]@posteriorEstimates$bandle.allocation
        res2 <- summaries(params)[[2]]@posteriorEstimates$bandle.allocation
        cl1 <- colnames(summaries(params)[[1]]@bandle.joint)
        cl2 <- colnames(summaries(params)[[2]]@bandle.joint)
    }
    if (inherits(params, "list")) {
        stopifnot(unlist(lapply(params, function(z) inherits(z, "MSnSet"))))
        params <- commonFeatureNames(params) ## keep only intersection between datasets
        params <- list(params[[1]], params[[2]])
        if (missing(fcol)) stop(message("Missing fcol, please specify feature columns"))
        if (length(fcol) == 1) {
            fcol <- rep(fcol, 2)
            message(
                c("------------------------------------------------",
                  "\nIf length(fcol) == 1 it is assumed that the",
                  "\nsame fcol is to be used for both datasets",
                  "\nsetting fcol = c(", fcol[1], ", ", fcol[2],")",
                  "\n----------------------------------------------")
            )
        }
        for (i in seq(fcol)) {
            if (!is.null(fcol[i]) && !fcol[i] %in% fvarLabels(params[[i]]))
                stop("No fcol found in MSnSet, please specify a valid fcol ",
                     immediate. = TRUE)
        }
        res1 <- fData(params[[1]])[, fcol[1]]
        res2 <- fData(params[[2]])[, fcol[2]]
        cl1 <- names(table(res1))
        cl2 <- names(table(res2))
    }
    fct.lev <- union(cl1, cl2)
    res1_lev <- factor(res1, fct.lev)
    res2_lev <- factor(res2, fct.lev)
    
    # create data frame of translocations
    dat <- data.frame(x = res1_lev, y = res2_lev)
    dat$z <- 1
    datdf <- dat %>% group_by(x, y, .drop = FALSE) %>%
        dplyr:::summarise(count=sum(z), .groups = "keep")
    if (!all) {
        torm <- which(datdf$x == datdf$y)
        datdf <- datdf[-torm, ]
    }
    df <- as.data.frame(datdf)
    
        # remove zero counts
        torm <- which(df$count == 0)
        df <- df[-torm, ]
        
        # set colours for alluvial plot (this is a little tricky as a specific 
        # ordering is required for ggalluvial)
        names(df) <- c("Condition1", "Condition2", "value")
        levs1 <- levels(df$Condition1) 
        levs2 <- levels(df$Condition2)
        # levs1 <- levs1[table(df$Condition1)!=0]
        # levs2 <- levs2[table(df$Condition2)!=0]
        res1 <- unique(df$Condition1)
        res2 <- unique(df$Condition2)
        cond1_col <- grid.col[levs1[levs1 %in% res1]]
        cond2_col <- grid.col[levs2[levs2 %in% res2]]
        columncol <- c(cond1_col, cond2_col)
        stratcol <- c(rev(cond1_col), rev(cond2_col))

        # return table
        return(df = df)
}
