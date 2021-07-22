##' These functions implement plotting functions for bandle objects
##' 
##' 
##' @title Generate a violin plot showing the probabilitiy of protein localisation
##'  to different organelles
##' @param params An instance of class `bandleParams`
##' @param fname The name of the protein to plot
##' @param cond Which conditions do we want to plot. Must be `1` or `2`. Default
##'  is `1`
##' @param n The chain from which we plot the probability distribution. Default is 1.
##' @param bw The bandwidth use in probability distribution smoothing of geom_violin
##'  Default is 0.05.
##' @param scale Scaling of geom_violin. Defaults to width.
##' @param trim trim parameter of geom_violin. Defaults to true. 
##' @return  returns a named vector of differential localisation probabilities
##' @md
##' 
##' @rdname bandle-plots
mcmc_plot_probs <- function(params,
                            fname,
                            cond = 1,
                            n = 1,
                            bw = 0.05,
                            scale = "width",
                            trim = TRUE) {
    
    stopifnot(class(params) == "bandleParams")
    stopifnot(class(fname) == "character")
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
##' @param cov.function The covariance function for the smoothing kernel. Defaults
##'  to wendland.cov
##' @param n The chain from which we plot the probability distribution. Default is 1.
##' @param theta The theta parameter of the wendland.cov. Defaults to 2.
##' @param derivative The derivative paramter of the wendland.cov. Defaults to 2.
##' @param k The k parameter of the wendland.cov 
##' @param breaks The levels at which to plot the contours. Defaults to c(0.99, 0.95, 0.9, 0.85, 0.8, 0.75, 0.7)
##' @param aspect The aspect ratio of the pca plots. Defaults to 0.5.
##' @param cond Which conditions do we want to plot. Must be `1` or `2`. Default
##'  is `1`
##' @return  returns a named vector of differential localisation probabilities
##' @md
##' 
##' @rdname bandle-plots

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
    
    stopifnot(class(object) == "MSnSet")
    stopifnot(class(params) == "bandleParams")
    organelle <- x <- y <- z <- level <- NULL
    
    if (is.null(cov.function)){
        cov.function <- fields::wendland.cov
    }
    
    ## Compute mean probabilities 
    ch <- chains(params)[[n]]
    probs <- apply(ch@nicheProb[[cond]], c(1, 3), mean)
    
    ## create allocation matrix for markers
    markerSet <- markerMSnSet(object, fcol = fcol)
    .probmat <- matrix(0, nrow = nrow(markerSet), ncol = length(getMarkerClasses(object)))
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
        coords[[j]] <- akima::interp(x = probs.lst.df$x[idxOrg],
                                     y = probs.lst.df$y[idxOrg],
                                     z = probs.lst.df$probability[idxOrg],
                                     extrap=FALSE, linear = TRUE, duplicate = TRUE) # interpolate onto appropriate grid
        coords[[j]]$z[is.na(coords[[j]]$z)] <- 0 # NaNs beyond data set to 0
        locations[[j]] <- cbind(rep(coords[[j]]$x, 40), rep(coords[[j]]$y, each = 40)) # get grid
        smoothedprobs <- fields::smooth.2d(coords[[j]]$z, x = locations[[j]],
                                           cov.function = cov.function,
                                           theta = theta,
                                           derivative = derivative, k = k) # spatial smoothing of probabilities
        # normalisation and formatting
        zmat <- matrix(smoothedprobs$z, ncol = 1)
        zmat <- zmat/max(zmat)
        df[[j]] <- data.frame(x = rep(smoothedprobs$x, 64), y = rep(smoothedprobs$y, each = 64), z = zmat)
    }
    # format data
    df.lst <- plyr::ldply(df, .fun = function(x) x, .id = "organelle") 
    df.lst <- df.lst %>%
        mutate(organelle = factor(organelle)) 
    K <- length(getMarkerClasses(object))
    cols <- getStockcol()[1:K] # get appropriate colours
    
    
    gg <- ggplot(
        data = df.lst,
        aes(x = x, y = y, z = z, color = organelle)) +
        coord_fixed() + 
        geom_contour(breaks = breaks, size = 1.2, aes(alpha = stat(level))) + 
        geom_point(alpha = 0) + 
        xlab(paste0(eigs[1])) + 
        ylab(paste0(eigs[2])) +
        scale_alpha(guide = "none") + 
        theme(legend.position = "right", 
              text = element_text(size = 12)) +
        scale_color_manual(values = cols) +
        scale_fill_manual(values = cols) +
        theme_minimal() + 
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              aspect.ratio = aspect,
              panel.border = element_rect(colour = "black", fill = NA, size = 1),
              plot.title = element_text(hjust = 0.5, size = 20),
              legend.text=element_text(size = 14)) +
        ggtitle(label = "Spatial variation of localisation probabilities") 
    gg
    return(gg)
}

# ##' @title Generate a circos/chord plot for visualising changes in localisation
# ##' between two conditions/datasets
# ##' @param params An instance of class `bandleParams` or `MSnSetList` 
# ##' @param fcol To be specified if input is a `MSnSetList`. `fcol` must be a 
# ##' `character` of length 2 specifying the feature column names for the first and second
# ##' `MSnSet` that contain the allocated class labels (see example below).
# ##' @param all
# ##' @param cols
# ##' @param ...
# ##' @return A `ggplot` object
# ##' @md
# ##' 
# ##' @rdname bandle-plots
# 
# chordplot <- function(params, 
#                        all = FALSE, 
#                        cols, 
#                        labels = TRUE, 
#                        ...) {
#     
#     stopifnot(inherits(params, "bandleParams"))
#     
#     # get results from params
#     res1 <- summaries(params)[[1]]@posteriorEstimates$bandle.allocation
#     res2 <- summaries(params)[[2]]@posteriorEstimates$bandle.allocation
#     
#     # get all possible allocation classes
#     cl1 <- colnames(summaries(params)[[1]]@bandle.joint)
#     cl2 <- colnames(summaries(params)[[2]]@bandle.joint)
#     fct.lev <- union(cl1, cl2)
#     res1_lev <- factor(res1, fct.lev)
#     res2_lev <- factor(res2, fct.lev)
#     
#     # create data frame of translocations
#     dat <- data.frame(x = res1_lev, y = res2_lev)
#     dat$z <- 1
#     datdf <- dat %>% group_by(x, y, .drop = FALSE) %>% 
#         dplyr:::summarise(count=sum(z), .groups = "keep")
#     if (!all) {
#         torm <- which(datdf$x == datdf$y)
#         datdf <- datdf[-torm, ]
#     }
#     df <- as.data.frame(datdf)
#     
#     # add colour scheme if not provided
#     if (missing(cols)) {
#         grid.col <- segcols <- setNames(getStockcol()[seq(fct.lev)], fct.lev)
#         if (length(fct.lev) > length(getStockcol()))
#             grid.col <- segcols <- setNames(rainbow(length(fct.lev)), fct.lev)
#     } else {
#         if (length(fct.lev) > length(cols))
#             stop(message("Not enough colours specified for subcellular classes"))
#         grid.col <- cols
#     }
#     
#     # create circos
#     par(mar = c(1,1,1,1)*12, cex = 0.6, xpd=NA)
#     circos.par(gap.degree = 4)
#     chordDiagram(df, annotationTrack = "grid", 
#                  preAllocateTracks = 1, 
#                  grid.col = grid.col,
#                  directional = 1, 
#                  direction.type = c("diffHeight", "arrows"), 
#                  link.arr.type = "big.arrow", ...)
#     
#     # annotate tracking regions and customise
#     if (labels) {
#         circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
#             xlim = get.cell.meta.data("xlim")
#             ylim = get.cell.meta.data("ylim")
#             sector.name = get.cell.meta.data("sector.index")
#             circos.text(mean(xlim), 
#                         ylim[1] + .1, 
#                         sector.name, 
#                         facing = "clockwise",
#                         niceFacing = TRUE, 
#                         adj = c(-0.5, 0.1),
#                         cex = 1,
#                         col=grid.col[sector.name],
#                         font = 2)
#             circos.axis(h = "top",
#                         labels.cex = .6,
#                         major.tick.length = 1,
#                         sector.index = sector.name,
#                         track.index = 2)
#         }, bg.border = NA)
#     } else {
#         circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
#             xlim = get.cell.meta.data("xlim")
#             ylim = get.cell.meta.data("ylim")
#             sector.name = get.cell.meta.data("sector.index")
#             circos.axis(h = "top",
#                         labels.cex = .6,
#                         major.tick.length = 1,
#                         sector.index = sector.name,
#                         track.index = 2)
#         }, bg.border = NA)
#     }
#     circos.clear()
# }

##' @title Generate an alluvial/sankey riverplot for visualising changes in localisation
##' between two conditions
##' @param params An instance of class `bandleParams` or a `MSnSetList` of length 2.
##' @param all A logical specifying whether to count all proteins or only show those that
##' have changed in location between conditions. Default is `FALSE`.
##' @param cols A list of colours to define the classes in the data
##' @param labels Choice of "adj", "repel" or "none". Default is "adj".
##' @return A `ggplot` object
##' @md
##' 
##' @rdname bandle-plots
##' @author Lisa M Breckels

riverplot <- function(params, 
                      all = FALSE, 
                      cols, 
                      labels = "adj"
                      ) {

    stopifnot(inherits(params, "bandleParams") | inherits(params, "MSnSetList"))
    
    # get results from params
    res1 <- summaries(params)[[1]]@posteriorEstimates$bandle.allocation
    res2 <- summaries(params)[[2]]@posteriorEstimates$bandle.allocation
    
    # get all possible allocation classes
    cl1 <- colnames(summaries(params)[[1]]@bandle.joint)
    cl2 <- colnames(summaries(params)[[2]]@bandle.joint)
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
    torm <- which(df$count == 0)
    df <- df[-torm, ]
    
    # add colour scheme if not provided
    if (missing(cols)) {
        grid.col <- segcols <- setNames(getStockcol()[seq(fct.lev)], fct.lev)
        if (length(fct.lev) > length(getStockcol()))
            grid.col <- segcols <- setNames(rainbow(length(fct.lev)), fct.lev)
    } else {
        if (length(fct.lev) > length(cols))
            stop(message("Not enough colours specified for subcellular classes"))
        grid.col <- cols
    }
    
    # set colours for alluvial plot (this is a little tricky as a specific 
    # ordering is required for ggalluvial)
    names(df) <- c("Condition1", "Condition2", "value")
    levs1 <- levels(df$Condition1) 
    levs2 <- levels(df$Condition2)
    # levs1 <- levs1[table(df$Condition1)!=0]
    # levs2 <- levs2[table(df$Condition2)!=0]
    res1 <- unique(df$Condition1)
    res2 <- unique(df$Condition2)
    cond1_cols <- grid.col[levs1[levs1 %in% res1]]
    cond2_cols <- grid.col[levs2[levs2 %in% res2]]
    columnCols <- c(cond1_cols, cond2_cols)
    stratCols <- c(rev(cond1_cols), rev(cond2_cols))
    
    # convert to long format
    df_expanded <- df[rep(row.names(df), df$value), ]
    df_expanded <- df_expanded %>%
        mutate(id = row_number()) %>%
        pivot_longer(-c(value, id), names_to = "Condition", values_to = "label")
    
    # plot alluvial diagram
    q <- ggplot(df_expanded, aes(x = Condition, stratum = label, alluvium = id, fill = label)) +
        geom_flow(width = 0) +
        scale_fill_manual(values = columnCols) +
        scale_color_manual(values = stratCols) +
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
    
    if (labels == "adj") {
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
    if (labels == "repel") {
        q <- q +
            ggrepel::geom_text_repel(
                aes(label = ifelse(after_stat(x) == 1, as.character(after_stat(stratum)), "")),
                stat = "stratum", size = 4, direction = "y", nudge_x = -.6) +
            ggrepel::geom_text_repel(
                aes(label = ifelse(after_stat(x)  == 2, as.character(after_stat(stratum)), "")),
                stat = "stratum", size = 4, direction = "y", nudge_x = .6)
    }
    q
    
}