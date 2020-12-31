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
    
    ch <- chains(params)[[n]]
    dfr <- as.data.frame(ch@nicheProb[[cond]][fname, , ])
    dfr_long <- data.frame(Organelle = rep(names(dfr), each = nrow(dfr)),
                           Probability = unlist(dfr, use.names = FALSE),
                           row.names = NULL,
                           stringsAsFactors = FALSE)
    gg2 <- ggplot(dfr_long,
                  aes(Organelle, Probability,
                      width = (Probability))) +
        geom_violin(aes(fill = Organelle), scale = scale, bw = bw)
    gg2 <- gg2 + theme_bw() +
        scale_fill_manual(values = pRoloc::getStockcol()[seq_len(nrow(dfr))]) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1),
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
                      cov.function = wendland.cov,
                      theta = 2,
                      derivative = 2,
                      k = 1,
                      cond = 1,
                      n = 1,
                      breaks = c(0.99, 0.95, 0.9, 0.85, 0.8, 0.75, 0.7),
                      aspect = 0.5) {
    
    stopifnot(class(object) == "MSnSet")
    stopifnot(class(params) == "bandleParams")
    
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
