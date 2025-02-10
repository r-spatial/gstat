#' Parallelized vgmArea
#'
#' This function provides a parallelized analogue of
#' \code{\link[gstat]{vgmArea}}.
#'
#' @param cl a cluster object created by \code{\link{parallel}} or
#'   \code{\link{snow}}
#' @param x a SpatialPolygons object
#' @param y a SpatialPolygons or SpatialPoints object
#' @param vgm a variogram model created by \code{\link[gstat]{vgm}}
#' @param ndiscr number of points into which SpatialPolygons will be discretized
#' @param covariance logical; if true returns covariance values, otherwise
#'   return semivariance values.
#' @param ... additional arguments passed to \code{\link[gstat]{vgmArea}}
#' @export
parVgmArea <- function(x, y = x, vgm, nnodes=NULL, ndiscr = 16, covariance = TRUE, ...) {
    if(is.null(nnodes)) { ## stop('parVgmArea must have nnodes argument')
        return(gstat::vgmArea(x, y, vgm, ndiscr, covariance, ...))
    }
    if(gridded(x)) x <- as(x, "SpatialPolygons")
    if(gridded(y)) y <- as(y, "SpatialPolygons")

    stopifnot(is(x, "SpatialPolygons") || is(x, "SpatialPoints"))
    stopifnot(is(y, "SpatialPolygons") || is(y, "SpatialPoints"))
    stopifnot(is(vgm, "variogramModel"))

    nx <- length(x)
    ny <- length(y)

    ## initialize cluster
    cl <- parallel::makeCluster(nnodes, ...)
    parallel::clusterEvalQ(cl, {
        library(gstat)
        library(sp)
    })
    parallel::clusterExport(
                  cl, c('x', 'y', 'vgm', 'ndiscr', 'covariance', 'nx', 'ny'), envir=environment())

    V <- parallel::parLapply(cl, 1:nx, function(i) { # outer loop
        if (is(x, "SpatialPolygons")) {
            px <- sp::spsample(x[i,], ndiscr, "regular", offset = c(.5,.5))
        } else px <- x[i,]
        sapply(1:ny, function(j) { # inner loop
            if (is(y, "SpatialPolygons")) {
                py <- sp::spsample(y[j,], ndiscr, "regular", offset = c(.5,.5))
            } else py <- y[j,]
            D <- sp::spDists(px, py)
            D[D == 0] <- 1e-10
            mean(gstat::variogramLine(vgm, dist_vector = D, covariance = covariance))
        })
    })

    ## clean up and return
    parallel::stopCluster(cl)
    do.call(rbind, V)
}
