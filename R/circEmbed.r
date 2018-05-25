## Benedikt Gräler (52North), 2018-05-25:
## circulant embedding following: Davies, Tilman M., and David Bryant. "On
## circulant embedding for Gaussian random fields in R." Journal of Statistical
## Software 55.9 (2013): 1-21.
## See i.e. the suplementary files at (retreived 2018-05-25): 
## https://www.jstatsoft.org/index.php/jss/article/downloadSuppFile/v055i09/v55i09.R

library(spatstat)
data("chorley")

library(sp)
sppix <- SpatialPixels()
chorley$window
spgrid <- SpatialGrid(GridTopology(cellcentre.offset = c(chorley$window$xrange[1],
                                                         chorley$window$yrange[1]), 
                                   cellsize = c(diff(chorley$window$xrange)/29,
                                                diff(chorley$window$yrange)/29), 
                                   cells.dim = c(29,29)))

plot(spgrid, axes=TRUE)

sp:::plot.Spatial

# extend the grid 
# input: SpatialGrid/SpatialPixels
# output: SpatialGrid

fullgrid(spgrid)

grid <- spgrid

ceExtGrid <- function(grid, ext=2) {
  stopifnot(gridded(grid))
  
  if (!fullgrid(grid)) {
    warning("SpatialPixels have been coerced to a SpatialGrid")
    fullgrid(grid) <- TRUE
  }
  
  SpatialGrid(GridTopology(grid@grid@cellcentre.offset,
                           grid@grid@cellsize,
                           grid@grid@cells.dim*ext))
}

# expand and wrap a grid on a torus + calc distances
# input: SpatialGrid
# output: distance matrix

ceWrapOnTorusCalcDist <- function(grid, ext=2) {
  grid <- ceExtGrid(grid, ext)
  
  rangeX <- diff(grid@bbox[1,])
  rangeY <- diff(grid@bbox[2,])
  MN.ext <- prod(grid@grid@cells.dim)
  gridCoords <- coordinates(grid)
  
  mmat.ext <- matrix(rep(gridCoords[, 1], MN.ext), MN.ext, MN.ext)
  nmat.ext <- matrix(rep(gridCoords[, 2], MN.ext), MN.ext, MN.ext)
  
  mmat.diff <- mmat.ext - t(mmat.ext)
  nmat.diff <- nmat.ext - t(nmat.ext)
  
  mmat.torus <- pmin(abs(mmat.diff), rangeX - abs(mmat.diff))
  nmat.torus <- pmin(abs(nmat.diff), rangeY - abs(nmat.diff))
  
  sqrt(mmat.torus^2 + nmat.torus^2)
}

## FFT with only first row of cov-matrix
vgmModel <- vgm(25, "Exp", 1)

ceWrapOnTorusCalcDistRow1 <- function(grid, vgmModel, ext=2) {
  grid <- ceExtGrid(grid, ext)
  
  stopifnot("variogramModel" %in% class(vgmModel))
  
  rangeX <- diff(grid@bbox[1,])
  rangeY <- diff(grid@bbox[2,])

  cenX <- seq(from = grid@grid@cellcentre.offset[1],
              by = grid@grid@cellsize[1],
              length.out = grid@grid@cells.dim[1])
  cenY <- seq(from = grid@grid@cellcentre.offset[2],
              by = grid@grid@cellsize[2],
              length.out = grid@grid@cells.dim[2])
  
  m.diff.row1 <- abs(cenX[1] - cenX)
  m.diff.row1 <- pmin(m.diff.row1, rangeX - m.diff.row1)
  
  n.diff.row1 <- abs(cenY[1] - cenY)
  n.diff.row1 <- pmin(n.diff.row1, rangeY - n.diff.row1)

  cent.ext.row1 <- expand.grid(m.diff.row1, n.diff.row1)
  D.ext.row1 <- matrix(sqrt(cent.ext.row1[, 1]^2 + cent.ext.row1[, 2]^2), 
                       grid@grid@cells.dim[1], 
                       grid@grid@cells.dim[2])

  variogramLine(vgmModel, dist_vector = D.ext.row1, covariance = T)
}

covMatRow1 <- ceWrapOnTorusCalcDistRow1(spgrid, vgmModel)
dim(covMatRow1) 

ext.eigs.row1 <- rev(sort(Re(fft(covMatRow1, inverse = TRUE))))
head(ext.eigs.row1)
 
# sim

ceSim <- function(covMatRow1, cells.dim=c(29, 29)) {
t1 <- Sys.time()
d <- dim(covMatRow1)
dp <- prod(d)
sdp <- sqrt(dp)
prefix <- sqrt(Re(fft(covMatRow1, TRUE)))
t2 <- Sys.time()
std <- rnorm(dp)
realz <- prefix * (fft(matrix(std, d[1], d[2]))/sdp)
realz <- as.vector(Re(fft(realz, TRUE)/sdp)[1:cells.dim[1], 1:cells.dim[2]])
realz[!inside.owin(x = cent[, 1], y = cent[, 2], w = W)] <- NA
realization.fft.29 <- matrix(realz, M, N, byrow = TRUE)
t3 <- Sys.time()
timings(t1, t2, t3)

}

library(lattice)
levelplot(t(realization.fft.29), col.regions=heat.colors)

krige(obs~1, )


m.abs.diff.row1 <- abs(mygrid$mcens.ext[1] - mygrid$mcens.ext)
m.diff.row1 <- pmin(m.abs.diff.row1, Rx - m.abs.diff.row1)
n.abs.diff.row1 <- abs(mygrid$ncens.ext[1] - mygrid$ncens.ext)
n.diff.row1 <- pmin(n.abs.diff.row1, Ry - n.abs.diff.row1)
cent.ext.row1 <- expand.grid(m.diff.row1, n.diff.row1)
D.ext.row1 <- matrix(sqrt(cent.ext.row1[, 1]^2 + cent.ext.row1[, 2]^2), mygrid$M.ext, 
                     mygrid$N.ext)
SIGMA.Y.ext.row1 <- sigma^2 * r(D.ext.row1, phi)
dim(SIGMA.Y.ext.row1)

ext.eigs.row1 <- rev(sort(Re(fft(SIGMA.Y.ext.row1, inverse = TRUE))))
head(ext.eigs.row1)



# Rx <- mygrid$M.ext * mygrid$cell.width
# Ry <- mygrid$N.ext * mygrid$cell.height
# MN.ext <- mygrid$M.ext * mygrid$N.ext
# cent.ext <- expand.grid(mygrid$mcens.ext, mygrid$ncens.ext)
# mmat.ext <- matrix(rep(cent.ext[, 1], MN.ext), MN.ext, MN.ext)
# nmat.ext <- matrix(rep(cent.ext[, 2], MN.ext), MN.ext, MN.ext)
# mmat.diff <- mmat.ext - t(mmat.ext)
# nmat.diff <- nmat.ext - t(nmat.ext)
# mmat.torus <- pmin(abs(mmat.diff), Rx - abs(mmat.diff))
# nmat.torus <- pmin(abs(nmat.diff), Ry - abs(nmat.diff))
# D.ext <- sqrt(mmat.torus^2 + nmat.torus^2)




