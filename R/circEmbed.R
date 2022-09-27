## Benedikt Gräler (52North), 2018-05-25:
## circulant embedding following: Davies, Tilman M., and David Bryant. "On
## circulant embedding for Gaussian random fields in R." Journal of Statistical
## Software 55.9 (2013): 1-21.
## See i.e. the suplementary files at (retreived 2018-05-25): 
## https://www.jstatsoft.org/index.php/jss/article/downloadSuppFile/v055i09/v55i09.R


# extend the grid 
# input: SpatialGrid/SpatialPixels/GridTopology
# output: extended grid of class GridTopology
ceExtGrid <- function(grid, ext=2) {
  if (!inherits(grid, "GridTopology")) {
    stopifnot(gridded(grid))
    
    grid <- grid@grid
  }
  
  GridTopology(grid@cellcentre.offset,
               grid@cellsize,
               grid@cells.dim*ext)
}

# only for comaprission following the above paper
# expand and wrap a grid on a torus + calc distances
# input: SpatialGrid
# output: distance matrix
ceWrapOnTorusCalcDist <- function(grid, ext=2) {
  grid <- ceExtGrid(grid, ext)
  
  rangeXY <- grid@cellsize * grid@cells.dim
  
  MN.ext <- prod(grid@cells.dim)
  gridCoords <- coordinates(grid)
  
  mmat.ext <- matrix(rep(gridCoords[, 1], MN.ext), MN.ext, MN.ext)
  nmat.ext <- matrix(rep(gridCoords[, 2], MN.ext), MN.ext, MN.ext)
  
  mmat.diff <- mmat.ext - t(mmat.ext)
  nmat.diff <- nmat.ext - t(nmat.ext)
  
  mmat.torus <- pmin(abs(mmat.diff), rangeXY[1] - abs(mmat.diff))
  nmat.torus <- pmin(abs(nmat.diff), rangeXY[2] - abs(nmat.diff))
  
  sqrt(mmat.torus^2 + nmat.torus^2)
}

## FFT preparation with only first row of cov-matrix
ceWrapOnTorusCalcCovRow1 <- function(grid, vgmModel, ext=2) {
  grid <- ceExtGrid(grid, ext)
  
  stopifnot("variogramModel" %in% class(vgmModel))
  
  rangeXY <- grid@cellsize * grid@cells.dim
  
  cenX <- seq(from = grid@cellcentre.offset[1],
              by = grid@cellsize[1],
              length.out = grid@cells.dim[1])
  cenY <- seq(from = grid@cellcentre.offset[2],
              by = grid@cellsize[2],
              length.out = grid@cells.dim[2])
  
  m.diff.row1 <- abs(cenX[1] - cenX)
  m.diff.row1 <- pmin(m.diff.row1, rangeXY[1] - m.diff.row1)
  
  n.diff.row1 <- abs(cenY[1] - cenY)
  n.diff.row1 <- pmin(n.diff.row1, rangeXY[2] - n.diff.row1)
  
  cent.ext.row1 <- expand.grid(m.diff.row1, n.diff.row1)
  D.ext.row1 <- matrix(sqrt(cent.ext.row1[, 1]^2 + cent.ext.row1[, 2]^2), 
                       grid@cells.dim[1], 
                       grid@cells.dim[2])
  
  variogramLine(vgmModel, dist_vector = D.ext.row1, covariance = T)
}

# simulate GRF with given covariance structure using fft
# @input
#   covMatRow1: the first row of the covariance matrix for the fft
#   n: number of simulations
#   cells.dim: the original dimrensions of the grid to clip from the larger embedded simulation
#   grid.index: grid.index of a SpatialPixels object to select the right pixels from the larger square-grid
# @output
#   matrix where each column holds one simulated GRF corresponding to cells.dim (and grid.index if appropriate)
ceSim <- function(covMatRow1, n=1, cells.dim, grid.index) {
  d <- dim(covMatRow1)
  dp <- prod(d)
  sdp <- sqrt(dp)
  prefix <- sqrt(Re(fft(covMatRow1, TRUE)))
  
  simFun <- function(x) {
    std <- rnorm(dp)
    realz <- prefix * (fft(matrix(std, d[1], d[2]))/sdp)
    as.numeric(Re(fft(realz, TRUE)/sdp)[1:cells.dim[1], 1:cells.dim[2]])
  }
  
  simFunGridIndex <- function(x) {
    std <- rnorm(dp)
    realz <- prefix * (fft(matrix(std, d[1], d[2]))/sdp)
    as.numeric(Re(fft(realz, TRUE)/sdp)[1:cells.dim[1], 1:cells.dim[2]])[grid.index]
  }
  
  if (missing(grid.index))
    do.call(cbind, lapply(1:n, simFun))
  else
    do.call(cbind, lapply(1:n, simFunGridIndex))
}

# computes the covariance matrixes and weights once, applied to series of
# variables/simulations where each variable/simulation is stored in one column of
# the multiVarMatrix copied from krige0 to avoid repeted calls to krige with
# multiple, identical inversions of the weights matrix

krigeMultiple <- function(formula, from, to, model, multiVarMatrix) {
  lst = extractFormula(formula, from, to)
  X = lst$X
  x0 = lst$x0
  ll = (!is.na(is.projected(from)) && !is.projected(from))
  s = coordinates(from)
  s0 = coordinates(to)
  
  V = variogramLine(model,
                    dist_vector = spDists(s, s, ll),
                    covariance = TRUE)
  v0 = variogramLine(model,
                     dist_vector = spDists(s, s0, ll),
                     covariance = TRUE)
  
  skwts = CHsolve(V, cbind(v0, X))
  ViX = skwts[, -(1:nrow(s0))]
  skwts = skwts[, 1:nrow(s0)]
  
  idPredFun <- function(sim) {
    sim <- matrix(sim, ncol = 1)
    beta = solve(t(X) %*% ViX, t(ViX) %*% sim)
    x0 %*% beta + t(skwts) %*% (sim - X %*% beta)
  }
  
  apply(multiVarMatrix, 2, idPredFun)
}

### pubic function ###
# @input: 
# formula: definition of the dependent variable
# data:    optional Spatial*DataFrame for conditional simulation
# newdata: SpatialGrid or SpatialPixels
# model:   variogram model of the GRF
# n:       number of desired simulations
# ext:     extension degree of the circulant embedding, default to 2
# @output
# SpatialPixels or SpatailGridDataFrame with (additional) n columns holding one (un)conditional simulation each

krigeSimCE <- function(formula, data, newdata, model, n = 1, ext = 2) {
  stopifnot(is(model, "variogramModel"))
  stopifnot(gridded(newdata))
  if (!missing(data))
    stopifnot(identical(data@proj4string@projargs, newdata@proj4string@projargs))
  
  varName <- all.vars(formula[[2]])
  
  condSim <- TRUE
  if (missing(data)) {
    condSim <- FALSE
    message("[No data provided: performing unconditional simulation.]")
  } else {
    message("[Performing conditional simulation.]")
  }
  
  # prepare covariance matrix
  covMat <- ceWrapOnTorusCalcCovRow1(newdata, model, ext = ext)
  
  # simulate
  sims <- ceSim(covMat, n, newdata@grid@cells.dim, newdata@grid.index)
  colnames(sims) <- paste0(varName, ".sim", 1:n)
  
  # bind simulations to newdata geometry
  if (!condSim) {
    if ("data" %in% slotNames(newdata))
      newdata@data <- cbind(newdata@data, sims)
    else
      addAttrToGeom(newdata, as.data.frame(sims))
    return(newdata)
  }
  
  # function call ends here if no data has been provided -> unconditional case
  
  ## conditioning
  # interpolate the observations to the simulation grid
  obsMeanField <- krige(formula, data, newdata, model)
  
  # interpolate to observation locations from the simulated grids for each simulation
  simMeanObsLoc <- krigeMultiple(as.formula(paste0("var1.pred ~", formula[[3]])),
                                 obsMeanField, data, model, sims)
  
  # interpolate from kriged mean sim at observed locations back to the grid for mean surface of the simulations
  simMeanFields <- krigeMultiple(as.formula(paste0(varName, "~", formula[[3]])),
                                 data, newdata, model, simMeanObsLoc)
  
  # add up the mean field and the corrected data
  sims <- obsMeanField@data$var1.pred + sims - simMeanFields
  
  # bind simulations to newdata geometry
  if ("data" %in% slotNames(newdata)) {
    newdata@data <- cbind(newdata@data, sims)
    return(newdata)
  }
  
  addAttrToGeom(newdata, as.data.frame(sims))
}

###
# Note: to avoid the smoothing effect, the irregular observation locations could also be independently simulated by e.g. their surrounding grid values
###

## circulant embedding for ST-ocvariance functions with grids along one spatial and one temporal axis


# inputs: 
#   hDiscrete = c(hStep, hn): spatial step width and number of steps
#   tDiscrete = c(tStep, tn): temporal step length and number of steps

# CAVE: hDiscrete and tDiscrete must have the correct spatial and temporal metrics

ceWrapSpaceTimeOnTorusCalcCovRow1 <- function(hDiscrete, tDiscrete, vgmStModel, ext=2, turningLayers=TRUE) {
  stopifnot(is(vgmStModel)  == "StVariogramModel")
  
  hDiscrete[2] <- hDiscrete[2]*ext
  tDiscrete[2] <- tDiscrete[2]*ext
  
  rangeST <- c(prod(hDiscrete), prod(tDiscrete))
  
  cenX <- seq(from = 0, by = hDiscrete[1], length.out = hDiscrete[2])
  cenY <- seq(from = 0, by = tDiscrete[1], length.out = tDiscrete[2])
  
  m.diff.row1 <- abs(cenX[1] - cenX)
  m.diff.row1 <- pmin(m.diff.row1, rangeST[1] - m.diff.row1)
  
  n.diff.row1 <- abs(cenY[1] - cenY)
  n.diff.row1 <- pmin(n.diff.row1, rangeST[2] - n.diff.row1)
  
  cent.ext.row1 <- expand.grid(m.diff.row1, n.diff.row1)
  D.ext.row1 <- matrix(sqrt(cent.ext.row1[, 1]^2 + cent.ext.row1[, 2]^2), 
                       hDiscrete[2], 
                       tDiscrete[2])
  colnames(cent.ext.row1) <- c("spacelag", "timelag")
  
  if (turningLayers) {
    return(matrix(tbOperator(vgmStModel, dist_grid = cent.ext.row1)$gamma,
                  hDiscrete[2], tDiscrete[2]))
  } 
  
  matrix(variogramSurface(vgmStModel, dist_grid = cent.ext.row1, 
                          covariance = TRUE)$gamma,
         hDiscrete[2], tDiscrete[2])
}


     
