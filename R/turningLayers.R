# 
# library(gstat)
# # separable model: spatial and temporal sill will be ignored
# # and kept constant at 1-nugget respectively. A joint sill is used.
# separableModel <- vgmST("separable", 
#                         space=vgm(0.9,"Exp", 147, 0.1),
#                         time =vgm(0.9,"Exp", 3.5, 0.1),
#                         sill=40)
# 
# covSTMat <- ceWrapSpaceTimeOnTorusCalcCovRow1(c(250,20), c(1,16), separableModel)
# dim(covSTMat)
# 
# image(matrix(ceSim(covSTMat, 1, c(20, 16)),
#              nrow=20, ncol=16), asp=1, ylab="time", xlab="space")
# 
# # product sum model: spatial and temporal nugget will be ignored and kept
# # constant at 0. Only a joint nugget is used.
# prodSumModel <- vgmST("productSum",
#                       space=vgm(39, "Sph", 343, 0),
#                       time= vgm(36, "Exp",   3, 0), 
#                       k=15)
# 
# covSTMat <- ceWrapSpaceTimeOnTorusCalcCovRow1(c(250,20), c(1,16), prodSumModel)
# dim(covSTMat)
# 
# image(matrix(ceSim(covSTMat, 1, c(20, 16)),
#              nrow=20, ncol=16), asp=1, ylab="time", xlab="space")
# 
# # sum metric model: spatial, temporal and joint nugget will be estimated
# sumMetricModel <- vgmST("sumMetric",
#                         space=vgm( 6.9, "Lin", 200, 3.0),
#                         time =vgm(10.3, "Lin",  15, 3.6),
#                         joint=vgm(37.2, "Exp",  84,11.7),
#                         stAni=77.7)
# 
# covSTMat <- ceWrapSpaceTimeOnTorusCalcCovRow1(c(250,20), c(1,16), sumMetricModel)
# dim(covSTMat)
# 
# image(matrix(ceSim(covSTMat, 1, c(20, 16)),
#              nrow=20, ncol=16), asp=1, ylab="time", xlab="space")
# 
# # simplified sumMetric model, only a overall nugget is fitted. The spatial, 
# # temporal and jont nuggets are set to 0.
# simpleSumMetricModel <- vgmST("simpleSumMetric",
#                               space=vgm(20,"Lin", 150, 0),
#                               time =vgm(20,"Lin", 10,  0),
#                               joint=vgm(20,"Exp", 150, 0),
#                               nugget=1, stAni=15)
# 
# covSTMat <- ceWrapSpaceTimeOnTorusCalcCovRow1(c(250,20), c(1,16), simpleSumMetricModel)
# dim(covSTMat)
# 
# image(matrix(ceSim(covSTMat, 1, c(20, 16)),
#              nrow=20, ncol=16), asp=1, ylab="time", xlab="space")
# 
# # metric model
# metricModel <- vgmST("metric",
#                      joint=vgm(60, "Exp", 150, 10),
#                      stAni=60)
# 
# covSTMat <- ceWrapSpaceTimeOnTorusCalcCovRow1(c(250,20), c(1,16), metricModel, turningLayers = FALSE)
# dim(covSTMat)
# 
# image(matrix(ceSim(covSTMat, 1, c(20, 16)),
#              nrow=20, ncol=16), asp=1, ylab="time", xlab="space")

## turning bands:
# for cov-functions valid in 3-dim 

# adjust the covariance matrix (Schlahter, Chap. 2, eq 2.25)
# simulate several layers (up-right hyper-planes in the 3D+T-cube) with random directions on a unit hemi-sphere
# orthogonally project target points onto rndm lyrs and average the values for all rndm lyrs

# -> one simulation of a 3d+T Gaussian random field, for 2D use a hyperplane in 3D

# the tunring bands operator (numerical, 3-dimensional case)
# @input:
#  model target covariance function 
#  dist_grid: data.frame wtih columns spacelag and timelag for which the covariances are calculated
#
# @value: data.frame with columns spacelag, timelag and gamma where the latter contains the covariances for the 1-dim spatial and temporal layers

# model <- separableModel
# dist_grid <- as.data.frame(cbind("spacelag" = rep(1:150*1., each=4),
#                                  "timelag" = rep(1:4, 150)))

# CAVE: distgrid must ahve the correct spatial and temporal metrics

tbOperator <- function(model, dist_grid) {
  r <- dist_grid$spacelag
  
  derFun <- function(r) r * variogramSurface(model, dist_grid = data.frame(spacelag=r, timelag=dist_grid$timelag), 
                                             covariance = TRUE)$gamma
  
  cbind(dist_grid, "gamma" = diag(attr(numericDeriv(quote(derFun(r)), "r"), "gradient")))
}

randomDirections <- function(n) {
  u <- runif(n, 0, 2*pi)
  v <- runif(n, 0, 1)
  s <- sqrt(1-v^2)
  
  cbind(s*cos(u), s*sin(u), v)
}

# library(rgl)
# plot3d(randomDirections(100), aspect = c(1,1,0.5))


# nLyrs <- 500
# coordMat <- matrix(c(1,5,5,
#                      1,7,7,
#                      1,5,7,
#                      1,7,5), byrow=T, ncol=3)
# 
# coordMat <- matrix(runif(3*1001,1,5), byrow=T, ncol=3)
# 
# # random directions in 3D
# rndDir <- randomDirections(nLyrs)
# 
# # how much does each direction contribute to the point
# cntrbtn <- rndDir %*% t(coordMat)
# # -> each column corresponds to one location, each row i to the "contribution" of the i-th rndDir
# 
# cl_cntrbtn <- ceiling(cntrbtn)+10 # centre for grid; mind spatial index
# fl_cntrbtn <- floor(cntrbtn) +10
# 
# covRow1 <- ceWrapSpaceTimeOnTorusCalcCovRow1(c(250,20), c(1,16), metricModel, turningLayers = TRUE)
# 
# sTime <- Sys.time()
# simLyrs <- ceSim(covRow1, nLyrs, c(20,16))
# simLyrs <- lapply(1:nLyrs, function(col) matrix(simLyrs[,col], nrow=20, ncol=16))
# 
# lambda <- cl_cntrbtn - cntrbtn -10
# 
# cntrbSngSimLyr <- function(lyrId) {
#   simLyrs[[lyrId]][cl_cntrbtn[lyrId,],] * (1-lambda)[lyrId,] + simLyrs[[lyrId]][fl_cntrbtn[lyrId,],] * lambda[lyrId,]
# }
# 
# # reduce to the one realisation based on the combination of nLyrs turning bands
# simTs <- Reduce('+', lapply(1:nLyrs, cntrbSngSimLyr))/sqrt(nLyrs)
# eTime <- Sys.time()
# 
# eTime - sTime # 0.7 secs/simulation of 100 random points
# 
# image(simTs)
# 
# dim(simTs)
# # 1001 locations at 16 time stamps
# 
# plot(simTs[500,])

# computes the covariance matrixes and weights once, applied to series of
# variables/simulations where each variable/simulation is stored in one column of
# the multiVarMatrix copied from krigeST to avoid repeted calls to krige with
# multiple, identical inversions of the weights matrix

# TODO: add functionality for temporal sandwich-wise processing: i.e. use the 
# +/- nmaxTime time slices to predict one time slice

krigeSTMultiple <- function(formula, from, to, modelList, multiVarMatrix, nmaxTime=Inf) {
  lst = extractFormula(formula, from, to)

  separate <- length(from) > 1 && length(to) > 1 &&
    inherits(from, "STF") && inherits(to, "STF")
  
  X = lst$X
  x0 = lst$x0
  
  V = covfn.ST(from, model = modelList, separate=separate)
  v0 = covfn.ST(from, to, modelList)
  
  if (modelList$stModel == "separable" & separate)
    skwts <- STsolve(V, v0, X) # use Kronecker trick
  else 
    skwts <- CHsolve(V, cbind(v0, X))

  npts = length(to)
  ViX = skwts[,-(1:npts)]
  skwts = skwts[,1:npts]
  
  idPredFun <- function(sim) {
    sim <- matrix(sim, ncol = 1)
    beta = solve(t(X) %*% ViX, t(ViX) %*% sim)
    x0 %*% beta + t(skwts) %*% (sim - X %*% beta)
  }
  
  apply(multiVarMatrix, 2, idPredFun)
}



# unconcitional
## STFDF
# library(sp)
# library(spacetime)
# library(zoo)
# library(xts)
# data(meuse, package = "sp")
# coordinates(meuse) <- ~x+y
# proj4string(meuse) <- CRS("+init=epsg:28992")

# krigeST <- function(formula, data, newdata, modelList, y, beta, nmax=Inf, stAni=NULL,
#                     computeVar = FALSE, fullCovariance = FALSE,
#                     bufferNmax=2, progress=TRUE)

## 
krigeSTSimTB <- function(formula, data, newdata, modelList, nsim, 
                            progress=TRUE, nLyrs=500, tGrid=NULL, sGrid=NULL, ceExt=2,
                            nmax=Inf) {
  stopifnot(zoo::is.regular(newdata@time))
  
  condSim <- TRUE
  if (missing(data)) {
    condSim <- FALSE
    message("[No data provided: performing unconditional simulation.]")
  } else {
    message("[Performing conditional simulation.]")
  }
  
  pb <- txtProgressBar(0,nsim,style=3)
  
  # ST-simulation grid
  if (is.null(tGrid)) {
    tDis <- diff(c(index(newdata@time[1]), newdata@endTime[1]))
    if (!is.null(attr(modelList, "temporal unit"))) {
      units(tDis) <- attr(modelList, "temporal unit")
    } else {
      message("[The spatio-temporal variogram model does not carry a time unit attribute: krigeST cannot check whether the temporal distance metrics coincide.]")
    }
    
    tGrid <- c(as.numeric(tDis), length(newdata@time))
    attr(tGrid, "units") <- c("", units(tDis))
    
    debug_time_unit(units(tDis))
  }
  
  if (is.null(sGrid)) {
    if (gridded(newdata@sp)) {
      # SpatialPixels/SpatialGrid:
      # based on GridTopology: use minimal cellsize of both directions; take enough to cover the diagonal
      sDis <- min(newdata@sp@grid@cellsize)
      sDim <- ceiling(sqrt(sum((newdata@sp@grid@cellsize * newdata@sp@grid@cells.dim)^2))/sDis)
      sGrid <-  c(sDis, sDim)
    } else {
      # treat (as) SpatialPoints:
      # Average area per location --assuming-a-regular-squared-outline-taking-the-sqrt--> length/location --take-later-twice-as-many-->
      bboxExt <- apply(newdata@sp@bbox, 1, diff)
      sDis <- sqrt(prod(bboxExt)/length(newdata@sp))
      sDim <- ceiling(sqrt(sum((bboxExt/sDis)^2)))
      sGrid <- c(sDis, sDim)
    }
  }
    
  # random directions in 3D
  rndDir <- randomDirections(nLyrs)
  
  # coordinates 
  coordMat <- coordinates(newdata@sp)

  # coordinates (embedded in 3D) shifted + scaled to the grid index of the spatial gridding sGrid
  coordMat[,1] <- (coordMat[,1] - newdata@sp@bbox[1,1])/sGrid[1]
  coordMat[,2] <- (coordMat[,2] - newdata@sp@bbox[2,1])/sGrid[1]
  
  if (ncol(coordMat) == 2) {
    coordMat <- cbind(coordMat,1)
  } else {
    coordMat[,3] <- (coordMat[,3] - newdata@sp@bbox[3,1])/sGrid[1]
  }

  # how much does each direction contribute to the point
  cntrbtn <- rndDir %*% t(coordMat)
  # -> each column corresponds to one location, each row i to the "contribution" of the i-th rndDir
  
  # ceiling and floor + shifted to avoid negative indices
  cl_cntrbtn <- ceiling(cntrbtn) + sGrid[2]
  fl_cntrbtn <- floor(cntrbtn) + sGrid[2]
  
  covRow1 <- ceWrapSpaceTimeOnTorusCalcCovRow1(c(sGrid[1], 2*sGrid[2]), tGrid, modelList, turningLayers = TRUE, ext=ceExt)
  
  origDim <- c(2*sGrid[2], tGrid[2])
  sims <- list()

  ##
  for (i in 1:nsim) {
    setTxtProgressBar(pb, i)
    simLyrs <- ceSim(covRow1, nLyrs, origDim)
    simLyrs <- lapply(1:nLyrs, function(col) matrix(simLyrs[,col], 
                                                    nrow=origDim[1], 
                                                    ncol=origDim[2]))
    lambda <- cl_cntrbtn - cntrbtn - sGrid[2]
    
    cntrbSngSimLyr <- function(lyrId) {
      simLyrs[[lyrId]][cl_cntrbtn[lyrId,],] * (1-lambda)[lyrId,] + simLyrs[[lyrId]][fl_cntrbtn[lyrId,],] * lambda[lyrId,]
    }
    
    # reduce to the one realisation based on the combination of nLyrs turning bands
    sims[[paste0("sim",i)]] <- Reduce('+', lapply(1:nLyrs, cntrbSngSimLyr))/sqrt(nLyrs)
  }
  close(pb)
  
  sims <- do.call(cbind, lapply(sims, as.numeric))

  # bind simulations to newdata geometry
  if (!condSim) {
    if ("data" %in% slotNames(newdata))
      newdata@data <- cbind(newdata@data, sims)
    else
      newdata <- addAttrToGeom(newdata, as.data.frame(sims))
    return(newdata)
  }
  
  # function call ends here if no data has been provided -> unconditional case
  varName <- all.vars(formula[[2]])
  
  ## conditioning
  # interpolate the observations to the simulation grid
  obsMeanField <- krigeST(formula=formula, data=data, newdata=newdata, modelList=modelList)
  
  # interpolate to observation locations from the simulated grids for each simulation
  simMeanObsLoc <- krigeSTMultiple(as.formula(paste0("var1.pred ~", formula[[3]])),
                                   obsMeanField, data, modelList, sims)
  
  # interpolate from kriged mean sim at observed locations back to the grid for mean surface of the simulations
  simMeanFields <- krigeSTMultiple(as.formula(paste0(varName, "~", formula[[3]])),
                                 data, newdata, modelList, simMeanObsLoc)
  
  # add up the mean field and the corrected data
  sims <- obsMeanField@data$var1.pred + sims - simMeanFields
  
  # bind simulations to newdata geometry
  if ("data" %in% slotNames(newdata)) {
    newdata@data <- cbind(newdata@data, sims)
    return(newdata)
  }
  
  addAttrToGeom(newdata, as.data.frame(sims))
}


# 
# sTime <- Sys.time()
# krigedSim <- krigeSTUncSimTB(stf, metricModel, 100)
# Sys.time() - sTime
# 
# # 27 secs for 100 simulated ST fields of 155 locations and 21 time steps: 325500 values
# 
# # plot one simulation along time
# stplot(krigedSim[,1:12])
# 
# # plot one simulation along time
# stplot(krigedSim[1:12,,"sim1"], mode="ts")
# 
# # plot the ten simulations of the first day
# spplot(krigedSim[,1], paste0("sim",1:10), as.table=TRUE)
