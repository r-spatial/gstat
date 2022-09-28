# constructiong spatio-temporal variogram models
vgmST <- function(stModel, ..., space, time, joint, sill, k, nugget, stAni, 
                  temporalUnit) {
  stopifnot(is.character(stModel) && length(stModel)==1)
  
  old.stModel <- stModel
  stModel <- strsplit(stModel, "_")[[1]][1]
  
  if (stModel == "productSum" && !missing(sill))
    stop("The sill argument for the product-sum model has been removed 
due a change in notation of the spatio-temporal models. This 
affects as well how the spatial and temporal variograms are parameterised. 
Re-fit your model or use \"productSumOld\" instead.")
  
  if(!missing(sill))
    if(sill <= 0) stop("\"sill\" must be positive.")
  if(!missing(k))
    if(k <= 0) stop("\"k\" must be positive.")
  if(!missing(nugget))
    if(nugget < 0) stop("\"nugget\" must be non-negative.")
  if(!missing(stAni))
    if(stAni <= 0) stop("\"stAni\" must be positive.")
  
  vgmModel <- switch(stModel,
                     separable = list(space = space, time = time, sill = sill),
                     productSum = list(space = space, time = time, k = k),
                     productSumOld = list(space = space, time = time,
                                          sill = sill, nugget = nugget),
                     sumMetric = list(space = space, time = time, 
                                      joint = joint, stAni = stAni),
                     simpleSumMetric = list(space = space, time = time, 
                                            joint = joint, nugget = nugget, 
                                            stAni = stAni),
                     metric = list(joint = joint, stAni = stAni),
                     stop(paste("model", stModel, "unknown")))
  
  vgmModel$stModel <- old.stModel
  
  if (!missing(temporalUnit))
    attr(vgmModel, "temporal unit") = temporalUnit
  
  class(vgmModel) <- c("StVariogramModel", "list")
  
  vgmModel
}

# calculating spatio-temporal variogram surfaces
variogramSurface <- function(model, dist_grid, covariance=FALSE) {
  stopifnot(inherits(model, "StVariogramModel"))
  stopifnot(all(c("spacelag", "timelag") %in% colnames(dist_grid)))
  
  if (covariance) {
    switch(strsplit(model$stModel, "_")[[1]][1],
           separable=covSurfSeparable(model, dist_grid),
           productSum=covSurfProdSum(model, dist_grid),
           productSumOld=covSurfProdSumOld(model, dist_grid),
           sumMetric=covSurfSumMetric(model, dist_grid),
           simpleSumMetric=covSurfSimpleSumMetric(model, dist_grid),
           metric=covSurfMetric(model, dist_grid),
           stop("Only \"separable\", \"productSum\", \"sumMetric\", \"simpleSumMetric\" and \"metric\" are implemented."))
  } else {
    switch(strsplit(model$stModel, "_")[[1]][1],
           separable=vgmSeparable(model, dist_grid),
           productSum=vgmProdSum(model, dist_grid),
           productSumOld=vgmProdSumOld(model, dist_grid),
           sumMetric=vgmSumMetric(model, dist_grid),
           simpleSumMetric=vgmSimpleSumMetric(model, dist_grid),
           metric=vgmMetric(model, dist_grid),
           stop("Only \"separable\", \"productSum\", \"sumMetric\", \"simpleSumMetric\" and \"metric\" are implemented."))
  }
}

################################
## separable model: C_s * C_t ##
################################

vgmSeparable <- function(model, dist_grid) {
  vs = variogramLine(model$space, dist_vector=dist_grid$spacelag)$gamma
  vt = variogramLine(model$time,  dist_vector=dist_grid$timelag)$gamma
  
  cbind(dist_grid, "gamma" = model$sill*(vs+vt-vs*vt))
}

covSeparable <- function(x, y, model, separate) {  
  if(missing(separate))
    separate <- inherits(x, "STF") && inherits(y, "STF") && length(x) > 1 && length(y) > 1
  
  # the STF case
  if (inherits(x, "STF") && inherits(y, "STF")) {
    # calculate all spatial and temporal distances
    ds = spDists(x@sp, y@sp)
    dt = abs(outer(index(x@time), index(y@time), "-"))
    if(!is.null(attr(model,"temporal unit")))
      units(dt) <- attr(model, "temporal unit") # ensure the same temporal metric as in the variogram definition
    debug_time_unit(units(dt))
    dt <- as(dt, "matrix")
    
    # compose the cov-matrix
    Sm = variogramLine(model$space, covariance = TRUE, dist_vector = ds)*model$sill
    Tm = variogramLine(model$time, covariance = TRUE, dist_vector = dt)
    
    if (separate)
      return(list(Sm = Sm, Tm = Tm))
    else
      return(Tm %x% Sm) # kronecker product
  } 
  
  # separate makes only sense if both of x and y inherit STF
  if (separate)
    stop("An efficient inversion by separating the covarinace model is only possible if both of \"x\" and \"y\" inherit \"STF\"")
  
  # the STI case
  if (inherits(x, c("STI", "sftime")) || inherits(y, c("STI", "sftime"))) {
    # make sure that now both are of type STI
    x <- as(x, "STI")
    y <- as(y, "STI")
    
    # calculate all spatial and temporal distances
    ds = spDists(x@sp, y@sp)
    dt = abs(outer(index(x@time), index(y@time), "-"))
    if(!is.null(attr(model,"temporal unit")))
      units(dt) <- attr(model, "temporal unit") # ensure the same temporal metric as in the variogram definition
    debug_time_unit(units(dt))
    dt <- as(dt, "matrix")
    
    # compose the cov-matrix
    Sm = variogramLine(model$space, covariance = TRUE, dist_vector = ds)*model$sill
    Tm = variogramLine(model$time, covariance = TRUE, dist_vector = dt)
    
    return(Sm * Tm)
  }
  
  # the remaining cases, none of x and y is STI nor are both STF
  # make sure both are of type STS
  x <- as(x, "STS")
  y <- as(y, "STS")
  
  # calculate all spatial and temporal distances
  ds = spDists(x@sp, y@sp)
  dt = abs(outer(index(x@time), index(y@time), "-"))
  if(!is.null(attr(model,"temporal unit")))
    units(dt) <- attr(model, "temporal unit") # ensure the same temporal metric as in the variogram definition
  debug_time_unit(units(dt))
  dt <- as(dt, "matrix")
  
  # re-arrange the spatial and temporal distances
  sMat <- matrix(NA, nrow(x@index), nrow(y@index))
  tMat <- matrix(NA, nrow(x@index), nrow(y@index))
  for(r in 1:nrow(x@index)) {
    sMat[r,] <- ds[x@index[r,1], y@index[,1]]
    tMat[r,] <- dt[x@index[r,2], y@index[,2]]
  }
  
  # compose the cov-matrix
  Sm = variogramLine(model$space, covariance = TRUE, dist_vector = sMat)*model$sill
  Tm = variogramLine(model$time, covariance = TRUE, dist_vector = tMat)
  
  return(Sm * Tm)  
}

# covariance for the circulant embedding in ST
covSurfSeparable <- function(model, dist_grid) {
  Sm = variogramLine(model$space, covariance = TRUE, dist_vector = dist_grid$spacelag)$gamma*model$sill
  Tm = variogramLine(model$time, covariance = TRUE, dist_vector = dist_grid$timelag)$gamma
  
  cbind(dist_grid, "gamma" = Tm * Sm)
}

###########################################
## productSum model: C_s*C_t + C_s + C_t ##
###########################################

vgmProdSumOld <- function(model, dist_grid) {
  .Deprecated("vgmProdSum", package = "gstat", 
              msg="The former product-sum model is dprecited, consider to refit the new model specification",
              old = "vgmProdSumOld")
  
  vs = variogramLine(model$space, dist_vector=dist_grid$spacelag)$gamma
  vt = variogramLine(model$time, dist_vector=dist_grid$timelag)$gamma
  vn <- rep(model$nugget, length(vs))
  vn[vs == 0 & vt == 0] <- 0
  
  k <- (sum(model$space$psill)+sum(model$time$psill)-(model$sill+model$nugget))/(sum(model$space$psill)*sum(model$time$psill))
  
  if (k <= 0 || k > 1/max(rev(model$space$psill)[1], rev(model$time$psill)[1])) 
    k <- 10^6*abs(k) # distorting the model to let optim "hopefully" find suitable parameters
  
  cbind(dist_grid, "gamma" = as.vector(vs+vt-k*vs*vt+vn))
}

covProdSumOld <- function(x, y, model) {
  .Deprecated("covProdSum", package = "gstat", 
              msg="The former product-sum model is dprecited, consider to refit the new model specification",
              old = "covProdSumOld")
  
  stopifnot(inherits(x, c("STF", "STS", "STI", "sftime")) &&
			inherits(y, c("STF", "STS", "STI", "sftime")))
  
  # double check model for validity, i.e. k:
  k <- (sum(model$space$psill)+sum(model$time$psill)-model$sill)/(sum(model$space$psill)*sum(model$time$psill))
  if (k <= 0 || k > 1/max(model$space$psill[model$space$model!="Nug"], 
                         model$time$psill[model$time$model!="Nug"]))
    stop(paste("k (",k,") is non-positive or too large: no valid model!",sep=""))
  
  # the STF case
  if (inherits(x, "STF") && inherits(y, "STF")) {
    # calculate all spatial and temporal distances
    ds = spDists(x@sp, y@sp)
    dt = abs(outer(index(x@time), index(y@time), "-"))
    if(!is.null(attr(model,"temporal unit")))
      units(dt) <- attr(model, "temporal unit") # ensure the same temporal metric as in the variogram definition
    debug_time_unit(units(dt))
    dt <- as(dt, "matrix")
    
    # compose the cov-matrix
    vs = variogramLine(model$space, dist_vector = ds, covariance = TRUE)
    vt = variogramLine(model$time, dist_vector = dt, covariance = TRUE)
    
    return(model$sill-(vt %x% matrix(1,nrow(vs),ncol(vs)) + matrix(1,nrow(vt),ncol(vt)) %x% vs - k * vt %x% vs))
  } 
  
  # the STI case
  if(inherits(x, c("STI", "sftime")) || inherits(y, c("STI", "sftime"))) {
    # make sure that now both are of type STI
    x <- as(x, "STI")
    y <- as(y, "STI")
    
    # calculate all spatial and temporal distances
    ds = spDists(x@sp, y@sp)
    dt = abs(outer(index(x@time), index(y@time), "-"))
    if(!is.null(attr(model,"temporal unit")))
      units(dt) <- attr(model, "temporal unit") # ensure the same temporal metric as in the variogram definition
    debug_time_unit(units(dt))
    dt <- as(dt, "matrix")
    
    # compose the cov-matrix
    vs = variogramLine(model$space, dist_vector = ds, covariance = TRUE)
    vt = variogramLine(model$time, dist_vector = dt, covariance = TRUE)
    
    return(model$sill-(vt + vs - k * vt * vs))
  }
  
  # the remaining cases, none of x and y is STI nor are both STF
  # make sure both are of type STS
  x <- as(x, "STS")
  y <- as(y, "STS")
  
  # calculate all spatial and temporal distances
  ds = spDists(x@sp, y@sp)
  dt = abs(outer(index(x@time), index(y@time), "-"))
  if(!is.null(attr(model,"temporal unit")))
    units(dt) <- attr(model, "temporal unit") # ensure the same temporal metric as in the variogram definition
  debug_time_unit(units(dt))
  dt <- as(dt, "matrix")
  
  # re-arrange the spatial and temporal distances
  sMat <- matrix(NA, nrow(x@index), nrow(y@index))
  tMat <- matrix(NA, nrow(x@index), nrow(y@index))
  for(r in 1:nrow(x@index)) {
    sMat[r,] <- ds[x@index[r,1], y@index[,1]]
    tMat[r,] <- dt[x@index[r,2], y@index[,2]]
  }
  
  # compose the cov-matrix
  vs = variogramLine(model$space, dist_vector = sMat, covariance = TRUE)
  vt = variogramLine(model$time, dist_vector = tMat, covariance = TRUE)
  
  return(model$sill-(vt + vs - k * vt * vs))
}

# covariance for the circulant embedding in ST
covSurfProdSumOld <- function(model, dist_grid) {
  .Deprecated("covSurfProdSum", package = "gstat", 
              msg="The former product-sum model is dprecited, consider to refit the new model specification",
              old = "covSurfProdSumOld")
  
  vs = variogramLine(model$space, dist_vector = dist_grid$spacelag, covariance = TRUE)$gamma
  vt = variogramLine(model$time,  dist_vector = dist_grid$timelag, covariance = TRUE)$gamma
  
  k <- (sum(model$space$psill)+sum(model$time$psill)-(model$sill+model$nugget))/(sum(model$space$psill)*sum(model$time$psill))
  
  cbind(dist_grid, "gamma" = model$sill-(vt + vs - k * vt * vs))
}

vgmProdSum <- function(model, dist_grid) {
  if(!is.null(model$sill)) # backwards compatibility
    vgmProdSumOld(model, dist_grid)
  vs = variogramLine(model$space, dist_vector=dist_grid$spacelag)$gamma
  vt = variogramLine(model$time, dist_vector=dist_grid$timelag)$gamma
  
  sill_s <- sum(model$space$psill)
  sill_t <- sum(model$time$psill)
  k <- model$k
  
  cbind(dist_grid, "gamma" = as.vector((k*sill_t+1)*vs + (k*sill_s+1)*vt-k*vs*vt))
}

covProdSum <- function(x, y, model) {
  stopifnot(inherits(x, c("STF", "STS", "STI", "sftime")) &&
			inherits(y, c("STF", "STS", "STI", "sftime")))
  if(!is.null(model$sill)) # backwards compatibility
    covProdSumOld(x, y, model)
  
  # the STF case
  if (inherits(x, "STF") && inherits(y, "STF")) {
    # calculate all spatial and temporal distances
    ds = spDists(x@sp, y@sp)
    dt = abs(outer(index(x@time), index(y@time), "-"))
    if(!is.null(attr(model,"temporal unit")))
      units(dt) <- attr(model, "temporal unit") # ensure the same temporal metric as in the variogram definition
    debug_time_unit(units(dt))
    dt <- as(dt, "matrix")
    
    # compose the cov-matrix
    vs = variogramLine(model$space, dist_vector = ds, covariance =TRUE)
    vt = variogramLine(model$time, dist_vector = dt, covariance =TRUE)
    
    return(vt %x% matrix(1,nrow(vs),ncol(vs)) + matrix(1,nrow(vt),ncol(vt)) %x% vs + model$k * vt %x% vs)
  } 
  
  # the STI case
  if(inherits(x, c("STI", "sftime")) || inherits(y, c("STI", "sftime"))) {
    # make sure that now both are of type STI
    x <- as(x, "STI")
    y <- as(y, "STI")
    
    # calculate all spatial and temporal distances
    ds = spDists(x@sp, y@sp)
    dt = abs(outer(index(x@time), index(y@time), "-"))
    if(!is.null(attr(model,"temporal unit")))
      units(dt) <- attr(model, "temporal unit") # ensure the same temporal metric as in the variogram definition
    debug_time_unit(units(dt))
    dt <- as(dt, "matrix")
    
    # compose the cov-matrix
    vs = variogramLine(model$space, dist_vector = ds, covariance=TRUE)
    vt = variogramLine(model$time, dist_vector = dt, covariance=TRUE)
    
    return(vt + vs + model$k * vt * vs)
  }
  
  # the remaining cases, none of x and y is STI nor are both STF
  # make sure both are of type STS
  x <- as(x, "STS")
  y <- as(y, "STS")
  
  # calculate all spatial and temporal distances
  ds = spDists(x@sp, y@sp)
  dt = abs(outer(index(x@time), index(y@time), "-"))
  if(!is.null(attr(model,"temporal unit")))
    units(dt) <- attr(model, "temporal unit") # ensure the same temporal metric as in the variogram definition
  debug_time_unit(units(dt))
  dt <- as(dt, "matrix")
  
  # re-arrange the spatial and temporal distances
  sMat <- matrix(NA, nrow(x@index), nrow(y@index))
  tMat <- matrix(NA, nrow(x@index), nrow(y@index))
  for(r in 1:nrow(x@index)) {
    sMat[r,] <- ds[x@index[r,1], y@index[,1]]
    tMat[r,] <- dt[x@index[r,2], y@index[,2]]
  }
  
  # compose the cov-matrix
  vs = variogramLine(model$space, dist_vector = sMat, covariance = TRUE)
  vt = variogramLine(model$time, dist_vector = tMat, covariance = TRUE)
  
  return(vt + vs + model$k * vt * vs)
}

# covariance for the circulant embedding in ST
covSurfProdSum <- function(model, dist_grid) {
  vs = variogramLine(model$space, dist_vector = dist_grid$spacelag, covariance = TRUE)$gamma
  vt = variogramLine(model$time, dist_vector = dist_grid$timelag, covariance = TRUE)$gamma
  
  cbind(dist_grid, "gamma" = vt + vs + model$k * vt * vs)
}

#########################################################
# sumMetric model: C_s + C_t + C_st (Gerard Heuvelink) ##
#########################################################

vgmSumMetric <- function(model, dist_grid) {
  vs = variogramLine(model$space, dist_vector=dist_grid$spacelag)$gamma
  vt = variogramLine(model$time,  dist_vector=dist_grid$timelag)$gamma
  h = sqrt(dist_grid$spacelag^2 + (model$stAni * as.numeric(dist_grid$timelag))^2)
  vst = variogramLine(model$joint, dist_vector=h)$gamma
  
  cbind(dist_grid, "gamma" = vs + vt + vst)
}

covSumMetric <- function(x, y, model) {
  stopifnot(inherits(x, c("STF", "STS", "STI", "sftime")) &&
			inherits(y, c("STF", "STS", "STI", "sftime")))
  
  # the STF case
  if (inherits(x, "STF") && inherits(y, "STF")) {
    # calculate all spatial and temporal distances
    ds = spDists(x@sp, y@sp)
    dt = abs(outer(index(x@time), index(y@time), "-"))
    if(!is.null(attr(model,"temporal unit")))
      units(dt) <- attr(model, "temporal unit") # ensure the same temporal metric as in the variogram definition
    debug_time_unit(units(dt))
    dt <- as(dt, "matrix")
    
    # compose the cov-matrix
    Sm = variogramLine(model$space, covariance = TRUE, dist_vector = ds)
    Tm = variogramLine(model$time, covariance = TRUE, dist_vector = dt)
    
    h  = sqrt((matrix(1,nrow(dt),ncol(dt)) %x% ds)^2 
              + (model$stAni * dt %x% matrix(1,nrow(ds),ncol(ds)))^2)
    Mm = variogramLine(model$joint, covariance = TRUE, dist_vector = h)
    
    return(matrix(1,nrow(Tm),ncol(Tm)) %x% Sm + Tm %x% matrix(1,nrow(Sm),ncol(Sm)) + Mm)
  } 
  
  # the STI case
  if(inherits(x, c("STI", "sftime")) || inherits(y, c("STI", "sftime"))) {
    # make sure that now both are of type STI
    x <- as(x, "STI")
    y <- as(y, "STI")
    
    # calculate all spatial and temporal distances
    ds = spDists(x@sp, y@sp)
    dt = abs(outer(index(x@time), index(y@time), "-"))
    if(!is.null(attr(model,"temporal unit")))
      units(dt) <- attr(model, "temporal unit") # ensure the same temporal metric as in the variogram definition
    debug_time_unit(units(dt))
    dt <- as(dt, "matrix")
    
    # compose the cov-matrix
    Sm = variogramLine(model$space, covariance = TRUE, dist_vector = ds)
    Tm = variogramLine(model$time, covariance = TRUE, dist_vector = dt)
    
    h  = sqrt(ds^2 + (model$stAni * dt)^2)
    Mm = variogramLine(model$joint, covariance = TRUE, dist_vector = h)
    
    return(Sm + Tm + Mm)
  }
  
  # the remaining cases, none of x and y is STI nor are both STF
  # make sure both are of type STS
  x <- as(x, "STS")
  y <- as(y, "STS")
  
  # calculate all spatial and temporal distances
  ds = spDists(x@sp, y@sp)
  dt = abs(outer(index(x@time), index(y@time), "-"))
  if(!is.null(attr(model,"temporal unit")))
    units(dt) <- attr(model, "temporal unit") # ensure the same temporal metric as in the variogram definition
  debug_time_unit(units(dt))
  dt <- as(dt, "matrix")
  
  # re-arrange the spatial and temporal distances
  sMat <- matrix(NA, nrow(x@index), nrow(y@index))
  tMat <- matrix(NA, nrow(x@index), nrow(y@index))
  for(r in 1:nrow(x@index)) {
    sMat[r,] <- ds[x@index[r,1], y@index[,1]]
    tMat[r,] <- dt[x@index[r,2], y@index[,2]]
  }
  
  # compose the cov-matrix
  Sm = variogramLine(model$space, covariance = TRUE, dist_vector = sMat)
  Tm = variogramLine(model$time, covariance = TRUE, dist_vector = tMat)
  
  h  = sqrt(sMat^2 + (model$stAni * tMat)^2)
  Mm = variogramLine(model$joint, covariance = TRUE, dist_vector = h)
  
  return(Sm + Tm + Mm)
}

# covariance for the circulant embedding in ST
covSurfSumMetric <- function(model, dist_grid) {
  Sm = variogramLine(model$space, covariance = TRUE, dist_vector = dist_grid$spacelag)$gamma
  Tm = variogramLine(model$time, covariance = TRUE, dist_vector = dist_grid$timelag)$gamma
  
  h  = sqrt(dist_grid$spacelag^2 + (model$stAni * dist_grid$timelag)^2)
  Mm = variogramLine(model$joint, covariance = TRUE, dist_vector = h)$gamma
  
  cbind(dist_grid, "gamma" = Sm + Tm + Mm)
}

################################
## simplified sumMetric model ##
################################

vgmSimpleSumMetric <- function(model, dist_grid) {
  vs = variogramLine(model$space, dist_vector=dist_grid$spacelag)$gamma
  vt = variogramLine(model$time,  dist_vector=dist_grid$timelag)$gamma
  
  h = sqrt(dist_grid$spacelag^2 + (model$stAni * as.numeric(dist_grid$timelag))^2)
  
  vm = variogramLine(model$joint, dist_vector=h)$gamma
  vn <- variogramLine(vgm(model$nugget, "Nug", 0), dist_vector=h)$gamma
  
  cbind(dist_grid, "gamma" = vs + vt + vm + vn)
}

covSimpleSumMetric <- function(x, y, model) {
  modelNew <- vgmST("sumMetric", 
                    space=model$space, 
                    time=model$time,
                    joint=vgm(model$joint$psill[model$joint$model != "Nug"],
                              model$joint$model[model$joint$model != "Nug"],
                              model$joint$range[model$joint$model != "Nug"], 
                              model$nugget),
                    stAni=model$stAni)
  if (!is.null(attr(model,"temporal unit")))
    attr(modelNew,"temporal unit") <- attr(model,"temporal unit")
  covSumMetric(x, y, modelNew) 
}

# covariance for the circulant embedding in ST
covSurfSimpleSumMetric <- function(model, dist_grid) {
  modelNew <- vgmST("sumMetric", 
                    space=model$space, 
                    time=model$time,
                    joint=vgm(model$joint$psill[model$joint$model != "Nug"],
                              model$joint$model[model$joint$model != "Nug"],
                              model$joint$range[model$joint$model != "Nug"], 
                              model$nugget),
                    stAni=model$stAni)
  
  if (!is.null(attr(model,"temporal unit")))
    attr(modelNew,"temporal unit") <- attr(model,"temporal unit")
  
  covSurfSumMetric(modelNew, dist_grid) 
}

##################
## metric model ##
##################

vgmMetric <- function(model, dist_grid) {
  h = sqrt(dist_grid$spacelag^2 + (model$stAni * as.numeric(dist_grid$timelag))^2)

  cbind(dist_grid, "gamma" = variogramLine(model$joint, dist_vector=h)$gamma)
}

covMetric <- function(x, y, model) {
  stopifnot(inherits(x, c("STF", "STS", "STI", "sftime")) &&
			inherits(y, c("STF", "STS", "STI", "sftime")))
  
  # the STF case
  if (inherits(x, "STF") && inherits(y, "STF")) {
    # calculate all spatial and temporal distances
    ds = spDists(x@sp, y@sp)
    dt = abs(outer(index(x@time), index(y@time), "-"))
    if(!is.null(attr(model,"temporal unit")))
      units(dt) <- attr(model, "temporal unit") # ensure the same temporal metric as in the variogram definition
    debug_time_unit(units(dt))
    dt <- as(dt, "matrix")
    
    # compose the cov-matrix
    h  = sqrt((matrix(1,nrow(dt),ncol(dt)) %x% ds)^2
              + (model$stAni * dt %x% matrix(1,nrow(ds),ncol(ds)))^2)
    Mm = variogramLine(model$joint, covariance = TRUE, dist_vector = h)
    
    return(Mm)
  } 
  
  # the STI case
  if (inherits(x, c("STI", "sftime")) || inherits(y, c("STI", "sftime"))) {
    # make sure that now both are of type STI
    x <- as(x, "STI")
    y <- as(y, "STI")
    
    # calculate all spatial and temporal distances
    ds = spDists(x@sp, y@sp)
    dt = abs(outer(index(x@time), index(y@time), "-"))
    if(!is.null(attr(model,"temporal unit")))
      units(dt) <- attr(model, "temporal unit") # ensure the same temporal metric as in the variogram definition
    debug_time_unit(units(dt))
    dt <- as(dt, "matrix")
    
    # compose the cov-matrix
    h  = sqrt(ds^2 + (model$stAni * dt)^2)
    Mm = variogramLine(model$joint, covariance = TRUE, dist_vector = h)
    
    return(Mm)
  }
  
  # the remaining cases, none of x and y is STI nor are both STF
  # make sure both are of type STS
  x <- as(x, "STS")
  y <- as(y, "STS")
  
  # calculate all spatial and temporal distances
  ds = spDists(x@sp, y@sp)
  dt = abs(outer(index(x@time), index(y@time), "-"))
  if(!is.null(attr(model,"temporal unit")))
    units(dt) <- attr(model, "temporal unit") # ensure the same temporal metric as in the variogram definition
  debug_time_unit(units(dt))
  dt <- as(dt, "matrix")
  
  # re-arrange the spatial and temporal distances
  sMat <- matrix(NA, nrow(x@index), nrow(y@index))
  tMat <- matrix(NA, nrow(x@index), nrow(y@index))
  for(r in 1:nrow(x@index)) {
    sMat[r,] <- ds[x@index[r,1], y@index[,1]]
    tMat[r,] <- dt[x@index[r,2], y@index[,2]]
  }
  
  # compose the cov-matrix
  h  = sqrt(sMat^2 + (model$stAni * tMat)^2)
  Mm = variogramLine(model$joint, covariance = TRUE, dist_vector = h)
  
  return(Mm)
}

# covariance for the circulant embedding in ST
covSurfMetric <- function(model, dist_grid) {
  h  = sqrt(dist_grid$spacelag^2 + (model$stAni * dist_grid$timelag)^2)
  
  cbind(dist_grid, "gamma" = variogramLine(model$joint, covariance = TRUE, dist_vector = h)$gamma)
}

###########################
## fitting ST variograms ##
###########################

fit.StVariogram <- function(object, model, ..., method = "L-BFGS-B", lower, upper, fit.method = 6, 
                            stAni=NA, wles) {
  stopifnot(inherits(object, "StVariogram"), 
			inherits(model, "StVariogramModel"))
  
  sunit <- attr(object$spacelag, "units")
  tunit <- attr(object$timelag, "units")
  tu.obj = attr(model, "temporal unit")
  if (!is.null(tu.obj))
    stopifnot(identical(tunit, tu.obj))
  
  object$timelag = as.numeric(object$timelag) # needed for R 4.1
  object <- na.omit(object)
  
  ret <- model
  
  if(!missing(wles)) {
    if (wles)
      fit.method = 1
    else
      fit.method = 6
  }
  
  if (fit.method == 0) {
    attr(ret,"optim.output") <- "no fit"
    attr(ret, "MSE") <- mean((object$gamma - variogramSurface(model,
                                                              data.frame(spacelag=object$dist, timelag=object$timelag))$gamma)^2)
    attr(ret, "spatial unit")  <- sunit
    attr(ret, "temporal unit") <- tunit
    
    return(ret)
  }
  
  if ((fit.method == 7 || fit.method == 11) && is.null(model$stAni) && is.na(stAni)) {
    message("[An uninformed spatio-temporal anisotropy value of '1 (spatial unit)/(temporal unit)' is automatically selected. Consider providing a sensible estimate for stAni or using a different fit.method.]")
    stAni <- 1
  }
  
  weightingFun <- switch(fit.method,
                         function(obj, ...) obj$np, # 1
                         function(obj, gamma, ...) obj$np/gamma^2, # 2
                         function(obj, ...) obj$np, # 3
                         function(obj, gamma, ...) obj$np/gamma^2, # 4
                         function(obj, ...) stop("fit.method = 5 (REML), is not yet implemented"), # 5
                         function(obj, ...) 1, # 6
                         function(obj, curStAni, ...) 
                           if(is.na(stAni))
                             obj$np/(obj$dist^2+(curStAni*obj$timelag)^2)
                         else
                           obj$np/(obj$dist^2+(stAni*obj$timelag)^2), # 7
                         function(obj, ...) {
                           dist <- obj$dist
                           dist[dist == 0] <- min(dist[dist != 0], na.rm = TRUE)
                           obj$np/dist^2 # 8, pure space, 0 dist = min (dist > 0)
                         },
                         function(obj, ...) {
                           dist <- obj$timelag
                           dist[dist == 0] <- min(dist[dist != 0], na.rm = TRUE)
                           obj$np/dist^2
                         }, # 9, pure time
                         function(obj, gamma, ...) 1/gamma^2, # 10
                         function(obj, curStAni, ...) {
                           if(is.na(stAni))
                             1/(obj$dist^2+(curStAni*obj$timelag)^2)
                           else
                             1/(obj$dist^2+(stAni*obj$timelag)^2)
                         }, # 11
                         function(obj, ...) {
                           dist <- obj$dist
                           dist[dist == 0] <- min(dist[dist != 0], na.rm = TRUE)
                           1/(obj$dist^2) # 12, pure space
                         },
                         function(obj, ...) {
                           dist <- obj$timelag
                           dist[dist == 0] <- min(dist[dist != 0], na.rm = TRUE)
                           1/(obj$timelag^2)
                         }) # 13, pure time
  
  if(is.null(weightingFun))
    stop(paste("fit.method =", fit.method, "is not implementend"))
  
  fitFun = function(par, trace = FALSE, ...) {
    curModel <- insertPar(par, model)
    gammaMod <- variogramSurface(curModel,
                                 data.frame(spacelag=object$dist,
                                            timelag=object$timelag))$gamma
    resSq <- (object$gamma - gammaMod)^2
    resSq <- resSq * weightingFun(object, gamma=gammaMod, curStAni=curModel$stAni)
    if (trace)
      print(c(par, MSE = mean(resSq)))
    mean(resSq) # seems numerically more well behaved
  }
  
  if(missing(lower)) {
    min.s <- min(object$dist[object$dist>0])*0.05 # 5 % of the minimum distance larger 0
    min.t <- min(object$dist[object$timelag>0])*0.05 # 5 % of the minimum time lag 0),
    pos <- sqrt(.Machine$double.eps) # at least positive
    lower <- switch(strsplit(model$stModel, "_")[[1]][1],
                    separable=c(min.s, 0, min.t, 0, 0),
                    productSum=c(0, min.s, 0, 
                                 0, min.t, 0,
                                 pos),
                    productSumOld=c(0, min.s, 0, 
                                    0, min.t, 0, 0),
                    sumMetric=c(0, min.s, 0, 
                                0, min.t, 0,
                                0, pos, 0, pos),
                    simpleSumMetric=c(0, min.s,
                                      0, min.t,
                                      0, pos, 0, 0, pos),
                    metric=c(0, pos, 0, pos),
                    stop("Only \"separable\", \"productSum\", \"sumMetric\", \"simpleSumMetric\" and \"metric\" are implemented."))
  }
  if(missing(upper))
    upper <- switch(strsplit(model$stModel, "_")[[1]][1],
                    separable=c(Inf, 1, Inf, 1, Inf),
                    productSum=Inf,
                    productSumOld=Inf,
                    sumMetric=Inf,
                    simpleSumMetric=Inf,
                    metric=Inf,
                    stop("Only \"separable\", \"productSum\", \"sumMetric\", \"simpleSumMetric\" and \"metric\" are implemented."))
  
  
  
  pars.fit <- optim(extractPar(model), fitFun, ..., method = method, lower = lower, upper = upper)
  
  ret <- insertPar(pars.fit$par, model)
  attr(ret,"optim.output") <- pars.fit
  attr(ret, "MSE") <- mean((object$gamma - variogramSurface(insertPar(pars.fit$par, model),
                                                            data.frame(spacelag=object$dist, timelag=object$timelag))$gamma)^2)
  attr(ret, "spatial unit")  <- sunit
  attr(ret, "temporal unit") <- tunit
  
  return(ret)
}

###########
## tools ##
###########

# insert parameters into models
insertPar <- function(par, model) {
  switch(strsplit(model$stModel, "_")[[1]][1],
         separable=insertParSeparable(par, model),
         productSum=insertParProdSum(par, model),
         productSumOld=insertParProdSumOld(par, model),
         sumMetric=insertParSumMetric(par, model),
         simpleSumMetric=insertParSimpleSumMetric(par,model),
         metric=insertParMetric(par,model),
         stop("Only \"separable\", \"productSum\", \"sumMetric\", \"simpleSumMetric\" and \"metric\" are implemented."))
}

# extract parameters from models
extractPar <- function(model) {
  switch(strsplit(model$stModel, "_")[[1]][1],
         separable=c(range.s=model$space$range[2], nugget.s=model$space$psill[1],
                     range.t=model$time$range[2],  nugget.t=model$time$psill[1],
                     sill= model$sill[[1]]),
         productSumOld=c(sill.s = rev(model$space$psill)[1], range.s = rev(model$space$range)[1],
                         sill.t = rev(model$time$psill)[1],  range.t = rev(model$time$range)[1], 
                         sill=model$sill[[1]], nugget=model$nugget[[1]]),
         productSum=c(sill.s = model$space$psill[2], range.s = model$space$range[2], nugget.s = model$space$psill[1],
                      sill.t = model$time$psill[2],  range.t = model$time$range[2],  nugget.t = model$time$psill[1],
                      k=model$k),
         sumMetric=c(sill.s = model$space$psill[2], range.s = model$space$range[2], nugget.s = model$space$psill[1], 
                     sill.t = model$time$psill[2], range.t = model$time$range[2], nugget.t = model$time$psill[1],
                     sill.st = model$joint$psill[2], range.st = model$joint$range[2], nugget.st = model$joint$psill[1],
                     anis = model$stAni[[1]]),
         # simplified sumMetric model
         simpleSumMetric=c(sill.s = rev(model$space$psill)[1], range.s = rev(model$space$range)[1], 
                           sill.t = rev(model$time$psill)[1], range.t = rev(model$time$range)[1],
                           sill.st = rev(model$joint$psill)[1], range.st = rev(model$joint$range)[1], 
                           nugget = model$nugget[[1]], anis = model$stAni[[1]]),
         metric=c(sill = model$joint$psill[2], range = model$joint$range[2], nugget = model$joint$psill[1],
                  anis = model$stAni[[1]]),
         stop("Only \"separable\", \"productSum\", \"sumMetric\", \"simpleSumMetric\" and \"metric\" are implemented."))
}

# extract names
extractParNames <- function(model) {
  names(extractPar(model))
}

## dedicated insertion functions
################################

# separable model
insertParSeparable <- function(par, model) {
  vgmST("separable",
        space=vgm(1-par[2],as.character(model$space$model[2]),par[1],par[2],
                  kappa=model$space$kappa[2]),
        time= vgm(1-par[4],as.character(model$time$model[2]),par[3],par[4],
                  kappa=model$time$kappa[2]),
        sill=par[5])
}

# product sum model
insertParProdSumOld <- function(par, model) {
  vgmST("productSumOld",
        space=vgm(par[1],as.character(rev(model$space$model)[1]),par[2],
                  kappa=rev(model$space$kappa)[1]),
        time= vgm(par[3],as.character(rev(model$time$model)[1]),par[4],
                  kappa=rev(model$time$kappa)[1]),
        sill=par[5], nugget=par[6])
}

insertParProdSum <- function(par, model) {
  vgmST("productSum",
        space=vgm(par[1],as.character(model$space$model[2]),par[2],par[3],
                  kappa=model$space$kappa[2]),
        time= vgm(par[4],as.character(model$time$model[2]),par[5], par[6],
                  kappa=model$time$kappa[2]),
        k=par[7])
}

# sum metric model
insertParSumMetric <- function(par, model) {
  vgmST("sumMetric",
        space=vgm(par[1],as.character(model$space$model[2]),par[2],par[3],
                  kappa=model$space$kappa[2]),
        time= vgm(par[4],as.character(model$time$model[2]),par[5],par[6],
                  kappa=model$time$kappa[2]),
        joint=vgm(par[7],as.character(model$joint$model[2]),par[8],par[9],
                  kappa=model$joint$kappa[2]),
        stAni=par[10])
}

# simplified sum metric model
insertParSimpleSumMetric <- function(par, model) {
  vgmST("simpleSumMetric",
        space=vgm(par[1],as.character(rev(model$space$model)[1]),par[2],
                  kappa=rev(model$space$kappa)[1]),
        time= vgm(par[3],as.character(rev(model$time$model)[1]),par[4],
                  kappa=rev(model$time$kappa)[1]),
        joint=vgm(par[5],as.character(rev(model$joint$model)[1]),par[6],
                  kappa=rev(model$joint$kappa)[1]),
        nugget=par[7], stAni=par[8])
}

# metric model
insertParMetric <- function(par, model) {
  vgmST("metric",
        joint=vgm(par[1], as.character(model$joint$model[2]), par[2], par[3],
                  kappa=model$joint$kappa[2]),
        stAni=par[4])
}
