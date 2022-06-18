#############################
## spatio-temporal kriging ##
#############################

debug_time_unit = function(tUnit) {
    # message("Using the following time unit: ", tUnit)
}

STsolve = function(A, b, X) {
  # V = A$T %x% A$S -- a separable covariance; solve A x = b for x
  # kronecker: T %x% S vec(L) = vec(c0) <-->  S L T = c0
  # solve for L: use Y = L T 
  # S Y = c0 -> Y = solve(S, c0)
  # L T = Y -> Tt Lt = Yt -> Lt = solve(Tt, Yt)
  #Tm = chol(A$Tm, LINPACK=TRUE)
  Tm = chol(A$Tm)
  #Sm = chol(A$Sm, LINPACK=TRUE)
  Sm = chol(A$Sm)
  STbacksolve = function(Tm, Cm, Sm) { 
    MyChSolve = function(A, b)
      backsolve(A, forwardsolve(A, b, upper.tri = TRUE, transpose = TRUE))
    # Y = MyChSolve(Sm, Cm)
    # L = MyChSolve(Tm, t(Y))
    # as.vector(t(L))
    as.vector(t(MyChSolve(Tm, t(MyChSolve(Sm, Cm)))))
  }
  # b comes separated:
  ret1 = apply(b$T, 2, function(x1) 
    apply(b$S, 2, function(x2)
      STbacksolve(Tm, matrix(x1 %x% x2, nrow(Sm), nrow(Tm)), Sm)))
  d = dim(ret1)
  dim(ret1) = c(d[1] / ncol(b$S), d[2] * ncol(b$S))
  # X comes full:
  ret2 = apply(X, 2, function(x) 
    STbacksolve(Tm, matrix(x, nrow(Sm), nrow(Tm)), Sm))
  cbind(ret1, ret2)
}

covfn.ST = function(x, y = x, model, ...) {
  switch(strsplit(model$stModel, "_")[[1]][1],
         separable=covSeparable(x, y, model, ...),
         productSumOld=covProdSumOld(x, y, model),
         productSum=covProdSum(x, y, model),
         sumMetric=covSumMetric(x, y, model),
         simpleSumMetric=covSimpleSumMetric(x, y, model),
         metric=covMetric(x, y, model),
         stop(paste("Provided spatio-temporal model (",model$stModel,") is not supported.",sep="")))
}


## krigeST
krigeST <- function(formula, data, newdata, modelList, beta, y, ...,
                    nmax=Inf, stAni=NULL,
                    computeVar = FALSE, fullCovariance = FALSE,
                    bufferNmax=2, progress=TRUE) {
  stopifnot(inherits(modelList, "StVariogramModel") || is.function(modelList))
  to_sftime = FALSE
  return_stars = if (inherits(data, c("stars"))) {
    if (!requireNamespace("sf", quietly = TRUE))
      stop("sf required: install that first") # nocov
    if (!requireNamespace("stars", quietly = TRUE))
      stop("stars required: install that first") # nocov
    if (sf::st_crs(data) != sf::st_crs(newdata))
      warning("CRS for data and newdata are not identical; assign CRS or use st_transform to correct")
    data = as(data, "STFDF")
    newdata = as(newdata, "STFDF")
    TRUE
  } else if (inherits(data, "sftime")) {
    if (sf::st_crs(data) != sf::st_crs(newdata))
      warning("CRS for data and newdata are not identical; assign CRS or use st_transform to correct")
    data = as(data, "STIDF")
    if (inherits(newdata, "stars")) {
		if (length(newdata) == 0)
			newdata$._dummy = 0.
		newdata = as(newdata, "STFDF")
	}
	if (inherits(newdata, "sftime")) {
  		to_sftime = TRUE
		newdata = as(newdata, "STIDF")
	}
    TRUE
  } else {
    if (!identical(data@sp@proj4string@projargs, newdata@sp@proj4string@projargs))
  	  message("please verify that the CRSs of data and newdata are identical, or transform them first to make them identical")
    FALSE
  }
  stopifnot(inherits(data, c("STF", "STS", "STI", "sftime")) && 
			inherits(newdata, c("STF", "STS", "STI", "sftime"))) 
  if (inherits(data, "sftime"))
    data = as(data, "STIDF")
  if (inherits(newdata, "sftime"))
    data = as(data, "STI")
  stopifnot(class(data@time) == class(newdata@time))
  stopifnot(nmax > 0)
  
  tUnitModel <- attr(modelList, "temporal unit")
  tUnitData <- units(abs(outer(index(data@time[1]), index(newdata@time[1]), "-")))
  
  if (is.null(tUnitModel)) {
    warning("The spatio-temporal variogram model does not carry the strongly recommended attribute 'temporal unit'.\n The unit '", tUnitData,
            "' has been assumed. krigeST could not check whether the temporal distances between locations and in the variogram coincide.")
    tUnit <- tUnitData
    attr(modelList, "temporal unit") <- tUnit
  } else {
    tUnit <- tUnitModel
    debug_time_unit(tUnit)
  }
  
  if(nmax < Inf) { # local neighbourhood ST kriging:
    ret = krigeST.local( formula = formula, data = data, 
                         newdata = newdata, modelList = modelList, beta=beta, # y=y, # for later use
                         nmax = nmax, stAni = stAni, 
                         computeVar = computeVar, fullCovariance = fullCovariance, 
                         bufferNmax = bufferNmax, progress = progress)

    if (return_stars) {
	  ret$._dummy = NULL
  	  if (to_sftime)
		  sftime::st_as_sftime(ret)
	  else 
		  stars::st_as_stars(as(ret, "STFDF"))
	} else
      ret
  } else {
    df <- krigeST.df(formula = formula, data = data, newdata = newdata, 
                   modelList = modelList, beta = beta, y = y, 
                   ..., 
                   nmax=nmax, stAni=stAni,
                   computeVar = computeVar, fullCovariance = fullCovariance,
                   bufferNmax = bufferNmax, progress = progress)
  
    # wrapping the predictions in ST*DF again
    if (!fullCovariance) {
      ret = addAttrToGeom(geometry(newdata), df)
      if (return_stars) {
	    ret$._dummy = NULL
  	  	if (to_sftime)
		  ret = sftime::st_as_sftime(ret)
  		else 
		  ret = stars::st_as_stars(as(ret, "STFDF"))
	  }
      ret
    } else
      df
  }
}

krigeST.df <- function(formula, data, newdata, modelList, beta, y, ...,
                       nmax = Inf, stAni = NULL,
                       computeVar = FALSE, fullCovariance = FALSE,
                       bufferNmax = 2, progress = TRUE) {
  
  separate <- length(data) > 1 && length(newdata) > 1 && inherits(data, "STF") && inherits(newdata, "STF")
  
  lst = extractFormula(formula, data, newdata)
  X = lst$X
  x0 = lst$x0
  if (missing(y))
    y = lst$y
  
  if (inherits(modelList, "StVariogramModel")) {
    V = covfn.ST(data, model = modelList, separate=separate)
    v0 = covfn.ST(data, newdata, modelList)
  
    if (is(data,"STSDF"))
      d0 <- data[data@index[1,1], data@index[1,2], drop = FALSE]
    else
      d0 = data[1, 1, drop=FALSE]
    
    c0 = as.numeric(covfn.ST(d0, d0, modelList, separate = FALSE))
  } else {
    V = modelList(data, data, ...)
    v0 = modelList(data, newdata, ...)
    
    if (computeVar) {
      if (is(newdata@sp, "SpatialLines") || is(newdata@sp, "SpatialPolygons"))
        stop("Varying target support (SpatialLines, SpatialPolygons) for kriging variance is not implemented.")
      c0 = as.numeric(modelList(newdata[1, drop=FALSE],
                                newdata[1, drop=FALSE]))
    }
  }
  
  if (!missing(beta)) { # sk:
    skwts = CHsolve(V, v0)
    npts = length(newdata)
    ViX = skwts[,-(1:npts)]
    skwts = skwts[,1:npts]
    if (computeVar)
      var <- c0 - apply(v0*skwts, 2, sum)
  } else {
    if (!is.function(modelList) && (modelList$stModel == "separable" & separate))
      skwts <- STsolve(V, v0, X) # use Kronecker trick
    else 
      skwts <- CHsolve(V, cbind(v0, X))
    npts = length(newdata)
    ViX = skwts[,-(1:npts)]
    skwts = skwts[,1:npts]
    beta = solve(t(X) %*% ViX, t(ViX) %*% y)
    if (computeVar) {
      # get (x0-X'C-1 c0)'(X'C-1X)-1 (x0-X'C-1 c0) -- precompute term 1+3:
      if(is.list(v0)) # in the separable case
        v0 = v0$Tm %x% v0$Sm
      Q = t(x0) - t(ViX) %*% v0
      # suggested by Marius Appel
      var = c0 - apply(v0 * skwts, 2, sum) + apply(Q * CHsolve(t(X) %*% ViX, Q), 2, sum)
      if (fullCovariance) {
        corMat <- cov2cor(covfn.ST(newdata, newdata, modelList))
        var <- corMat*matrix(sqrt(var) %x% sqrt(var), nrow(corMat), ncol(corMat))
        # var = c0 - t(v0) %*% skwts + t(Q) %*% CHsolve(t(X) %*% ViX, Q)
        # return(list(pred=pred, var=var))
      }
    }
  }
  
  pred = x0 %*% beta + t(skwts) %*% (y - X %*% beta)
  
  if(computeVar) {
    if (fullCovariance)
      list(pred=pred, var=var)
    else
      data.frame(var1.pred = pred, var1.var = var)
  } else
    data.frame(var1.pred = pred)
}

# local spatio-temporal kriging
krigeST.local <- function(formula, data, newdata, modelList, beta, nmax, stAni=NULL,
                          computeVar=FALSE, fullCovariance=FALSE, 
                          bufferNmax=2, progress=TRUE) {
  dimGeom <- ncol(coordinates(data))
  if (fullCovariance)
    stop("fullCovariance cannot be returned for local ST kriging")
  
  if(is.null(stAni) & !is.null(modelList$stAni)) {
    stAni <- modelList$stAni
    
    # scale stAni [spatial/temporal] to seconds
    if(!is.null(attr(modelList,"temporal unit")))
      stAni <- stAni/switch(attr(modelList, "temporal unit"),
                            secs=1,
                            mins=60,
                            hours=3600,
                            days=86400,
                            stop("Temporal unit",attr(modelList, "temporal unit"),"not implemented."))
  }
  
  if(is.null(stAni))
    stop("The spatio-temporal model does not provide a spatio-temporal 
         anisotropy scaling nor is the parameter stAni provided. One of 
         these is necessary for local spatio-temporal kriging.")
  
  # check whether the model meets the coordinates' unit
  if(!is.null(attr(modelList, "spatial unit")))
    stopifnot((is.projected(data) & (attr(modelList, "spatial unit") %in% c("km","m"))) | (!is.projected(data) & !(attr(modelList, "spatial unit") %in% c("km","m"))))
  
  if (inherits(data, c("STFDF", "STSDF", "sftime")))
    data <- as(data, "STIDF")
  
  clnd <- class(newdata)
  
  if(inherits(newdata, c("STFDF", "STSDF", "sftime")))
    newdata <- as(newdata, "STIDF")
  if(inherits(newdata, c("STF", "STS")))
    newdata <- as(newdata, "STI")
  
  # from here on every data set is assumed to be STI*
  
  if(dimGeom == 2) {
    df = as(data, "data.frame")[,c(1,2,4)]
    df$time = as.numeric(df$time)*stAni
    
    query = as(newdata, "data.frame")[,c(1,2,4)]
    query$time = as.numeric(query$time)*stAni
  } else {
    df = as(data, "data.frame")[,c(1,2,3,5)]
    df$time = as.numeric(df$time)*stAni
    
    query = as(newdata, "data.frame")[,c(1,2,3,5)]
    query$time = as.numeric(query$time)*stAni
  }
  
  df <- as.matrix(df)
  query <- as.matrix(query)
  
  if (computeVar) {
    res <- data.frame(var1.pred = rep(NA, nrow(query)),
                      var1.var = rep(NA, nrow(query)))
  } else {
    res <- data.frame(var1.pred = rep(NA, nrow(query)))
  }
  
  if(progress)
    pb = txtProgressBar(style = 3, max = nrow(query))
  
  nb <- t(apply(get.knnx(df, query, ceiling(bufferNmax*nmax))[[1]],1,sort))
  
  for (i in 1:nrow(query)) {
    nghbrData <- data[nb[i, ], , drop = FALSE]
    
    if(bufferNmax > 1) {
      nghbrCov <- covfn.ST(nghbrData, newdata[i, , drop = FALSE], modelList)
      nghbrData <- nghbrData[sort(order(nghbrCov, decreasing=T)[1:nmax]), , drop = FALSE]
    }
    
    res[i,] <- krigeST.df(formula, nghbrData, newdata[i, , drop = FALSE],
                          modelList, computeVar=computeVar, beta = beta,
                          fullCovariance=fullCovariance)
    if(progress)
      setTxtProgressBar(pb, i)  
  }
  if(progress)
    close(pb)
  
  if (clnd %in% c("STI", "STS", "STF")) {
    newdata <- addAttrToGeom(newdata, as.data.frame(res))
    newdata <- switch(clnd, 
                      STI = as(newdata, "STIDF"),
                      STS = as(newdata, "STSDF"),
                      STF = as(newdata, "STFDF"))
  } else {
    newdata@data <- cbind(newdata@data, res)
    newdata <- as(newdata, clnd)
  }
  newdata
}

## Area to point kriging

# define a spatio-temporal variogram model FUNCTION that can deal with x and y
# being of ST* with @sp of class SpatialPolygons OR SpatialPoints;
# SpatialGrid/Pixels are coerced to SpatialPolygons -- as in the so case, see
# krige0.R

# sts <- STS(meuse.grid, stf@time, cbind(sort(sample(length(meuse.grid), 500, replace = T)),
#                                        sample(21,500, replace=TRUE)))
# 
# x <- sts
# y <- stf
# model <- separableModel

vgmAreaST = function(x, y = x, model, ndiscrSpace = 16, verbose = FALSE, covariance = TRUE) {
  if (gridded(x@sp))
    x@sp = as(x@sp, "SpatialPolygons")
  if (gridded(y@sp))
    y@sp = as(y@sp, "SpatialPolygons")
  stopifnot(is(x@sp, "SpatialPolygons") || is(x@sp, "SpatialPoints"))
  stopifnot(is(y@sp, "SpatialPolygons") || is(y@sp, "SpatialPoints"))
  stopifnot(is(model, "StVariogramModel"))
  
  # make x and y both of type STS/STSDF
  if ("data" %in% slotNames(x)) {
    x <- as(x, "STSDF") 
  } else {
    x <- as(x, "STS") 
  }

  if ("data" %in% slotNames(y)) {
    y <- as(y, "STSDF") 
  } else {
    y <- as(y, "STS") 
  }
  
  nx = length(x)
  ny = length(y)
  
  V <- matrix(NA, nx, ny)
  
  # switch cases for polygons in x and y 
  if (is(x@sp, "SpatialPolygons")) {
    # x contains polygons -> loop and sample
    
    indexX <- x@index
    if (verbose)
      pb <- txtProgressBar(style = 3, max = nx)
    for (ix in 1:nx) { # ix <- 1
      px <- x[indexX[ix,,drop=F],drop=F]
      ptsx = STF(spsample(px@sp, ndiscrSpace, "regular", offset = c(.5,.5)), 
                 px@time, px@endTime)
      
      # does y also contain polygons?
      if (is(y@sp, "SpatialPolygons")) {
        # yes, also y contains polygons -> loop and sample
        
        indexY <- y@index
        for (iy in 1:ny) { # iy <- 1
          py <- y[indexY[iy,,drop=F],drop=F]
          ptsy = STF(spsample(py@sp, ndiscrSpace, "regular", offset = c(.5,.5)), 
                     py@time, py@endTime)
          
          suppressMessages(subV <- covfn.ST(ptsx, ptsy, model, separate=FALSE))
          
          V[ix, iy] <- mean(subV)
        }
      } else {
        # no, y contains points -> no second loop, calc covariance
        suppressMessages(subV <- covfn.ST(ptsx, y, model, separate=FALSE))
        V[ix, ] <- apply(subV, 2, mean)
      }
      
      if (verbose)
        setTxtProgressBar(pb, ix)
    }
    if (verbose)
      close(pb)
  } else {
    # none of x and y contain polygons -> no loops, calc covariance
    suppressMessages(V <- covfn.ST(x, y, model))
  }

  V
}

####################
## trans Gaussian ##
####################

krigeSTTg <- function(formula, data, newdata, modelList, y, nmax=Inf, stAni=NULL,
                      bufferNmax=2, progress=TRUE, lambda = 0) {
  if(!is.infinite(nmax))
    return(krigeSTTg.local(formula, data, newdata, modelList, y, nmax, stAni,
                           bufferNmax, progress, lambda))
  
  lst <- extractFormula(formula, data, newdata)
  
  Y <- lst$y
  X <- lst$X
  
  if (ncol(X) > 1)
    stop("only formula with intercept allowed, e.g. y ~ 1")
  
  data$value = phiInv(Y, lambda)
  data$value1 = rep(1, length(data$value))
  
  OK = krigeST(value ~ 1, data, newdata, modelList, 
               nmax = nmax, stAni = stAni, 
               computeVar=TRUE, bufferNmax = bufferNmax, progress = progress)
  
  separate <- length(data) > 1 && length(newdata) > 1 && 
    inherits(data, "STF") && inherits(newdata, "STF")
  
  V <- covfn.ST(data, model = modelList, separate=separate)
  Vi <- solve(V)
  
  muhat <- sum(Vi %*% data$value)/sum(Vi)
  
  # find m:
  v0 <- covfn.ST(data, newdata, model = modelList, separate=separate)
  m <- (1 - apply(Vi %*% v0,2,sum))/sum(Vi)
  
  # compute transGaussian kriging estimate & variance:
  OK$var1TG.pred = phi(OK$var1.pred, lambda) + 
    phiDouble(muhat, lambda) * (OK$var1.var/2 - m)
  OK$var1TG.var = phiPrime(muhat, lambda)^2 * OK$var1.var
  OK
}



krigeSTTg.local <- function(formula, data, newdata, modelList, y, nmax=Inf, stAni=NULL,
                            bufferNmax=2, progress=TRUE, lambda = 0) {
  stopifnot(!is.infinite(nmax))
  dimGeom <- ncol(coordinates(data))
  
  if(is.null(stAni) & !is.null(modelList$stAni)) {
    stAni <- modelList$stAni
    
    # scale stAni [spatial/temporal] to seconds
    if(!is.null(attr(modelList,"temporal unit")))
      stAni <- stAni/switch(attr(modelList, "temporal unit"),
                            secs=1,
                            mins=60,
                            hours=3600,
                            days=86400,
                            stop("Temporal unit",attr(modelList, "temporal unit"),"not implemented."))
  }
  
  if(is.null(stAni))
    stop("The spatio-temporal model does not provide a spatio-temporal 
         anisotropy scaling nor is the parameter stAni provided. One of 
         these is necessary for local spatio-temporal kriging.")
  
  # check whether the model meets the coordinates' unit
  if(!is.null(attr(modelList, "spatial unit")))
    stopifnot((is.projected(data) & (attr(modelList, "spatial unit") %in% c("km","m"))) | (!is.projected(data) & !(attr(modelList, "spatial unit") %in% c("km","m"))))
  
  if(is(data, "STFDF") || is(data, "STSDF"))
    data <- as(data, "STIDF")
  
  clnd <- class(newdata)
  
  if(is(newdata, "STFDF") || is(newdata, "STSDF"))
    newdata <- as(newdata, "STIDF")
  if(is(newdata, "STF") || is(newdata, "STS"))
    newdata <- as(newdata, "STI")
  
  # from here on every data set is assumed to be STI*
  
  if(dimGeom == 2) {
    df = as(data, "data.frame")[,c(1,2,4)]
    df$time = as.numeric(df$time)*stAni
    
    query = as(newdata, "data.frame")[,c(1,2,4)]
    query$time = as.numeric(query$time)*stAni
  } else {
    df = as(data, "data.frame")[,c(1,2,3,5)]
    df$time = as.numeric(df$time)*stAni
    
    query = as(newdata, "data.frame")[,c(1,2,3,5)]
    query$time = as.numeric(query$time)*stAni
  }
  
  df <- as.matrix(df)
  query <- as.matrix(query)
  
  res <- data.frame(var1.pred = rep(NA, nrow(query)),
                    var1.var = rep(NA, nrow(query)),
                    var1TG.pred = rep(NA, nrow(query)),
                    var1TG.var = rep(NA, nrow(query)))
  
  if(progress)
    pb = txtProgressBar(style = 3, max = nrow(query))
  
  nb <- t(apply(get.knnx(df, query, ceiling(bufferNmax*nmax))[[1]],1,sort))
  
  for (i in 1:nrow(query)) {
    nghbrData <- data[nb[i, ], , drop = FALSE]
    
    if(bufferNmax > 1) {
      nghbrCov <- covfn.ST(nghbrData, newdata[i, , drop = FALSE], modelList)
      nghbrData <- nghbrData[sort(order(nghbrCov, decreasing=T)[1:nmax]), , drop = FALSE]
    }
    
    res[i,] <- krigeSTTg(formula, nghbrData, newdata[i, , drop = FALSE],
                         modelList, lambda)@data
    if(progress)
      setTxtProgressBar(pb, i)  
  }
  if(progress)
    close(pb)
  
  if (clnd %in% c("STI", "STS", "STF")) {
    newdata <- addAttrToGeom(newdata, as.data.frame(res))
    newdata <- switch(clnd, 
                      STI = as(newdata, "STIDF"),
                      STS = as(newdata, "STSDF"),
                      STF = as(newdata, "STFDF"))
  } else {
    newdata@data <- cbind(newdata@data, res)
    newdata <- as(newdata, clnd)
  }
  newdata
}
