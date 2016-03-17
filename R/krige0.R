extractFormula = function(formula, data, newdata) {
	# extract y and X from data:
    m = model.frame(terms(formula), as(data, "data.frame"), na.action = na.fail)
    y = model.extract(m, "response")
    if (length(y) == 0)
        stop("no response variable present in formula")
    Terms = attr(m, "terms")
    X = model.matrix(Terms, m)

	# extract x0 from newdata:
    terms.f = delete.response(terms(formula))
    mf.f = model.frame(terms.f, newdata) #, na.action = na.action)
    x0 = model.matrix(terms.f, mf.f)
	list(y = y, X = X, x0 = x0)
}

idw0 = function(formula, data, newdata, y, idp = 2.0) {
	s = coordinates(data)
	s0 = coordinates(newdata)
	if (missing(y))
		y = extractFormula(formula, data, newdata)$y
	D = 1.0 / (spDists(s0, s) ^ idp)
	sumD = apply(D, 1, sum)
	D %*% y / sumD
}

CHsolve = function(A, b) {
	# solves A x = b for x if A is PD symmetric
	#A = chol(A, LINPACK=TRUE) -> deprecated
	A = chol(A) # but use pivot=TRUE?
	backsolve(A, forwardsolve(A, b, upper.tri = TRUE, transpose = TRUE))
}

krige0 <- function(formula, data, newdata, model, beta, y, ..., 
		computeVar = FALSE, fullCovariance = FALSE) {

	stopifnot(identical(proj4string(data), proj4string(newdata)))
	lst = extractFormula(formula, data, newdata)
	X = lst$X
	x0 = lst$x0
	if (missing(y))
		y = lst$y
	ll = (!is.na(is.projected(data)) && !is.projected(data))
	s = coordinates(data)
	s0 = coordinates(newdata)
	if (is(model, "variogramModel")) {
		require(gstat)
		V = variogramLine(model, dist_vector = spDists(s, s, ll),
			covariance = TRUE)
		v0 = variogramLine(model, dist_vector = spDists(s, s0, ll),
			covariance = TRUE)
		c0 = variogramLine(model, dist_vector = c(0), covariance = TRUE)$gamma
	} else {
		V = model(data, data, ...)
		v0 = model(data, newdata, ...)
		if (computeVar) {
			if (is(newdata, "SpatialLines") || is(newdata, "SpatialPolygons"))
				stop("varying target support (SpatialLines, SpatialPolygons) is not implemented")
			c0 = as.numeric(model(newdata[1, drop=FALSE],
				newdata[1, drop=FALSE]))
			# ?check this: provide TWO arguments, so model(x,y) can target
			# eventually Y, instead of measurements Z=Y+e
			# with e measurement error term e
		}
	}
	if (!missing(beta)) { # sk:
		skwts = CHsolve(V, v0)
		if (computeVar)
			var <- c0 - apply(v0*skwts, 2, sum)
	} else { # ok/uk -- need to estimate beta:
		skwts = CHsolve(V, cbind(v0, X))
		ViX = skwts[,-(1:nrow(s0))]
		skwts = skwts[,1:nrow(s0)]
		beta = solve(t(X) %*% ViX, t(ViX) %*% y)
		if (computeVar) { 
			Q = t(x0) - t(ViX) %*% v0
			var <- c0 - apply(v0*skwts, 2, sum) + 
				apply(Q * CHsolve(t(X) %*% ViX, Q), 2, sum)
		}
	}
	pred = x0 %*% beta + t(skwts) %*% (y - X %*% beta)
	if (computeVar) {
		if (fullCovariance) {
			corMat <- cov2cor(variogramLine(model, 
				dist_vector = spDists(s0, s0), covariance = TRUE))
			var <- corMat*matrix(sqrt(var) %x% sqrt(var), 
				nrow(corMat), ncol(corMat))
		}
		list(pred = pred, var = var)
	} else
		pred
}

krigeST <- function(formula, data, newdata, modelList, y, beta, nmax=Inf, stAni=NULL,
                       computeVar = FALSE, fullCovariance = FALSE,
                       bufferNmax=2, progress=TRUE) {
  stopifnot(inherits(modelList, "StVariogramModel"))
  stopifnot(inherits(data, c("STF", "STS", "STI")) & inherits(newdata, c("STF", "STS", "STI"))) 
  stopifnot(identical(proj4string(data@sp), proj4string(newdata@sp)))
  stopifnot(class(data@time) == class(newdata@time))
  stopifnot(nmax > 0)
  
  if(is.null(attr(modelList,"temporal unit")))
    warning("The spatio-temporal variogram model does not carry a time unit attribute: krisgeST cannot check whether the temporal distance metrics coincide.")
  
  if(nmax < Inf) # local neighbourhood ST kriging:
	return(krigeST.local(formula = formula, data = data, 
                         newdata = newdata, modelList = modelList, beta=beta, nmax = nmax, 
                         stAni = stAni, computeVar = computeVar, 
                         fullCovariance = fullCovariance, 
                         bufferNmax = bufferNmax, progress))
  
  df <- krigeST.df(formula=formula, data=data, newdata=newdata, 
                   modelList=modelList, y=y, beta=beta, nmax=nmax, stAni=stAni,
                   computeVar = computeVar, fullCovariance = fullCovariance,
                   bufferNmax=bufferNmax, progress=progress)
  
  
  # wrapping the predictions in ST*DF again
  if (!fullCovariance)
  	addAttrToGeom(geometry(newdata), df)
  else
  	df
}
  
  
krigeST.df <- function(formula, data, newdata, modelList, y, beta, nmax=Inf, stAni=NULL,
                    computeVar = FALSE, fullCovariance = FALSE,
                    bufferNmax=2, progress=TRUE) {

	separate <- length(data) > 1 && length(newdata) > 1 && 
		inherits(data, "STF") && inherits(newdata, "STF")
  
	lst = extractFormula(formula, data, newdata)
	X = lst$X
	x0 = lst$x0
	if (missing(y))
		y = lst$y

	V = covfn.ST(data, model = modelList, separate=separate)
	v0 = covfn.ST(data, newdata, modelList)

	if (is(data,"STSDF"))
		d0 <- data[data@index[1,1], data@index[1,2], drop = FALSE]
	else
    	d0 = data[1, 1, drop=FALSE]
	
	c0 = as.numeric(covfn.ST(d0, d0, modelList, separate = FALSE))

	if (!missing(beta)) { # sk:
	  skwts = CHsolve(V, v0)
	  npts = length(newdata)
	  ViX = skwts[,-(1:npts)]
	  skwts = skwts[,1:npts]
	  if (computeVar)
	    var <- c0 - apply(v0*skwts, 2, sum)
	} else {
  	if (modelList$stModel == "separable" & separate)
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
  	    return(list(pred=pred, var=var))
  	  }
  	}
	}
	
	pred = x0 %*% beta + t(skwts) %*% (y - X %*% beta)
	
  if(computeVar)
		return(data.frame(var1.pred = pred, var1.var = var))
  else
  	return(data.frame(var1.pred = pred))
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

subsetThroughDfInd <- function(data, ind) {
  switch(class(data),
         STSDF = data[cbind(data@index[ind,1],data@index[ind,2]),drop=F],
#          STFDF = data[(ind-1)%%length(data@sp)+1,(ind-1)%/%length(data@sp)+1,drop=F],
         STIDF = data[ind,,drop=F],
         STS = data[cbind(data@index[ind,1],data@index[ind,2]),drop=F],
#          STF = data[(ind-1)%%length(data@sp)+1,(ind-1)%/%length(data@sp)+1,drop=F],
         STI = data[ind,,drop=F])
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

## covariance models
####################

covSeparable <- function(x, y, model, separate) {  
  if(missing(separate))
    separate <- inherits(x, "STF") & inherits(y, "STF") & length(x) > 1 & length(y) > 1

  # the STF case
  if (inherits(x, "STF") && inherits(y, "STF")) {
    # calculate all spatial and temporal distances
    ds = spDists(x@sp, y@sp)
    dt = abs(outer(index(x@time), index(y@time), "-"))
    if(!is.null(attr(model,"temporal unit")))
      units(dt) <- attr(model, "temporal unit") # ensure the same temporal metric as in the variogram definition
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
  if (inherits(x, "STI") || inherits(y, "STI")) {
    # make sure that now both are of type STI
    x <- as(x, "STI")
    y <- as(y, "STI")
    
    # calculate all spatial and temporal distances
    ds = spDists(x@sp, y@sp)
    dt = abs(outer(index(x@time), index(y@time), "-"))
    if(!is.null(attr(model,"temporal unit")))
      units(dt) <- attr(model, "temporal unit") # ensure the same temporal metric as in the variogram definition
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

## old product-sum model, BG
covProdSumOld <- function(x, y, model) {
  stopifnot(inherits(x, c("STF", "STS", "STI")) & inherits(y, c("STF", "STS", "STI")))
  
  # double check model for validity, i.e. k:
  k <- (sum(model$space$psill)+sum(model$time$psill)-model$sill)/(sum(model$space$psill)*sum(model$time$psill))
  if (k <= 0 | k > 1/max(model$space$psill[model$space$model!="Nug"], 
                         model$time$psill[model$time$model!="Nug"]))
    stop(paste("k (",k,") is non-positive or too large: no valid model!",sep=""))
  
  # the STF case
  if (inherits(x, "STF") & inherits(y, "STF")) {
    # calculate all spatial and temporal distances
    ds = spDists(x@sp, y@sp)
    dt = abs(outer(index(x@time), index(y@time), "-"))
    if(!is.null(attr(model,"temporal unit")))
      units(dt) <- attr(model, "temporal unit") # ensure the same temporal metric as in the variogram definition
    dt <- as(dt, "matrix")
    
    # compose the cov-matrix
    vs = variogramLine(model$space, dist_vector = ds)
    vt = variogramLine(model$time, dist_vector = dt)
     
    return(model$sill-(vt %x% matrix(1,nrow(vs),ncol(vs)) + matrix(1,nrow(vt),ncol(vt)) %x% vs - k * vt %x% vs))
  } 

  # the STI case
  if(inherits(x, "STI") | inherits(y, "STI")) {
    # make sure that now both are of type STI
    x <- as(x, "STI")
    y <- as(y, "STI")
    
    # calculate all spatial and temporal distances
    ds = spDists(x@sp, y@sp)
    dt = abs(outer(index(x@time), index(y@time), "-"))
    if(!is.null(attr(model,"temporal unit")))
      units(dt) <- attr(model, "temporal unit") # ensure the same temporal metric as in the variogram definition
    dt <- as(dt, "matrix")
    
    # compose the cov-matrix
    vs = variogramLine(model$space, dist_vector = ds)
    vt = variogramLine(model$time, dist_vector = dt)
    
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
  dt <- as(dt, "matrix")
  
  # re-arrange the spatial and temporal distances
  sMat <- matrix(NA, nrow(x@index), nrow(y@index))
  tMat <- matrix(NA, nrow(x@index), nrow(y@index))
  for(r in 1:nrow(x@index)) {
    sMat[r,] <- ds[x@index[r,1], y@index[,1]]
    tMat[r,] <- dt[x@index[r,2], y@index[,2]]
  }
  
  # compose the cov-matrix
  vs = variogramLine(model$space, dist_vector = sMat)
  vt = variogramLine(model$time, dist_vector = tMat)
  
  return(model$sill-(vt + vs - k * vt * vs))
}

## new product-sum model

covProdSum <- function(x, y, model) {
  stopifnot(inherits(x, c("STF", "STS", "STI")) & inherits(y, c("STF", "STS", "STI")))
  if(!is.null(model$sill)) # backwards compatibility
    covProdSumOld(x, y, model)

  # the STF case
  if (inherits(x, "STF") & inherits(y, "STF")) {
    # calculate all spatial and temporal distances
    ds = spDists(x@sp, y@sp)
    dt = abs(outer(index(x@time), index(y@time), "-"))
    if(!is.null(attr(model,"temporal unit")))
      units(dt) <- attr(model, "temporal unit") # ensure the same temporal metric as in the variogram definition
    dt <- as(dt, "matrix")
    
    # compose the cov-matrix
    vs = variogramLine(model$space, dist_vector = ds, covariance =TRUE)
    vt = variogramLine(model$time, dist_vector = dt, covariance =TRUE)
    
    return(vt %x% matrix(1,nrow(vs),ncol(vs)) + matrix(1,nrow(vt),ncol(vt)) %x% vs + model$k * vt %x% vs)
  } 
  
  # the STI case
  if(inherits(x, "STI") | inherits(y, "STI")) {
    # make sure that now both are of type STI
    x <- as(x, "STI")
    y <- as(y, "STI")
    
    # calculate all spatial and temporal distances
    ds = spDists(x@sp, y@sp)
    dt = abs(outer(index(x@time), index(y@time), "-"))
    if(!is.null(attr(model,"temporal unit")))
      units(dt) <- attr(model, "temporal unit") # ensure the same temporal metric as in the variogram definition
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

## sumMetric model
covSumMetric <- function(x, y, model) {
  stopifnot(inherits(x, c("STF", "STS", "STI")) & inherits(y, c("STF", "STS", "STI")))
  
  # the STF case
  if (inherits(x, "STF") & inherits(y, "STF")) {
    # calculate all spatial and temporal distances
    ds = spDists(x@sp, y@sp)
    dt = abs(outer(index(x@time), index(y@time), "-"))
    if(!is.null(attr(model,"temporal unit")))
      units(dt) <- attr(model, "temporal unit") # ensure the same temporal metric as in the variogram definition
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
  if(inherits(x, "STI") | inherits(y, "STI")) {
    # make sure that now both are of type STI
    x <- as(x, "STI")
    y <- as(y, "STI")
    
    # calculate all spatial and temporal distances
    ds = spDists(x@sp, y@sp)
    dt = abs(outer(index(x@time), index(y@time), "-"))
    if(!is.null(attr(model,"temporal unit")))
      units(dt) <- attr(model, "temporal unit") # ensure the same temporal metric as in the variogram definition
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

## simple sumMetric model
covSimpleSumMetric <- function(x, y, model) {
  modelNew <- vgmST("sumMetric", space=model$space, time=model$time,
	joint=vgm(sum(model$joint$psill),
	model$joint$model[model$joint$model != "Nug"],
	model$joint$range, model$nugget), stAni=model$stAni)
  if (!is.null(attr(model,"temporal unit")))
    attr(modelNew,"temporal unit") <- attr(model,"temporal unit")
  covSumMetric(x, y, modelNew) 
}

## metric model
covMetric <- function(x, y, model) {
  stopifnot(inherits(x, c("STF", "STS", "STI")) & inherits(y, c("STF", "STS", "STI")))
  
  # the STF case
  if (inherits(x, "STF") & inherits(y, "STF")) {
    # calculate all spatial and temporal distances
    ds = spDists(x@sp, y@sp)
    dt = abs(outer(index(x@time), index(y@time), "-"))
    if(!is.null(attr(model,"temporal unit")))
      units(dt) <- attr(model, "temporal unit") # ensure the same temporal metric as in the variogram definition
    dt <- as(dt, "matrix")
    
    # compose the cov-matrix
    h  = sqrt((matrix(1,nrow(dt),ncol(dt)) %x% ds)^2
              + (model$stAni * dt %x% matrix(1,nrow(ds),ncol(ds)))^2)
    Mm = variogramLine(model$joint, covariance = TRUE, dist_vector = h)
    
    return(Mm)
  } 
  
  # the STI case
  if(inherits(x, "STI") | inherits(y, "STI")) {
    # make sure that now both are of type STI
    x <- as(x, "STI")
    y <- as(y, "STI")
    
    # calculate all spatial and temporal distances
    ds = spDists(x@sp, y@sp)
    dt = abs(outer(index(x@time), index(y@time), "-"))
    if(!is.null(attr(model,"temporal unit")))
      units(dt) <- attr(model, "temporal unit") # ensure the same temporal metric as in the variogram definition
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

# define variogram model FUNCTION that can deal with x and y
# being of class SpatialPolygons OR SpatialPoints; SpatialGrid/Pixels are coerced to SpatialPolygons
vgmArea = function(x, y = x, vgm, ndiscr = 16, verbose = FALSE, covariance = TRUE) {
	if (gridded(x))
		x = as(x, "SpatialPolygons")
	if (gridded(y))
		y = as(y, "SpatialPolygons")
	stopifnot(is(x, "SpatialPolygons") || is(x, "SpatialPoints"))
	stopifnot(is(y, "SpatialPolygons") || is(y, "SpatialPoints"))
	stopifnot(is(vgm, "variogramModel"))
	nx = length(x)
	ny = length(y)
	V = matrix(NA, nx, ny)
	if (verbose)
		pb = txtProgressBar(style = 3, max = nx)
	for (i in 1:nx) {
		if (is(x, "SpatialPolygons"))
			px = spsample(x[i,], ndiscr, "regular", offset = c(.5,.5))
		else
			px = x[i,]
		for (j in 1:ny) {
			if (is(y, "SpatialPolygons"))
				py = spsample(y[j,], ndiscr, "regular", offset = c(.5,.5))
			else
				py = y[j,]
			D = spDists(px, py)
			D[D == 0] = 1e-10
			V[i,j] = mean(variogramLine(vgm, dist_vector = D, 
				covariance = covariance))
		}
		if (verbose)
			setTxtProgressBar(pb, i)
	}
	if (verbose)
		close(pb)
	V
}

### trans Gaussian

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
