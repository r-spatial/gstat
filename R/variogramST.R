VgmFillNA <- function(x, boundaries) {
  # pads the sample variogram with NA rows where no data are available.
	n = length(boundaries) - 1
	ix = rep(NA, n)
	#ix[which(1:n %in% findInterval(x$dist, boundaries))] = 1:nrow(x)
	ix[ findInterval(x$dist, boundaries, rightmost.closed = TRUE) ] = 1:nrow(x)
	# x$b = boundaries[-1]
	x[ix,]
}

VgmAverage = function(ret, boundaries) {
	# take out NULL variograms:
	ret = ret[!sapply(ret, is.null)]
	# take care of missing rows...
	ret = lapply(ret, VgmFillNA, 
			boundaries = c(0, 1e-6 * boundaries[2], boundaries[-1]))
	# weighted average:
  # sum np, weighted sum of gamma and dist devided by summed np
	np = apply(do.call(cbind, lapply(ret, function(x) x$np)), 1, sum,
	           na.rm = TRUE)
	gamma = apply(do.call(cbind, lapply(ret, function(x) x$gamma*x$np)), 1, sum,
		na.rm = TRUE)/np
	dist = apply(do.call(cbind, lapply(ret, function(x) x$dist*x$np)), 1, sum,
		na.rm = TRUE)/np
	v = data.frame(np = np, dist = dist, gamma = gamma)
	class(v) = class(ret[[1]])
	attr(v, "boundaries") = attr(ret[[1]], "boundaries")
	v[is.na(v)] = NA
	v
}

StVgmLag = function(formula, data, dt, pseudo, ...) {
  dotLst <- list(...)
	.ValidObs = function(formula, data)
		!is.na(data[[as.character(as.list(formula)[[2]])]])
	d = dim(data)
	ret = vector("list", d[2] - dt)
	if (dt == 0) {
		for (i in 1:d[2]) {
			d0 = data[,i]
			valid = .ValidObs(formula, d0)
			if(sum(valid) <= 1)
				ret[[i]] <- NULL
			else {
				d0 = d0[valid,]
				ret[[i]] = variogram(formula, d0, ...)
			}
		}
	} else {
		for (i in 1:(d[2] - dt)) {
			d1 = data[, i]
			valid1 = .ValidObs(formula, d1)
			d2 = data[, i + dt]
			valid2 = .ValidObs(formula, d2)
			if(sum(valid1)==0 || sum(valid2)==0)
				ret[[i]] <- NULL
			else {
				d1 = d1[valid1,]
				d2 = d2[valid2,]
				obj = gstat(NULL, paste("D", i, sep=""), formula, d1, 
					set = list(zero_dist = 3), beta = 0)
				obj = gstat(obj, paste("D", i+dt, sep=""), formula, d2, 
					beta = 0)
				ret[[i]] = variogram(obj, cross = "ONLY", pseudo = pseudo, ...)
			}
		}
	}
  if(!is.null(dotLst$cloud)) {
    if(dotLst$cloud)
      ret <- do.call("rbind", ret)
      ret$id <- "var1"
      return(ret)
  } else {
	   return(VgmAverage(ret, dotLst$boundaries))
  }
}

variogramST = function(formula, locations, data, ..., tlags = 0:15, cutoff, 
                       width = cutoff/15, boundaries=seq(0,cutoff,width),
                       progress = interactive(), pseudo = TRUE, 
                       assumeRegular=FALSE, na.omit=FALSE, cores = 1) {
	if (missing(data))
		data = locations

	if (inherits(data, "stars")) {
		if (!requireNamespace("stars", quietly = TRUE))
			stop("stars required: install that first") # nocov
		data = as(data, "STFDF")
	}

	if(missing(cutoff)) {
		ll = !is.na(is.projected(data@sp)) && !is.projected(data@sp)
		cutoff <- spDists(t(data@sp@bbox), longlat = ll)[1,2]/3
	}
  
	if (formula[[3]] != 1) { # there is a regression model:
		data$resid = residuals(lm(formula, data, na.action = na.exclude))
		formula = resid ~ 1
	}
	if(is(data, "STIDF"))
		return(variogramST.STIDF(formula, data, tlags, cutoff, width, 
                             boundaries, progress, cores = cores, ...))
  
	stopifnot(is(data, "STFDF") || is(data, "STSDF"))
	it = index(data@time)
	if (assumeRegular || is.regular(zoo(matrix(1:length(it)), order.by = it),
                                  strict = TRUE)) {
		twidth = diff(it)[1]
		tlags = tlags[tlags <= min(max(tlags), length(unique(it)) - 1)]
	} else {
		warning("strictly irregular time steps were assumed to be regular")
		twidth = mean(diff(it))
	}
	obj = NULL
	t = twidth * tlags
	if (progress)
	  pb = txtProgressBar(style = 3, max = length(tlags))
	if (cores == 1) {
		ret = vector("list", length(tlags))
		for (dt in seq(along = tlags)) {
	  		ret[[dt]] = StVgmLag(formula, data, tlags[dt], pseudo = pseudo, 
						  	boundaries = boundaries, ...)
	  		ret[[dt]]$id = paste("lag", dt - 1, sep="")
	  		if (progress)
				setTxtProgressBar(pb, dt)
		}
	} else {
		if (!requireNamespace("future", quietly = TRUE) || 
				!requireNamespace("future.apply", quietly = TRUE))
	  		stop("For parallelization, future and future.apply packages are required")

	  	future::plan('multiprocess', workers = cores)
  		ret <- split(seq(along=tlags), seq(along=tlags))
		ret <- future.apply::future_lapply(X = ret,
				FUN = function(x){
					xx <- StVgmLag(formula, data, tlags[x], pseudo = pseudo, 
					boundaries = boundaries, ...)
					xx$id <- paste("lag", x - 1, sep="")
					if (progress)
						setTxtProgressBar(pb, x)
					return(xx)
				},
				future.seed = NULL # silence warning
			)
	}
	
	if (progress)
		close(pb)
	# add time lag:
	v = do.call(rbind, ret)
	v$timelag = rep(t, sapply(ret, nrow))
	if (is(t, "yearmon"))
		class(v$timelag) = "yearmon"

	b = attr(ret[[min(length(tlags),2)]], "boundaries")
	b = c(0, b[2]/1e6, b[-1])
	# ix = findInterval(v$dist, b) will use all spacelags
	b = b[-2]
	# spacelags = c(0, b[-length(b)] + diff(b)/2) will use all spacelags
	v$spacelag = c(0, b[-length(b)] + diff(b)/2) # spacelags[ix] will use all spacelags

	v$avgDist <- v$dist * v$np
	for (lagId in unique(v$spacelag)) {
	  bool <- v$spacelag == lagId
	  v$avgDist[bool] <- sum(v$avgDist[bool], na.rm = TRUE) / sum(v$np[bool], na.rm = TRUE)
	}

	class(v) = c("StVariogram", "data.frame")
	if(na.omit)
	v <- na.omit(v)

	# setting attributes to allow krigeST to check units
	attr(v$timelag, "units") <- attr(twidth,"units")
	if (isTRUE(!is.projected(data)))
		attr(v$spacelag, "units") = "km"
  
	return(v)
}

## very irregular data
variogramST.STIDF <- function (formula, data, tlags, cutoff, 
                               width, boundaries, progress, 
                               twindow, tunit, cores = 1) {
  ll = !is.na(is.projected(data@sp)) && !is.projected(data@sp)
  
  if (missing(cutoff))
    cutoff <- spDists(t(data@sp@bbox), longlat = ll)[1, 2]/3
  
  m = model.frame(terms(formula), as(data, "data.frame"))
  
  diffTime <- diff(index(data))
  timeScale <- units(diffTime)
  if(missing(tunit))
    warning(paste("The argument 'tunit' is missing: tlags are assumed to be given in ", timeScale, ".",sep=""))
  else {
    stopifnot(tunit %in% c("secs", "mins", "hours", "days", "weeks"))
    units(diffTime) <- tunit
    timeScale <- tunit
  }
  diffTime <- as.numeric(diffTime)
  if (missing(twindow)) {
    twindow <- round(2 * max(tlags, na.rm=TRUE)/mean(diffTime,na.rm=TRUE),0)
  }
    
  nData <- nrow(data)
  
  # re-using the order propertie of the time slot to only store the next "twindow" distances
  numTime <- as.numeric(index(data))
  diffTimeMat <- matrix(NA, nData, twindow)
  
  for (i in 1:nData) { # i <- 1
    diffTimeMat[i,1:min(nData,twindow)] <- cumsum(diffTime[i+0:(min(nData,twindow)-1)])
  }
  
  nSp <- length(boundaries)
  nTp <- length(tlags)
  
  distTp <- matrix(NA, nSp-1, nTp-1)
  distSp <- matrix(NA, nSp-1, nTp-1)
  gamma  <- matrix(NA, nSp-1, nTp-1)
  np     <- matrix(NA, nSp-1 ,nTp-1)
  
  # temporal selection
  if(progress)
    pb <- txtProgressBar(0, nTp-1, 0, style=3)
  for (i in 1:(nTp-1)) { # i <- 1
    ind <- which(diffTimeMat >= tlags[i] & diffTimeMat < tlags[i+1])
    if (length(ind) < 1) 
      next
    tmpInd <- matrix(NA,nrow=length(ind),4)
    tmpInd[,1] <- ind %% nData  # row number
    tmpInd[,2] <- (ind %/% nData)+1 # col number
    if (cores == 1){
      tmpInd[,3] <- apply(tmpInd[,1:2,drop=FALSE], 1, function(x) spDists(data@sp[x[1]], data@sp[x[2]+x[1],]))
    } else {
      if(!requireNamespace("future", quietly = TRUE) || !requireNamespace("future.apply", quietly = TRUE))
        stop("For parallelization, future and future.apply packages are required")
      future::plan("multiprocess", workers = cores)
      tmpInd[,3] <- future.apply::future_apply(X = tmpInd[,1:2,drop=FALSE], MARGIN = 1, 
                                 FUN = function(x) spDists(data@sp[x[1]], data@sp[x[2]+x[1],]),
                                 future.seed = NULL)
    }
    tmpInd[,4] <- diffTimeMat[tmpInd[,1:2, drop=FALSE]]
    
    # spatial selection
    for (j in 1:(nSp-1)) { # j <- 3
      indSp <- which(tmpInd[,3] >= boundaries[j] & tmpInd[,3] < boundaries[j+1])
      if (length(indSp) < 1)
        next
      distSp[j,i] <- mean(tmpInd[indSp,3])
      distTp[j,i] <- mean(tmpInd[indSp,4])
      
      indSp <- cbind(ind[indSp] %% nData, (ind[indSp] %/% nData)+1)
      np[j,i] <- nrow(indSp) # Issue #7, Thanks to Roelof.
      gamma[j,i] <- 0.5*mean((data[,,colnames(m)[1]]@data[indSp[,1],1] - data[,,colnames(m)[1]]@data[indSp[,1]+indSp[,2],1])^2,
                             na.rm=TRUE)
    }
    if(progress)
      setTxtProgressBar(pb, value=i)
  }
  if(progress)
    close(pb)
  
  res <- data.frame(np=as.vector(np),
                    dist=as.vector(distSp),
                    gamma=as.vector(gamma),
                    id=paste("lag",rep(1:(nTp-1),each=nSp-1), sep=""),
                    timelag=rep(tlags[-nTp]+diff(tlags)/2,each=nSp-1),
                    spacelag=rep(boundaries[-nSp]+diff(boundaries)/2, nTp-1))
  
  res$avgDist <- res$dist * res$np
  for (lagId in unique(res$spacelag)) {
    bool <- res$spacelag == lagId
    res$avgDist[bool] <- sum(res$avgDist[bool], na.rm = TRUE) / sum(res$np[bool], na.rm = TRUE)
  }
  
  attr(res$timelag, "units") <- timeScale
  attr(res$spacelag, "units") <- ifelse(ll, "km", "m")
  class(res) <- c("StVariogram", "data.frame")
  
  return(res)
}



## plotting
plot.StVariogram = function(x, model=NULL, ..., col = bpy.colors(), xlab, ylab,
                            map = TRUE, convertMonths = FALSE, as.table = TRUE,
                            wireframe = FALSE, diff = FALSE, all=FALSE) {
	lst = list(...)
	if (!is.null(lst$col.regions))
		col = lst$col.regions
	if (is(x$timelag, "yearmon")) {
		if (convertMonths) {
			x$timelag = as.numeric(x$timelag) * 12
			attr(x$timelag, "units") = "months"
		} else
			attr(x$timelag, "units") = "years"
	}
	if (missing(xlab)) {
		xlab = "distance"
		u =  attr(x$spacelag, "units")
		if (!is.null(u))
			xlab = paste(xlab, " (", u, ")", sep="")
	}
	if (missing(ylab)) {
		ylab = "time lag"
		u = attr(x$timelag, "units")
		if (!is.null(u))
			ylab = paste(ylab, " (", u, ")", sep="")
	}
	x$timelag = as.numeric(x$timelag)
  
  # check for older spatio-temporal variograms and compute avgDist on demand
  if(is.null(x$avgDist)) {
    x$avgDist <- x$dist * x$np
    for (lagId in unique(x$spacelag)) {
      bool <- x$spacelag == lagId
      x$avgDist[bool] <- sum(x$avgDist[bool] / sum(x$np[bool], na.rm = TRUE), na.rm = TRUE)
    }
    
  }
  
  if(!is.null(model)) {
    if (is(model,"StVariogramModel"))
      model <- list(model)
    for (mod in model) {
      slag <- x$avgDist
      slag[slag == 0 & x$timelag == 0] <- sqrt(.Machine$double.eps)
      x[[mod$stModel]] <- variogramSurface(mod, data.frame(spacelag = slag,
                                                           timelag = x$timelag))$gamma
      if (diff)
        x[[mod$stModel]] <- x[[mod$stModel]] - x$gamma
    }
  }
	x0 = x # needed by wireframe()
	
	if (!is.null(model)) {
    modelNames  <- sapply(model, function(x) x$stModel)
    
    if (all && !diff)
      v0 <- x[,c("dist", "id", "avgDist", "timelag")]
    else
      v0 <- NULL
    
    for (i in modelNames)
      v0 <- rbind(v0, x[,c("dist", "id", "avgDist", "timelag")])

    if(all & !diff) {# we also need the sample
      v0$what = factor(c(rep("sample", nrow(x)), rep(modelNames, each=nrow(x))),
                       levels=c("sample", modelNames), ordered = TRUE)
      v0$gamma = c(x$gamma, unlist(x[,modelNames]))
    } else {
      v0$what = factor(rep(modelNames, each=nrow(x)),
                       levels=modelNames, ordered=TRUE)
      v0$gamma = c(unlist(x[,modelNames]))
    }
		
		x = v0
	}
	if (wireframe) { 
		if (!is.null(model)) {
      if (length(model) > 1 || all)
        wireframe(gamma ~ avgDist*timelag | what, 
                  x, drape = TRUE, col.regions = col, 
                  xlab = xlab, ylab = ylab, as.table=as.table, ...)
      else 
        wireframe(as.formula(paste(model[[1]]$stModel,"~ avgDist*timelag")), 
                  x0, drape = TRUE, col.regions = col, 
                  xlab = xlab, ylab = ylab, as.table=as.table, ...)
		} else # without a model, plot only the sample variogram as a wireframe
			wireframe(gamma ~ avgDist * timelag, x0, drape = TRUE, col.regions = col,
				xlab = xlab, ylab = ylab, ...)
	} else if (map) {
		if (!is.null(model))
			f = gamma ~ avgDist + timelag | what
		else
			f = gamma ~ avgDist + timelag
		levelplot(f, x, xlab = xlab, ylab = ylab, col.regions = col, as.table=as.table, ...)
	} else { # not map, not wireplot
		if (!is.null(model))
			f = gamma ~ dist | what
		else
			f = gamma ~ dist
		x$id = factor(x$id, levels=unique(x$id))
		bp = bpy.colors(length(levels(x$id)))
		ps = list(superpose.line=list(col=bp), superpose.symbol=list(col=bp))
		xlim = c(0, max(x$dist) * 1.04)
		if (diff) {
		  xyplot(f, x, groups = x$id, type='b', xlim = xlim,
		         auto.key = list(space = "right"), xlab = xlab, 
		         par.settings = ps, as.table=as.table, ...)
		} else {
  		ylim = c(0, max(x$gamma) * 1.04)
  		xyplot(f, x, groups = x$id, type='b', ylim = ylim, xlim = xlim,
  		       auto.key = list(space = "right"), xlab = xlab, 
  		       par.settings = ps, as.table=as.table, ...)
    }
	}
}

print.StVariogramModel <- function(x, ...) {
  possComp <- c("space", "time", "joint")
  for(comp in possComp[possComp %in% names(x)]) {
    rownames(x[[comp]]) <- 1:nrow(x[[comp]])
    cat(paste(comp,"component: \n"))
    print(x[[comp]], ...)
  }
  
  possAddPar <- c("sill", "nugget", "stAni", "k")
  for(addPar in possAddPar[possAddPar %in% names(x)]) {
    cat(paste(addPar, ": ",x[[addPar]],"\n", sep=""))
  }
}


## guess the spatio-temporal anisotropy without spatio-temporal models
######################################################################

estiStAni <- function(empVgm, interval, method="linear", spatialVgm, temporalVgm, s.range=NA, t.range=NA) {
  if (!is.na(s.range))
    empVgm <- empVgm[empVgm$dist <= s.range,]
  if (!is.na(t.range))
    empVgm <- empVgm[empVgm$timelag <= t.range,]
  
  empVgm$timelag = as.numeric(empVgm$timelag) # in case it is of class difftime, messes up on R 4.1
  switch(method,
         linear = estiStAni.lin(empVgm, interval),
         range = estiAni.range(empVgm, spatialVgm, temporalVgm),
         vgm = estiAni.vgm(empVgm, spatialVgm, interval),
         metric = estiAni.metric(empVgm, spatialVgm, interval),
         stop(paste("Method", method,"is not implemented.")))
}

# linear
estiStAni.lin <- function(empVgm, interval) {
  lmSp <- lm(gamma~dist, empVgm[empVgm$timelag == min(empVgm$timelag),])
  
  optFun <- function(stAni) {
    sqrt(mean((predict(lmSp, newdata = data.frame(dist=empVgm[empVgm$spacelag == min(empVgm$spacelag),]$timelag*stAni)) - empVgm[empVgm$spacelag == min(empVgm$spacelag),]$gamma)^2, na.rm=TRUE))
  
  optimise(optFun, interval, empVgm = empVgm)$minimum  
}

}
# range
estiAni.range <- function(empVgm, spatialVgm, temporalVgm) {
  spEmpVgm <- empVgm[empVgm$timelag == min(empVgm$timelag),]
  class(spEmpVgm) <- c("gstatVariogram","data.frame")
  spEmpVgm <- spEmpVgm[-1,1:3]
  spEmpVgm$dir.hor <- 0
  spEmpVgm$dir.ver <- 0
  
  spatialVgm <- fit.variogram(spEmpVgm, spatialVgm)
  
  tmpEmpVgm <- empVgm[empVgm$spacelag == min(empVgm$timelag),]
  class(tmpEmpVgm) <- c("gstatVariogram","data.frame")
  tmpEmpVgm <- tmpEmpVgm[-1,c("np","timelag","gamma")]
  colnames(tmpEmpVgm) <- c("np", "dist", "gamma")
  tmpEmpVgm$dir.hor <- 0
  tmpEmpVgm$dir.ver <- 0
  
  temporalVgm <- fit.variogram(tmpEmpVgm, temporalVgm)
  
  spatialVgm$range[2]/temporalVgm$range[2]
}

# variograms
estiAni.vgm <- function(empVgm, spatialVgm, interval) {
  spEmpVgm <- empVgm[empVgm$timelag == min(empVgm$timelag),]
  class(spEmpVgm) <- c("gstatVariogram","data.frame")
  spEmpVgm <- spEmpVgm[-1,1:3]
  spEmpVgm$dir.hor <- 0
  spEmpVgm$dir.ver <- 0
  
  spatialVgm <- fit.variogram(spEmpVgm, spatialVgm)
  
  optFun <- function(stAni) {
    sqrt(mean((variogramLine(spatialVgm, dist_vector = empVgm[empVgm$spacelag == min(empVgm$spacelag),]$timelag*stAni)$gamma - empVgm[empVgm$spacelag == min(empVgm$spacelag),]$gamma)^2, na.rm=TRUE))
  }
  
  optimise(optFun, interval)$minimum
}

# metric variogram
estiAni.metric <- function(empVgm, spatialVgm, interval) {
  fit.StVariogram(empVgm, vgmST("metric", joint=spatialVgm, stAni=mean(interval)))$stAni[[1]]
}

