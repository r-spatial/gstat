# constructiong spatio-temporal variogram models
vgmST <- function(stModel, ..., space, time, joint, sill, k, nugget, stAni, 
		temporalUnits) {
	stopifnot(is.character(stModel) && length(stModel)==1)
  
  old.stModel <- stModel
  stModel <- strsplit(stModel, "_")[[1]][1]
  
  if (stModel == "productSum" & !missing(sill))
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
		productSumOld = list(space = space, time = time, sill = sill, nugget = nugget),
		sumMetric = list(space = space, time = time, joint = joint, stAni = stAni),
		simpleSumMetric = list(space = space, time = time, 
			joint = joint, nugget = nugget, stAni = stAni),
		metric = list(joint = joint, stAni = stAni),
		stop(paste("model", stModel, "unknown")))
	vgmModel$stModel <- old.stModel
	if (!missing(temporalUnits))
		attr(vgmModel, "temporal units") = temporalUnits
	class(vgmModel) <- c("StVariogramModel", "list")
	vgmModel
}

# calculating spatio-temporal variogram surfaces
variogramSurface <- function(model, dist_grid, ...) {
  if (!inherits(model, "StVariogramModel"))
    warning("\"model\" should be of class \"StVariogramModel\"; no further checks for a proper model will made.")
  
  switch(strsplit(model$stModel, "_")[[1]][1],
         separable=vgmSeparable(model, dist_grid, ...),
         productSum=vgmProdSum(model, dist_grid, ...),
         productSumOld=vgmProdSumOld(model, dist_grid, ...),
         sumMetric=vgmSumMetric(model, dist_grid, ...),
         simpleSumMetric=vgmSimpleSumMetric(model, dist_grid, ...),
         metric=vgmMetric(model, dist_grid, ...),
         stop("Only \"separable\", \"productSum\", \"sumMetric\", \"simpleSumMetric\" and \"metric\" are implemented."))
}

# separable model: C_s * C_t
vgmSeparable <- function(model, dist_grid) {
  vs = variogramLine(model$space, dist_vector=dist_grid$spacelag)[,2]
  vt = variogramLine(model$time,  dist_vector=dist_grid$timelag)[,2]

  data.frame(spacelag=dist_grid$spacelag, timelag=dist_grid$timelag, 
             model=model$sill*(vs+vt-vs*vt))
}

# productSum model: C_s*C_t + C_s + C_t
vgmProdSumOld <- function(model, dist_grid) {
  warning("Please consider to re-fit your model for the new product-sum notation.")
  vs = variogramLine(model$space, dist_vector=dist_grid$spacelag)[,2]
  vt = variogramLine(model$time, dist_vector=dist_grid$timelag)[,2]
  vn <- rep(model$nugget, length(vs))
  vn[vs == 0 & vt == 0] <- 0

  k <- (sum(model$space$psill)+sum(model$time$psill)-(model$sill+model$nugget))/(sum(model$space$psill)*sum(model$time$psill))
  
  if (k <= 0 | k > 1/max(rev(model$space$psill)[1], rev(model$time$psill)[1])) 
    k <- 10^6*abs(k) # distorting the model to let optim "hopefully" find suitable parameters
  data.frame(spacelag=dist_grid$spacelag, timelag=dist_grid$timelag, 
             model=as.vector(vs+vt-k*vs*vt+vn))
}

vgmProdSum <- function(model, dist_grid) {
  if(!is.null(model$sill)) # backwards compatibility
    vgmProdSumOld(model, dist_grid)
  vs = variogramLine(model$space, dist_vector=dist_grid$spacelag)[,2]
  vt = variogramLine(model$time, dist_vector=dist_grid$timelag)[,2]
  
  sill_s <- sum(model$space$psill)
  sill_t <- sum(model$time$psill)
  k <- model$k
  
  data.frame(spacelag=dist_grid$spacelag, timelag=dist_grid$timelag,
             model=as.vector((k*sill_t+1)*vs + (k*sill_s+1)*vt-k*vs*vt))
}

# sumMetric model: C_s + C_t + C_st (Gerard Heuvelink)
vgmSumMetric <- function(model, dist_grid) {
  vs = variogramLine(model$space, dist_vector=dist_grid$spacelag)[,2]
  vt = variogramLine(model$time,  dist_vector=dist_grid$timelag)[,2]
  h = sqrt(dist_grid$spacelag^2 + (model$stAni * as.numeric(dist_grid$timelag))^2)
  vst = variogramLine(model$joint, dist_vector=h)[,2]
  data.frame(spacelag=dist_grid$spacelag, timelag=dist_grid$timelag, model=(vs + vt + vst))
}

# simplified sumMetric model
vgmSimpleSumMetric <- function(model, dist_grid) {
  vs = variogramLine(model$space, dist_vector=dist_grid$spacelag)[,2]
  vt = variogramLine(model$time,  dist_vector=dist_grid$timelag)[,2]
  h = sqrt(dist_grid$spacelag^2 + (model$stAni * as.numeric(dist_grid$timelag))^2)
  vm = variogramLine(model$joint, dist_vector=h)[,2]
  vn <- variogramLine(vgm(model$nugget, "Nug", 0), dist_vector=h)[,2]
  data.frame(spacelag=dist_grid$spacelag, timelag=dist_grid$timelag, model=(vs + vt + vm + vn))
}

vgmMetric <- function(model, dist_grid) {
  h = sqrt(dist_grid$spacelag^2 + (model$stAni * as.numeric(dist_grid$timelag))^2)
  vm = variogramLine(model$joint, dist_vector=h)[,2]
  data.frame(spacelag=dist_grid$spacelag, timelag=dist_grid$timelag, model=vm)
}

fit.StVariogram <- function(object, model, ..., method = "L-BFGS-B", lower, upper, fit.method = 6, 
		stAni=NA, wles) {
  if (!inherits(object, "StVariogram"))
    stop("\"object\" must be of class \"StVariogram\"")
  if (!inherits(model, "StVariogramModel"))
    stop("\"model\" must be of class \"StVariogramModel\".")

  sunit <- attr(object$spacelag, "units")
  tunit <- attr(object$timelag, "units")
  tu.obj = attr(model, "temporal units")
  if (!is.null(tu.obj))
  	stopifnot(identical(tunit, tu.obj))
  
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
                                                              data.frame(spacelag=object$dist, timelag=object$timelag))$model)^2)
    attr(ret, "spatial unit")  <- sunit
    attr(ret, "temporal unit") <- tunit
    
    return(ret)
  }
    
  if ((fit.method == 7 | fit.method == 11) & is.null(model$stAni) & is.na(stAni)) {
    warning("An uninformed spatio-temporal anisotropy value of '1 (spatial unit)/(temporal unit)' is automatically selected. Consider providing a sensible estimate for stAni or using a different fit.method.")
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
    gammaMod <- variogramSurface(curModel, data.frame(spacelag=object$dist,
                                                      timelag=object$timelag))$model
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
                                                            data.frame(spacelag=object$dist, timelag=object$timelag))$model)^2)
  attr(ret, "spatial unit")  <- sunit
  attr(ret, "temporal unit") <- tunit

  return(ret)
}

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

extractParNames <- function(model) {
  names(extractPar(model))
}

insertParSeparable <- function(par, model) {
  vgmST("separable",
        space=vgm(1-par[2],as.character(model$space$model[2]),par[1],par[2],
                  kappa=model$space$kappa[2]),
        time= vgm(1-par[4],as.character(model$time$model[2]),par[3],par[4],
                  kappa=model$time$kappa[2]),
        sill=par[5])
}

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

# simplified sumMetric model
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

insertParMetric <- function(par, model) {
  vgmST("metric",
        joint=vgm(par[1], as.character(model$joint$model[2]), par[2], par[3],
                  kappa=model$joint$kappa[2]),
        stAni=par[4])
}

## guess the spatio-temporal anisotropy without spatio-temporal models

estiStAni <- function(empVgm, interval, method="linear", spatialVgm, temporalVgm, s.range=NA, t.range=NA) {
  if (!is.na(s.range))
    empVgm <- empVgm[empVgm$dist <= s.range,]
  if (!is.na(t.range))
    empVgm <- empVgm[empVgm$timelag <= t.range,]
  
  switch(method,
         linear = estiStAni.lin(empVgm, interval),
         range = estiAni.range(empVgm, spatialVgm, temporalVgm),
         vgm = estiAni.vgm(empVgm, spatialVgm, interval),
         metric = estiAni.metric(empVgm, spatialVgm, interval),
         stop(paste("Method", method,"is not implemented.")))
}

estiStAni.lin <- function(empVgm, interval) {
  lmSp <- lm(gamma~dist, empVgm[empVgm$timelag == 0,])
  
  optFun <- function(stAni) {
    sqrt(mean((predict(lmSp, newdata = data.frame(dist=empVgm[empVgm$spacelag == 0,]$timelag*stAni)) - empVgm[empVgm$spacelag == 0,]$gamma)^2, na.rm=TRUE))
  }
  
  optimise(optFun, interval)$minimum  
}


estiAni.range <- function(empVgm, spatialVgm, temporalVgm) {
  spEmpVgm <- empVgm[empVgm$timelag == 0,]
  class(spEmpVgm) <- c("gstatVariogram","data.frame")
  spEmpVgm <- spEmpVgm[-1,1:3]
  spEmpVgm$dir.hor <- 0
  spEmpVgm$dir.ver <- 0
  
  spatialVgm <- fit.variogram(spEmpVgm, spatialVgm)
  
  tmpEmpVgm <- empVgm[empVgm$spacelag == 0,]
  class(tmpEmpVgm) <- c("gstatVariogram","data.frame")
  tmpEmpVgm <- tmpEmpVgm[-1,c("np","timelag","gamma")]
  colnames(tmpEmpVgm) <- c("np", "dist", "gamma")
  tmpEmpVgm$dir.hor <- 0
  tmpEmpVgm$dir.ver <- 0
  
  temporalVgm <- fit.variogram(tmpEmpVgm, temporalVgm)
  
  spatialVgm$range[2]/temporalVgm$range[2]
}

estiAni.vgm <- function(empVgm, spatialVgm, interval) {
  spEmpVgm <- empVgm[empVgm$timelag == 0,]
  class(spEmpVgm) <- c("gstatVariogram","data.frame")
  spEmpVgm <- spEmpVgm[-1,1:3]
  spEmpVgm$dir.hor <- 0
  spEmpVgm$dir.ver <- 0
  
  spatialVgm <- fit.variogram(spEmpVgm, spatialVgm)
  
  optFun <- function(stAni) {
    sqrt(mean((variogramLine(spatialVgm, dist_vector = empVgm[empVgm$spacelag == 0,]$timelag*stAni)$gamma - empVgm[empVgm$spacelag == 0,]$gamma)^2, na.rm=TRUE))
  }
  
  optimise(optFun, interval)$minimum
}

estiAni.metric <- function(empVgm, spatialVgm, interval) {
  fit.StVariogram(empVgm, vgmST("metric", joint=spatialVgm, stAni=mean(interval)))$stAni[[1]]
}

