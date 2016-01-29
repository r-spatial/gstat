# $Id: vgm.q,v 1.12 2007-02-27 22:09:32 edzer Exp $

warn.angle3 = TRUE

"vgm" <-
function(psill = NA, model, range = NA, nugget, add.to, anis, kappa = 0.5,
		..., covtable, Err = 0) {
	add.to.df = function(x, y) {
		x = rbind(x, y)
		row.names(x) = 1:nrow(x)
		return(x)
	}
	if (is.character(psill)) { # as of ``vgm("Sph")''
		if (length(psill) > 1) {
			ret = lapply(psill, function(x) vgm(x))
			class(ret) = c("variogramModelList", "list")
			return(ret)
		} # else:
		if (psill == "Nug")
			return(vgm(NA, "Nug", NA))
		else
			return(vgm(NA, psill, NA, NA))
	}
	stopifnot(length(psill) == 1)
	stopifnot(length(range) == 1)
	stopifnot(missing(nugget) || length(nugget) == 1)
	stopifnot(length(kappa) == 1)
	m = .Call(gstat_get_variogram_models, as.integer(0))
	n = length(m)
	mf = factor(m, levels = m)
	if (missing(model)) {
		ml = .Call(gstat_get_variogram_models, as.integer(1))
		mlf = factor(ml, levels = ml)
		return(data.frame(short = mf, long = mlf))
	}
	if (length(model) > 1) {
		ret = lapply(model, function(x) vgm(,x))
		class(ret) = c("variogramModelList", "list")
		return(ret)
	} 
	table = NULL
	if (model == "Tab" && !missing(covtable)) {
		table = as.matrix(covtable)
		if (NCOL(table) != 2)
			stop("covtable should be a 2-column matrix with distance and cov.")
		range = max(table[,1])
		if (min(table[,1]) != 0.0)
			stop("the first covariance value should be at distance 0.0")
		table = table[,2]
		mf = factor(c(m, "Tab"), levels = c(m, "Tab"))
		if (!missing(add.to) || !missing(nugget) || Err > 0)
			stop("cannot add submodels or nugget to covariance Table model")
	} else if (!any(m == model)) 
		stop(paste("variogram model", model, "unknown\n"))
	if (missing(anis))
		anis = c(0,0,0,1,1)
	if (length(anis) == 2)
		anis = c(anis[1], 0, 0, anis[2], 1)
	else if (length(anis) != 5)
		stop("anis vector should have length 2 (2D) or 5 (3D)")
	if (warn.angle3 && anis[3] != 0.0) {
		warn.angle3 = FALSE
		warning("you are using the third rotation angle; this code is based on the GSLIB2 code\nand must contain the bug described at the end of http://pangea.stanford.edu/ERE/research/scrf/software/gslib/bug/")
	}
	if (!is.na(range)) {
		if (model != "Nug") {
			if (model != "Lin" && model != "Err" && model != "Int")
				if (range <= 0.0) stop("range should be positive")
			else if(range < 0.0) stop("range should be non-negative")
			} else {
			if (range != 0.0) stop("Nugget should have zero range")
			if (anis[4] != 1.0 || anis[5] != 1.0)
				stop("Nugget anisotropy is not meaningful")
		}
	}
	if (!missing(nugget)) {
		ret = data.frame(model=mf[mf==model], psill=psill, range=range,
			kappa = kappa, ang1=anis[1], ang2=anis[2], ang3=anis[3], 
			anis1=anis[4], anis2=anis[5])
		n.vgm = data.frame(model=mf[mf=="Nug"], psill=nugget, range=0,
			kappa = 0.0, ang1=0.0, ang2=0.0, ang3=0.0, anis1=1.0, anis2=1.0)
		ret = add.to.df(n.vgm, ret)
	} else
		ret = data.frame(model=mf[mf==model], psill=psill, range=range,
			kappa = kappa, ang1=anis[1], ang2=anis[2], ang3=anis[3], 
			anis1=anis[4], anis2=anis[5])
	if (!missing(add.to))
		ret = add.to.df(data.frame(add.to), ret)
	if (Err > 0)
		ret = add.to.df(data.frame(model=mf[mf=="Err"], psill=Err,
			range=0.0, kappa = kappa, ang1=anis[1], ang2=anis[2],
			ang3=anis[3], anis1=anis[4], anis2=anis[5]), ret)
	if (!is.null(table))
		attr(ret, "table") = table
	class(ret) = c("variogramModel", "data.frame")
	ret
}

as.vgm.variomodel = function(m) {
	model = NULL
	if (m$cov.model == "exponential")
		model = "Exp"
	else if (m$cov.model == "circular")
		model = "Cir"
	else if (m$cov.model == "gaussian")
		model = "Gau"
	else if (m$cov.model == "linear")
		# model = "Lin"
		stop("no correct conversion available; use power model with power 1?")
	else if (m$cov.model == "matern")
		model = "Mat"
	else if (m$cov.model == "wave")
		model = "Wav"
	else if (m$cov.model == "power")
		model = "Pow"
	else if (m$cov.model == "spherical")
		model = "Sph"
	else if (m$cov.model == "pure.nugget")
		return(vgm(m$nugget + m$cov.pars[1], "Nug", 0))
	else
		stop("variogram model not supported")
# "cauchy",
#,"cubic",
#            "gneiting",
#			"gneiting.matern",
#			"powered.exponential",
#			"wave") ) {
	vgm(m$cov.pars[1], model, m$cov.pars[2], m$nugget, kappa = m$kappa)
}

#$lambda
#[1] 1

#$trend
#[1] "cte"

#$max.dist
#[1] 1441.83

#attr(,"class")
#[1] "variomodel"

