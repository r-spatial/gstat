# $Id: krige.q,v 1.15 2009-02-20 13:53:38 edzer Exp $

if (!isGeneric("krige"))
	setGeneric("krige", function(formula, locations, ...)
		standardGeneric("krige"))

"krige.locations" <-
function (formula, locations, data = sys.frame(sys.parent()), 
	newdata, model = NULL, ..., beta = NULL, nmax = Inf, nmin = 0, omax = 0,
	maxdist = Inf, block = numeric(0), nsim = 0, indicators = FALSE, 
	na.action = na.pass, debug.level = 1)
{
    g = gstat(formula = formula, locations = locations, data = data, 
		model = model, beta = beta, nmax = nmax, nmin = nmin, omax = omax,
		maxdist = maxdist, ...)
    predict(g, newdata = newdata, block = block, nsim = nsim,
		indicators = indicators, na.action = na.action, debug.level = debug.level)
}
setMethod("krige", c("formula", "formula"), krige.locations)

krige.spatial <- function(formula, locations, newdata, model = NULL, ..., 
	beta = NULL, nmax = Inf, nmin = 0, omax = 0, maxdist = Inf, 
	block = numeric(0), nsim = 0, indicators = FALSE, 
	na.action = na.pass, debug.level = 1)
{
	# locations = coordinates(arg2)
    g = gstat(formula = formula, # locations = locations, 
		data = locations, 
		model = model, beta = beta, nmax = nmax, nmin = nmin, omax = omax,
		maxdist = maxdist, ...)
    predict(g, newdata = newdata, block = block, nsim = nsim,
		indicators = indicators, na.action = na.action, debug.level = debug.level)
}
setMethod("krige", c("formula", "Spatial"), krige.spatial)

setMethod("krige", c("formula", "NULL"),
	function(formula, locations, newdata, ...) { # manual dispatch based on newdata:
		if (inherits(newdata, c("sf", "sfc", "stars")))
			krige.sf(formula, locations, newdata = newdata, ...)
		else
			krige.spatial(formula, locations, newdata = newdata, ...)
	}
)

krige.sf <- function(formula, locations, newdata, ..., nsim = 0) {
	if (!requireNamespace("sf", quietly = TRUE))
		stop("sf required: install that first") # nocov
	if (!requireNamespace("stars", quietly = TRUE))
		stop("stars required: install that first") # nocov
	crs = sf::st_crs(newdata)
	if (!is.null(locations)) {
		stopifnot(sf::st_crs(locations) == sf::st_crs(newdata))
		crs = sf::st_crs(locations)
		if (!isTRUE(sf::st_is_longlat(locations))) {
			sf::st_crs(locations) = sf::NA_crs_
			sf::st_crs(newdata) = sf::NA_crs_# to avoid problems not handled by sp...
		}
		locations = as(locations, "Spatial")
	}
	ret = krige(formula, locations, as(newdata, "Spatial"), ..., nsim = nsim)
	if (gridded(ret)) {
		st = stars::st_as_stars(ret)
		if (nsim > 0)
			st = sim_to_dimension(st, nsim)
		sf::st_set_crs(st, crs)
	} else
		sf::st_set_crs(sf::st_as_sf(ret), crs)
}
setMethod("krige", c("formula", "sf"), krige.sf)

sim_to_dimension = function(st, nsim) {
	nms = names(stars::st_dimensions(st))
	if (length(st) > nsim) {
		nvars = length(st) / nsim
		l = vector("list", nvars)
		vars = unique(sub("([^.]+)\\.[[:alnum:]]+$", "\\1", names(st)))
		for (i in 1:nvars) {
			range = seq((i-1) * nsim + 1, length.out = nsim)
			m = setNames(st[range], paste0("sim", 1:nsim))
			l[[i]] = setNames(stars::st_set_dimensions(merge(m), 
				names = c(nms, "sample")), vars[i])
		}
		do.call(c, l)
	} else {
		if (nsim > 1)
			setNames(stars::st_set_dimensions(merge(st), names = c(nms, "sample")), "var1")
		else
			st
	}
}

setMethod(krige, signature("formula", "ST"),
	function(formula, locations, newdata, model, ...) {
		krigeST(formula, locations, newdata, model,...) 
	}
)

if (!isGeneric("idw"))
	setGeneric("idw", function(formula, locations, ...)
		standardGeneric("idw"))

idw.locations <- function (formula, locations, data = sys.frame(sys.parent()), 
		newdata, nmax = Inf, nmin = 0, omax = 0, maxdist = Inf, 
		block = numeric(0), 
		na.action = na.pass, idp = 2.0, debug.level = 1) {
	krige(formula, locations, data, newdata, nmax = nmax, nmin = nmin,
		omax = omax, maxdist = maxdist, block = block, na.action = na.action,
		set = list(idp = idp), debug.level = debug.level)
}
setMethod("idw", c("formula", "formula"), idw.locations)

idw.spatial <- function (formula, locations, 
		newdata, nmax = Inf, nmin = 0, omax = 0, 
		maxdist = Inf, block = numeric(0), 
		na.action = na.pass, idp = 2.0, debug.level = 1) {
	krige(formula, locations, newdata, nmax = nmax, nmin = nmin, omax = omax,
		maxdist = maxdist, block = block, na.action = na.action,
		set = list(idp = idp), debug.level = debug.level, model = NULL)
}
setMethod("idw", c("formula", "Spatial"), idw.spatial)

idw.sf <- function (formula, locations, 
		newdata, ..., idp = 2.0) {

	if (!requireNamespace("sf", quietly = TRUE))
		stop("sf required: install that first") # nocov
	if (!requireNamespace("stars", quietly = TRUE))
		stop("stars required: install that first") # nocov

	ret = krige(formula, locations, newdata, ..., set = list(idp = idp), model = NULL)
	if (inherits(newdata, c("sf", "sfc")))
		sf::st_as_sf(ret)
	else if (inherits(newdata, "stars"))
		stars::st_as_stars(ret)
	else stop("newdata should be of class sf or stars")
}
setMethod("idw", c("formula", "sf"), idw.sf)

STx2SpatialPoints = function(x, multiplyTimeWith = 1.0) { 
	x = as(geometry(x), "STI")
	t1 = as.numeric(as.POSIXct(index(x@time)))
	t2 = as.numeric(x@endTime)
	time = multiplyTimeWith * (t1 + t2) / 2
	cc = cbind(coordinates(x), time)
	SpatialPoints(cc, proj4string = x@sp@proj4string) 
}

STxDF2SpatialPointsDataFrame = function(x, multiplyTimeWith = 1.0) { 
	pts = STx2SpatialPoints(geometry(x), multiplyTimeWith)
	SpatialPointsDataFrame(pts, x@data)
}

SpatialPointsDataFrame2STxDF = function(x, class, tz = "", 
		origin = as.POSIXct("1970-01-01",tz=tz)) { 
	cc = coordinates(x)
	time = as.POSIXct(cc[,ncol(cc)], tz=tz, origin = origin)
	sp = SpatialPoints(cc[,-ncol(cc)], proj4string = x@sp@proj4string)
	st = as(STI(sp, time), class)
	addAttrToGeom(STI(sp, time), x@data)
}

idw.ST <- function (formula, locations, 
		newdata, nmax = Inf, nmin = 0, omax = 0, 
		maxdist = Inf, block = numeric(0), 
		na.action = na.pass, idp = 2.0, debug.level = 1, 
		multiplyTimeWith = 1.0) {
	stopifnot(ncol(coordinates(locations@sp)) == 2)
	ret = krige(formula, 
		STxDF2SpatialPointsDataFrame(locations, multiplyTimeWith),
		STx2SpatialPoints(newdata, multiplyTimeWith),
		nmax = nmax, nmin = nmin, omax = omax,
		maxdist = maxdist, block = block, na.action = na.action,
		set = list(idp = idp), debug.level = debug.level, model = NULL)
	SpatialPointsDataFrame2STxDF(ret, class(geometry(newdata)))
}
setMethod("idw", c("formula", "ST"), idw.ST)
