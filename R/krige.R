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
setMethod("krige", c("formula", "NULL"), krige.spatial)

setMethod(krige, signature("formula", "ST"),
          function(formula, locations, newdata, model, ...) {
            krigeST(formula, locations, newdata, model,...) 
          }
# 	function(formula, locations, newdata, model, ...) {
# 		d = data.frame(krigeST(formula, locations, newdata, model,...))
# 		if (ncol(d) == 1)
# 			names(d) = "var1.pred"
# 		if (ncol(d) == 2)
# 			names(d) = c("var1.pred", "var1.var")
# 		addAttrToGeom(geometry(newdata), d)
# 	}
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

STx2SpatialPoints = function(x, multiplyTimeWith = 1.0) { 
	x = as(geometry(x), "STI")
	t1 = as.numeric(as.POSIXct(index(x@time)))
	t2 = as.numeric(x@endTime)
	time = multiplyTimeWith * (t1 + t2) / 2
	cc = cbind(coordinates(x), time)
	SpatialPoints(cc, proj4string = CRS(proj4string(x))) 
}

STxDF2SpatialPointsDataFrame = function(x, multiplyTimeWith = 1.0) { 
	pts = STx2SpatialPoints(geometry(x), multiplyTimeWith)
	SpatialPointsDataFrame(pts, x@data)
}

SpatialPointsDataFrame2STxDF = function(x, class, tz = "", 
		origin = as.POSIXct("1970-01-01",tz=tz)) { 
	cc = coordinates(x)
	time = as.POSIXct(cc[,ncol(cc)], tz=tz, origin = origin)
	sp = SpatialPoints(cc[,-ncol(cc)], proj4string = CRS(proj4string(x)))
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
