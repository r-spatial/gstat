# $Id: gstat.formula.predict.q,v 1.14 2008-02-19 10:01:22 edzer Exp $

"gstat.formula.predict" <-
function (formula, newdata, na.action, BLUE.estimates = FALSE, xlev = NULL) 
{
	if (is(newdata, "SpatialPolygons")) {
		# locs = coordinates(getSpatialPolygonsLabelPoints(newdata)) -- deprecated, now use:

		locs = t(sapply(slot(newdata, "polygons"), function(x) slot(x, "labpt")))
		SpatialPoints(locs, CRS(proj4string(newdata)))
		locs = coordinates(locs)

		colnames(locs) = c("x", "y")
		if (is(newdata, "SpatialPolygonsDataFrame"))
			newdata = as.data.frame(newdata)
		else
			newdata = data.frame(a = rep(1, nrow(locs)))
	} else if (is(newdata, "SpatialLines")) {
		# locs = coordinates(getSpatialLinesMidPoints(newdata)) -- deprecated, now use:

		ret = lapply(newdata@lines,
	        	function(x) sapply(x@Lines,
				function(X) apply(X@coords, 2, mean)
			)
		)
		ret = t(sapply(ret, function(x) apply(x, 1, mean)))
		locs = coordinates(SpatialPoints(ret, CRS(proj4string(newdata))))
		colnames(locs) = c("x", "y")

		if (is(newdata, "SpatialLinesDataFrame"))
			newdata = as.data.frame(newdata)
		else
			newdata = data.frame(a = rep(1, nrow(locs)))
	} else {
		if (gridded(newdata))
			fullgrid(newdata) = FALSE
		locs = coordinates(newdata)
		newdata = as.data.frame(newdata)
	} 

	# resolve formula:
	terms.f = delete.response(terms(formula))
    mf.f = model.frame(terms.f, newdata, na.action = na.action, xlev = xlev)
    X = model.matrix(terms.f, mf.f)

	if (BLUE.estimates) { # fake the whole thing to get a matrix with BLUE parameter estimates:
		cnames = colnames(X)
		X = matrix(0, ncol(X), ncol(X))
		diag(X) = 1
		locs = locs[1,,drop=FALSE]
		if (ncol(X) > 1) {
			for (i in 2:ncol(X))
				locs = rbind(locs, locs[1,])
		}
		rownames(locs) = cnames
	}

	if (NROW(locs) != NROW(X)) { 
		# NA's were filtered in X, but not in coords:
    	mf.f =    model.frame(terms.f, newdata, na.action = na.pass)
		valid.pattern = !(apply(mf.f, 1, function(x) any(is.na(x))))
		X    = model.matrix(terms.f, mf.f[valid.pattern, , drop = FALSE])
		locs = locs[valid.pattern, ]
		if (NROW(locs) != NROW(X))
			stop("NROW(locs) != NROW(X): this should not occur")
	}
    list(locations = as.matrix(locs), X = as.matrix(X))
}
