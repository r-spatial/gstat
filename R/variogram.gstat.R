# $Id: variogram.gstat.q,v 1.9 2007-04-06 11:29:58 edzer Exp $

"variogram.gstat" = function (object, ...) {
	if (!inherits(object, "gstat"))
		stop("first argument should be of class gstat")
	y = list()
	locations = list()
	X = list()
	beta = list()
	grid = list()
	projected = TRUE
	for (i in seq(along = object$data)) {
		d = object$data[[i]]
		beta[[i]] = d$beta
		if (i > 1 && !identical(proj4string(object$data[[1]]$data), proj4string(d$data)))
			stop("data items in gstat object have different coordinate reference systems")
		raw = gstat.formula(d$formula, d$data)
		y[[i]] = raw$y
		locations[[i]] = raw$locations
		X[[i]] = raw$X
		grid[[i]] = raw$grid
		if (is(d$data, "Spatial"))
			projected = is.projected(d$data)
		if (d$degree != 0)
			stop("degree != 0: residual variograms wrt coord trend using degree not supported")
	}
	names(y) = names(locations) = names(X) = names(object$data)
	# call variogram.default() next:
	variogram(y, locations, X, trend.beta = beta, grid = grid, g = object, ...,
		projected = projected)
}
