# $Id: variogram.formula.q,v 1.8 2006-02-10 19:01:07 edzer Exp $

"variogram.formula" <-
function (object, locations = coordinates(data), data, ...) 
{
	if ((missing(locations) && inherits(data, c("sf", "stars"))) || (inherits(locations, c("sf", "stars")))) {
		if (!requireNamespace("sf", quietly = TRUE))
       		stop("sf required: install that first") # nocov
		if (missing(locations))
			data = as(data, "Spatial")
		else
			locations = as(locations, "Spatial")
	}

	# gstat.formula takes care of the case where locations contains
	# both data and coordinates --- see there.
	## ret = gstat.formula(object, locations, data)
	## variogram(object = ret$y, locations = ret$locations, X = ret$X, ...)
	if ((missing(locations) && is(data, "ST")) || (is(locations, "ST")))
		variogramST(formula = object, locations = locations, data = data, ...)
	else {
		g = gstat(formula = object, locations = locations, data = data)
		variogram(g, ...)
	}
}
