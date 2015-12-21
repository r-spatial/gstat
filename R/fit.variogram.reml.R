# $Id: fit.variogram.reml.q,v 1.14 2009-11-02 21:33:17 edzer Exp $

"fit.variogram.reml" <-
function (formula, locations, data, model, debug.level = 1, set, degree = 0)
{
    if (missing(formula)) 
        stop("nothing to fit to")
    if (!inherits(formula, "formula"))
		stop("first argument should be of class formula")
	if (!missing(locations)) {
    	if (inherits(locations, "formula"))
			coordinates(data) = locations
		else if (is(locations, "Spatial"))
			data = as(locations, "SpatialPointsDataFrame")
	}
    if (!is(data, "SpatialPointsDataFrame"))
        stop("data should (now) be of class SpatialPointsDataFrame")
    if (missing(model)) 
        stop("no model to fit")
    if (!inherits(model, "variogramModel"))
        stop("model should be of class variogramModel (use vgm)")
    fit.sills = rep(TRUE, length(model$model))
    fit.ranges = rep(FALSE, length(model$model))
    .Call(gstat_init, as.integer(debug.level))
    ret = gstat.formula(formula, data)
    ret$y <- residuals(lm(formula, data))
    .Call(gstat_new_data, as.double(ret$y), as.double(ret$locations),
		as.double(ret$X), as.integer(1), double(0), as.integer(-1), 
		as.integer(0), as.double(-1), as.integer(0), as.integer(1), 
		double(0), double(0), as.integer(degree),
		as.integer(is.projected(data)), as.integer(0),
		as.double(1.0), as.integer(0))
    load.variogram.model(model)
    if (!missing(set))
    	gstat.load.set(set)
    ret = .Call(gstat_fit_variogram, as.integer(5), 
        as.integer(fit.sills), as.integer(fit.ranges))
    .Call(gstat_exit, 0) 
    model$psill = ret[[1]]
    model$range = ret[[2]]
    model
}
