# $Id: fit.variogram.q,v 1.10 2008-12-15 14:27:29 edzer Exp $

"fit.variogram" <-
function (object, model, fit.sills = TRUE, fit.ranges = TRUE, 
    fit.method = 7, debug.level = 1, warn.if.neg = FALSE) 
{
    if (missing(object)) 
        stop("nothing to fit to")
	if (!inherits(object, "gstatVariogram") && !inherits(object, "variogramCloud"))
		stop("object should be of class gstatVariogram or variogramCloud")
	if (inherits(object, "variogramCloud"))
		object$np = rep(1, nrow(object))
	if (length(unique(object$id)) > 1)
		stop("to use fit.variogram, variogram object should be univariable")
    if (missing(model)) 
        stop("no model to fit")
    if (!inherits(model, "variogramModel"))
        stop("model should be of class variogramModel (use vgm)")
    if (fit.method == 5)
    	stop("use function fit.variogram.reml() to use REML")
    if (length(fit.sills) < length(model$model)) 
        fit.sills = rep(fit.sills, length(model$model))
    if (length(fit.ranges) < length(model$model)) 
        fit.ranges = rep(fit.ranges, length(model$model))
	if (fit.method == 7 && any(object$dist == 0))
		stop("fit.method 7 will not work with zero distance semivariances; use another fit.method value")
    fit.ranges = fit.ranges & !(model$model %in% c("Nug", "Err"))
	initialRange = model$range
    .Call(gstat_init, as.integer(debug.level))
    .Call(gstat_load_ev, object$np, object$dist, object$gamma)
    load.variogram.model(model)
    ret = .Call(gstat_fit_variogram, as.integer(fit.method), 
        as.integer(fit.sills), as.integer(fit.ranges))
    .Call(gstat_exit, 0)
    model$psill = ret[[1]]
    model$range = ret[[2]]
	attr(model, "singular") = as.logical(ret[[3]]);
	attr(model, "SSErr") = ret[[4]]
	direct = attr(object, "direct")
	if (!is.null(direct)) {
		id = unique(object$id)
		if (direct[direct$id == id, "is.direct"] && any(model$psill < 0)) {
			if (warn.if.neg)
				warning("partial sill or nugget fixed at zero value")
			fit.sills = model$psill > 0
			model$psill[model$psill < 0] = 0.0
			model$range = initialRange
			return(fit.variogram(object, model, fit.sills = fit.sills, fit.ranges =
				fit.ranges, fit.method = fit.method, debug.level = debug.level,
				warn.if.neg = warn.if.neg))
		}
	}
	if (attr(model, "singular") && debug.level) {
		rat = mean(object$gamma) / mean(object$dist) 
		if (rat > 1e6 || rat < 1e-6)
			print("a possible solution MIGHT be to scale semivariances and/or distances")
	}
    model
}
