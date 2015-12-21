# $Id: load.variogram.model.q,v 1.7 2008-11-12 10:04:22 edzer Exp $

"load.variogram.model" <- function(model, ids = c(0, 0), max_dist = rep(-1.0, 3)) {
	if (missing(model))
		stop("model is missing");
	if (!inherits(model, "variogramModel"))
		stop("model should be of mode variogramModel (use function vgm)")
	if (any(model$range < 0.0)) {
		print(model)
		stop("variogram range can never be negative")
	}
	stopifnot(length(max_dist) == 3)
	anis = c(model$ang1, model$ang2, model$ang3, model$anis1, model$anis2)
	if (is.null(attr(model, "table")))
		covtable = numeric(0)
	else  {
		covtable = attr(model, "table")
		if (dim(model)[1] > 1 || model$model != "Tab")
			stop("table can only have one single model")
	}
# max_dist hack here for Lin(0) models:
#	if (max_dist > 0) {
#		w = which(model$model %in% c("Lin") & model$range == 0)
#		if (length(w) > 0) {
#			model[w,"psill"] = max_dist * model[w,"psill"]
#			model[w,"range"] = max_dist
#			cat("Conversion into equivalent model:\n")
#			print(model)
#		}
#	}
	if (!any(model$model %in% c("Lin", "Pow")))
		max_dist = rep(-1.0, 3) # ignore
	.Call(gstat_load_variogram, 
		as.integer(ids),
		as.character(model$model),
		as.numeric(model$psill),
		as.numeric(model$range),
		as.numeric(model$kappa),
		as.numeric(anis),
		covtable,
		as.numeric(max_dist))
}
