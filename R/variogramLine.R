# $Id: variogramLine.q,v 1.4 2008-08-18 16:32:42 edzer Exp $

"variogramLine" <-
function(object, maxdist, n = 200, min=1.0e-6 * maxdist, dir = c(1,0,0), 
	covariance = FALSE,	..., dist_vector = numeric(0), debug.level = 0)
{
	if (missing(object))
		stop("model is missing");
	if (!inherits(object, "variogramModel"))
		stop("model should be of mode variogramModel (use function vgm)")
	if (length(dist_vector) > 0)
		maxdist = 0.0
	else if (missing(maxdist))
		stop("maxdist or dist_vector needs to be set");
	if (length(dir) != 3)
		stop("dir should be numeric vector of length 3")
	.Call(gstat_init, as.integer(debug.level))
	pars = c(min,maxdist,n,dir)
	load.variogram.model(object, c(0,0)) # loads object into gstat 
	ret = .Call(gstat_variogram_values, as.integer(c(0,0)),
		as.numeric(pars), as.integer(covariance), as.numeric(dist_vector))
	.Call(gstat_exit, 0);
	if (is.matrix(dist_vector))
		matrix(ret[[2]], nrow(dist_vector), ncol(dist_vector))
	else
		data.frame(dist=ret[[1]], gamma=ret[[2]])
}

# Sat Mar 14 15:11:55 CET 2015: removed this:
# "variogram.line" <- function(..., deprecate = TRUE) {
# 	if (deprecate)
# 		cat("variogram.line is DEPRECATED, please use variogramLine instead\n")
# 	variogramLine(...)
# }
