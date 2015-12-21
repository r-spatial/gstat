library(sp)
library(gstat)
data(meuse)
coordinates(meuse) = ~x+y
data(meuse.grid)
gridded(meuse.grid) = ~x+y

# NLKrige: non-linear kriging (e.g. log-normal kriging), simulation based.
# arguments:
# formula, data, newdata, vgm: see ?krige
# trans: transformation and back-transformation function;
# summarize: optional, summarize the simulations, see example below;
# nmax, nsim: see ?krige
# density: number of points to discretize each block PROVIDED BLOCK SIZES
#   ARE CONSTANT (note that newdata can also be SpatialPolygons -- if newdata is
#   a grid, cell size is used as block size)
NLKrige = function(formula, data, newdata, vgm, trans = c(log,exp), summarize,
		nmax = 50, nsim = 100, density = 16, ...) {
	# transform target:
	target = as.character(as.list(formula)[[2]])
	data[[target]] = trans[[1]](data[[target]])
	# conditional simulation sample at finer grid:
	finegrid = spsample(newdata, n = density * length(newdata), 
		type = "regular", offset = c(.5,.5))
	sim = krige(formula, data, finegrid, vgm, nmax = nmax, nsim = nsim)
	# back transform ALL simulations:
	sim@data = trans[[2]](sim@data)
	# spatial aggregation of each simulation, taking block MEAN:
	aggr = aggregate(sim, newdata, mean)
	# aggregation to summary:
	if (!missing(summarize)) {
		ret = apply(aggr@data, 1, summarize, ...)
		if (is.matrix(ret))
			aggr@data = data.frame(t(ret))
		else
			aggr@data = data.frame(ret)
	}
	aggr
}

aggr = NLKrige(zinc~1, meuse, meuse.grid, vgm(.5, "Sph", 900, .1), nmax = 10, 
	summarize = mean)
spplot(aggr, main = "expected value of block means")

aggr = NLKrige(zinc~1, meuse, meuse.grid, vgm(.5, "Sph", 900, .1), nmax = 10, 
	summarize = quantile, probs = c(0.025, 0.975))
spplot(aggr, main = "95% CI for block means")
