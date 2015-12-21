# $Id: ossfim.q,v 1.3 2006-02-10 19:01:07 edzer Exp $

"ossfim" <-
function(spacings = 1:5, block.sizes = 1:5, model, nmax = 25, debug = 0)
{
	n = floor(sqrt(nmax)) + 1
	x = 0:(n-1) + .5
	x = sort(c(-x, x))
	ret = matrix(NA, length(spacings) * length(block.sizes), 3)
	r = 1
	for (sp in spacings) {
		for (bl in block.sizes) {
			data.grid = data.frame(expand.grid(x * sp, x * sp),
				z = rep(1, length(x)^2))
			names(data.grid) = c("x", "y", "z")
			gridded(data.grid) = c("x", "y")
			x0 = SpatialPoints(matrix(0, 1, 2))
			kr = krige(z~1, data.grid, x0,
				block = c(bl, bl), model = model, nmax = nmax,
				set = list(debug = debug))
			ret[r, ] = c(sp, bl, sqrt(kr[["var1.var"]][1]))
			r = r + 1
		}
	}
	ret = data.frame(ret)
	names(ret) = c("spacing", "block.size", "kriging.se")
	ret
}
