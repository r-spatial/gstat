# $Id: uisim.R,v 1.5 2008-09-25 10:26:00 edzer Exp $
library(sp)
library(gstat)

# prediction grid:
data(meuse.grid)
gridded(meuse.grid) = ~x+y

# define variable as dummy data
v = vgm(.25, "Sph", 900)
g = gstat(NULL, "var1", x~1, beta = .5, nmax = 20, model = v, dummy = TRUE)

# simulation of a single variable
out = predict(g, meuse.grid, nsim = 20, indicators = TRUE)
spplot(out)

# simulation of two correlated variables:
v = vgm(.1, "Sph", 900)
g = gstat(g, "var2", x~1, beta = .25, nmax = 20, model = v, dummy = TRUE)

v = vgm(-.1, "Sph", 900)
g = gstat(g, c("var1", "var2"), model = v)
out = predict(g, meuse.grid, nsim = 10, indicators = TRUE, set = list(order = 2))
spplot(out)

# merge all 10 individual simulations into three-group factors:
for (i in 1:10) {
	v1 = paste("var1.sim", i, sep = "")
	v2 = paste("var2.sim", i, sep = "")
	m = cbind(out[[v1]], out[[v2]], 1 - (out[[v1]]+out[[v2]]))
	mout = factor(apply(m, 1, function(x) which(x == 1)))
	if (i == 1)
		out2 = SpatialPixelsDataFrame(as(out, "SpatialPixels"), data.frame(mout))
	else
		out2[[i]] = mout
}
names(out2) = paste("sim", 1:10, sep="")
spplot(out2)

require(RColorBrewer)
spplot(out2, col.regions=brewer.pal(3, "Set2"))
