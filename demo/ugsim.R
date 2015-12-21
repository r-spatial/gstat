# $Id: ugsim.R,v 1.2 2006-02-10 19:05:02 edzer Exp $
library(sp)
library(gstat)

# prediction grid:
data(meuse.grid)
gridded(meuse.grid) = ~x+y

# define variable as dummy data (parameters from log-zinc, meuse)
v = vgm(.55, "Sph", 900, .05)
g = gstat(NULL, "var1", lzn~1, beta = 5.9, nmax = 20, model = v, dummy = TRUE)

# simulation of a single variable
out = predict(g, meuse.grid, nsim = 20)
spplot(out)

# simulation of two negatively correlated variables:
v = vgm(.55, "Sph", 900, .05)
g = gstat(g, "var2", x~1, beta = 5.9, nmax = 20, model = v, dummy = TRUE)

v = vgm(-.3, "Sph", 900, 0.00001)
g = gstat(g, c("var1", "var2"), model = v)
out = predict(g, meuse.grid, nsim = 10)
spplot(out)
