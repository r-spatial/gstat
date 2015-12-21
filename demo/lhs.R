library(sp)
library(gstat)

# roughly follows the case presented in:
# E.J. Pebesma and G.B.M. Heuvelink, 1999. Latin hypercube sampling 
# of Gaussian random fields. Technometrics 41 (4), pp. 303-312.

data(meuse)
data(meuse.grid)
coordinates(meuse) = ~x+y
gridded(meuse.grid) = ~x+y

x <- variogram(log(zinc) ~ 1, meuse)
v <- vgm(.5, "Sph", 800, nug = .1)
v.fit = fit.variogram(x, model = v)
plot(x, model = v.fit)

n = 100
ok = krige(log(zinc)~1, meuse, meuse.grid, v.fit, nmax=40)
sim = krige(log(zinc)~1, meuse, meuse.grid, v.fit, nsim = n, nmax=40)
simo = t(apply(as.data.frame(sim)[1:n], 1, order)) # rank order
nr = nrow(simo) # number of prediction locations
simo = (simo - 1.0 + matrix(runif(n * nr), nr, n))/n
summary(simo) # LHS on uniform [0,1] distribution; back to Gaussian:
lhs = t(apply(cbind(ok$var1.pred, sqrt(ok$var1.var), simo), 1,
	function(x) qnorm(x[-(1:2)], x[1], x[2])))
sim2 = sim
sim2@data = data.frame(lhs)
spplot(sim2[1:10], main = 'lhs', col.regions=bpy.colors())
# verify that simulated and true mean/var are close:
m = apply(lhs, 1, mean)
v = apply(lhs, 1, var)
summary(m - ok$var1.pred)
summary(v - ok$var1.var)
