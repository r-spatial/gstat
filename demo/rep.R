library(sp)
library(gstat)
data(meuse)
coordinates(meuse) = ~x + y
data(meuse.grid)
coordinates(meuse.grid) = ~x + y

# Variogram log Zn

lzn.vgm = variogram(log(zinc) ~ 1, meuse)
lzn.fit = fit.variogram(lzn.vgm, model = vgm(1, "Sph", 900, 1))

#Conditional simulation

nsim = 100
lzn.sim = krige(log(zinc) ~ 1, meuse, meuse.grid, model = lzn.fit,
	nmax = 30, nsim = nsim)

# Variogram of all relizations

m = out = list()
for (i in 1:nsim) {
	s = paste("sim", i, sep="")
	f = as.formula(paste(s, "~1"))
	v = variogram(f, lzn.sim)
	v$id = rep(s, nrow(v))
	out[[s]] = v
	m[[s]] = fit.variogram(v, lzn.fit)
}
plot(do.call(rbind, out), m, layout=c(10,10), skip = FALSE, 
	scales = list(y = list(relation = "same")))
