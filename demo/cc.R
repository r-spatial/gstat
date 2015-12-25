library(sp)
demo(meuse, ask = FALSE, echo = FALSE)

library(gstat)

# use collocated data:
g = gstat(NULL, "lzinc", log(zinc)~1, meuse)
g.coll = gstat(g,    "dist", dist~1, meuse, nmax = 1, merge = c("lzinc", "dist"))

g.fit = fit.lmc(variogram(g.coll), g.coll, vgm(1, "Sph", 900, 1),
	correct.diagonal = 1.01)

g.non_coll = gstat(g,  "dist", dist~1, meuse.grid, nmax = 1, merge = c("lzinc", "dist"))
g.non_coll$model = g.fit$model

# collocated cokriging:

pr = predict(g.non_coll, meuse.grid)
spplot(pr[c(1,3)])
