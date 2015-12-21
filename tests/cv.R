# try bivariate cokriging; cross validate first variable
library(sp)
data(meuse)
library(gstat)
g=gstat(NULL, "log-zinc", log(zinc)~1,  ~x+y, meuse, nmax=10)
g=gstat(g, "log-lead", log(lead)~1,     ~x+y, meuse, nmax=10)
g=gstat(g, "log-copper", log(copper)~1, ~x+y, meuse, nmax=10)
v=variogram(g)
g=gstat(g, model=vgm(1, "Sph", 1000), fill.all=T)
g=fit.lmc(v, g)
g
set.seed(13131)
summary(gstat.cv(g, remove.all=TRUE, nfold=5))
summary(gstat.cv(g, remove.all=FALSE, nfold=5))
