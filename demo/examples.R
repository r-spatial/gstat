# $Id: examples.R,v 1.6 2006-02-10 19:05:02 edzer Exp $
## ex01.cmd, ex02.cmd:
##
## Two variables with (initial estimates of) variograms,
## calcute sample variogram and plot fitted model
##

library(sp)
par(ask = TRUE)
data(meuse)
coordinates(meuse)=~x+y
x <- variogram(zinc ~ 1, meuse)
v <- vgm(140000, "Sph", 800, nug = 10000)
plot(x, model = v)
plot(x, model = fit.variogram(x, model = v))
x <- variogram(log(zinc) ~ 1, meuse)
v <- vgm(.5, "Sph", 800, nug = .1)
plot(x, model = v)
plot(x, model = fit.variogram(x, model = v))

##
## ex03.cmd:
## Inverse distance interpolation on a mask map
##
data(meuse.grid)
gridded(meuse.grid) = ~x+y
x <- krige(zinc ~ 1, meuse, meuse.grid, model = NULL)

library(lattice)

spplot(x[1])

##
## ex04.cmd 
## Local ordinary block kriging at non-gridded locations
##
## the gstat "classic" radius maps into the gstat "S" maxdist argument
##
new.locs <- SpatialPoints(cbind(x = c(181170, 180310, 180205, 178673, 178770, 178270),
	y = c(333250, 332189, 331707, 330066, 330675, 331075)))
krige(zinc ~ 1, meuse, new.locs, 
		model = vgm(1.34e5, "Sph", 800, nug = 2.42e4), 
		block = c(40,40), nmax = 40, nmin = 20, 
		maxdist = 1000)

##
## ex05.cmd 
##
## Local simple point kriging on a mask map
##
v <- vgm(0.581, "Sph", 900, nug = 0.0554)
x <- krige(log(zinc) ~ 1, meuse, meuse.grid, model = v, 
	nmax = 40, nmin = 20, maxdist = 1000, beta = 5.9)
spplot(x[1], main = "log(zinc) simple kriging prediction")
x$se = sqrt(x$var1.var)
spplot(x["se"], main = "log(zinc) simple kriging standard errors")

##
## ex06.cmd 
##
## Unconditional Gaussian simulation on a mask
## (local neigbourhoods, simple kriging)
##
x <- krige(log(zinc) ~ 1, locations = NULL, newdata = meuse.grid, 
	model = v, nmax = 20, beta = c(5.9), nsim = 5, dummy = TRUE)
spplot(x, main = "five unconditional realisations of a correlated Gaussian field")

##
## ex07.cmd 
##
## Gaussian simulation, conditional upon data
## (local neighbourhoods, simple and ordinary kriging)
##
x <- krige(log(zinc) ~ 1, meuse, meuse.grid, 
	model = v, nmax = 20, beta = c(5.9), nsim = 5)
spplot(x, main = "five conditional realisations of a correlated Gaussian field")

##
## ex08.cmd 
##
## Change of support: local ordinary block kriging on a mask
##
x <- krige(log(zinc) ~ 1, meuse, meuse.grid, 
	model = v, nmax = 40, nmin = 20, maxdist = 1000,
	block = c(40,40))
spplot(x[1], main = "ordinary block kriging predictions")
x$se = sqrt(x$var1.var)
spplot(x["se"], main = "ordinary block kriging prediction standard errors")

##
## ex09.cmd 
##
## Obtain map values at data() locations
## (Point-map overlay)
##
# we trick here by using inv.weighted distance interpolation, using the
# single nearest observation. It will not fail on points outside the grid.
# Note that we reversed meuse.grid and meuse to get these results.
x <- krige(part.a ~ 1, meuse.grid, meuse, model = NULL, nmax = 1)
meuse$part.a = x$var1.pred
x <- krige(part.b ~ 1, meuse.grid, meuse, model = NULL, nmax = 1)
meuse$part.b = x$var1.pred

##
## ex10.cmd 
##
## Multiple kriging: prediction on more than one variable
## (ordinary kriging of two variables)
## (note that zinc_map.eas wass obtained through ex09.gst)
##
x <- variogram(dist~1,meuse)
v.dist <- fit.variogram(x, model = vgm(1,"Gau",100))
plot(x, model = v.dist)
g <- gstat(id = "ln.zinc", form = log(zinc) ~ 1,
	data = meuse, nmax = 40, nmin = 20, maxdist = 1000, model = v)
g <- gstat(g, id = "dist", form = dist ~ 1, 
	data = meuse, nmax = 40, nmin = 20, maxdist = 1000,
	model = vgm(.01, "Nug", 0, add.to = v.dist))
# the added nugget variance is necessary to avoid near-singular covariances
x <- predict(g, meuse.grid)
spplot(x["ln.zinc.pred"], main = "log(zinc) ordinary kriging predictions")
x$ln.zinc.se = sqrt(x$ln.zinc.var)
spplot(x["ln.zinc.se"], main = "log(zinc) ordinary kriging prediction standard errors")
spplot(x["dist.pred"], main = "dist ordinary kriging predictions")
x$dist.se = sqrt(x$dist.var)
spplot(x["dist.se"], main = "dist ordinary kriging prediction standard errors")

##
## ex11.cmd 
##
## Multivariable kriging: ordinary local cokriging of two variables
## For examples of fitting an LMC: see demo(cokriging)
##
g <- gstat(id = "ln.zinc", form = log(zinc) ~ 1, 
	data = meuse, nmax = 40, nmin = 20, maxdist = 1000,
	model = vgm(0.581, "Sph", 900, 0.0554))
g <- gstat(g, id = "sq.dist", form = sqrt(dist) ~ 1, 
	data = meuse, nmax = 40, nmin = 20, maxdist = 1000,
	model = vgm(0.0631, "Sph", 900, 0.0001))
g <- gstat(g, id = c("ln.zinc", "sq.dist"), 
	model = vgm(-0.156, "Sph", 900, 1e-5))
# small nugget necessary to let gstat recognize LMC
x <- predict(g, meuse.grid)
spplot(x["ln.zinc.pred"], main = "log(zinc) ordinary cokriging predictions")
x$ln.zinc.se = sqrt(x$ln.zinc.var)
spplot(x["ln.zinc.se"], main = "log(zinc) ordinary cokriging prediction standard errors")
spplot(x["sq.dist.pred"], main = "dist ordinary cokriging predictions")
x$sq.dist.se = sqrt(x$sq.dist.var)
spplot(x["sq.dist.se"], main = "dist ordinary cokriging prediction standard errors")

##
## ex12.cmd 
##
## Stratified ordinary kriging (within-category ordinary kriging)
##

# find out in which part the data are:
meuse$part.a = krige(part.a~1, meuse.grid, meuse, nmax=1)$var1.pred
x1 = krige(log(zinc)~1, meuse[meuse$part.a == 0,], 
	meuse.grid[meuse.grid$part.a == 0,], model = vgm(.548, "Sph", 900, .0654), 
	nmin = 20, nmax = 40, maxdist = 1000)
x2 = krige(log(zinc)~1, meuse[meuse$part.a == 1,], 
	meuse.grid[meuse.grid$part.a == 1,], model = vgm(.716, "Sph", 900), 
	nmin = 20, nmax = 40, maxdist = 1000)
x = rbind(as.data.frame(x1), as.data.frame(x2))
gridded(x) = ~x+y
spplot(x["var1.pred"], main = "stratified kriging predictions")

##
## ex13.cmd 
##
## Local universal kriging, using one continuous variable
###
## the variogram should be that of the residual:
x <- krige(log(zinc) ~ sqrt(dist), meuse, meuse.grid, 
	model = vgm(.149, "Sph", 700, .0674), 
	nmax = 40, nmin = 20, maxdist = 1000)
spplot(x["var1.pred"], main = "universal kriging predictions")
x$var1.se = sqrt(x$var1.var)
spplot(x["var1.se"], main = "universal kriging prediction standard errors")

##
## ex14.cmd 
##
## Universal kriging, using one continuous and
## two binary variables.
##
x <- krige(log(zinc) ~ -1 + sqrt(dist)+ part.a + part.b, 
	meuse, meuse.grid, model = vgm(.149, "Sph", 700, .0674))
spplot(meuse.grid["part.a"], main = "the areas defining part.a (1) and part.b (0)")
spplot(x["var1.pred"], main = "universal kriging predictions")
x$var1.se = sqrt(x$var1.var)
spplot(x["var1.se"], main = "universal kriging prediction standard errors")

##
## ex14a.cmd 
## 
## stratified universal kriging: 
## (again: not implemented)
##

## ex15.cmd 
##
## Local linear model, using one continuous variable
##
x <- krige(log(zinc) ~ sqrt(dist), meuse, meuse.grid, 
	model = NULL, nmax = 40, nmin = 20, maxdist = 1000)
spplot(x["var1.pred"], main = "IID local linear model kriging predictions")
x$var1.se = sqrt(x$var1.var)
spplot(x["var1.se"], main = "IID local linear model prediction standard errors")

##
## ex16.cmd 
##
## Multivariable indicator cosimulation 
## ==>> see demo(cosimulation) for an extended example how to do this
##

##
## ex17.cmd 
##
## global coordinate polynomial trend surfaces
## trend orders 0-3 ==>> better use lm() for this
##
