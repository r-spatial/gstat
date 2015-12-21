# $Id: cokriging.R,v 1.4 2006-02-10 19:05:02 edzer Exp $
library(sp)
data(meuse)
coordinates(meuse) = ~x+y
data(meuse.grid)
gridded(meuse.grid) = ~x+y

# cokriging of the four heavy metal variables
meuse.g <- gstat(id="zn", formula=log(zinc)~1, data=meuse, nmax = 10)
meuse.g <- gstat(meuse.g, "cu", log(copper)~1, meuse, nmax = 10)
meuse.g <- gstat(meuse.g, "cd", log(cadmium)~1, meuse, nmax = 10)
meuse.g <- gstat(meuse.g, "pb", log(lead)~1, meuse, nmax = 10)
meuse.g <- gstat(meuse.g, model=vgm(1, "Sph", 900, 1), fill.all=T)
x <- variogram(meuse.g, cutoff=1000)
meuse.fit = fit.lmc(x, meuse.g)
plot(x, model = meuse.fit)
z <- predict(meuse.fit, newdata = meuse.grid)

library(lattice)

pl1 <- spplot(z["zn.pred"], main="log-zinc predictions")
pl2 <- spplot(z["cu.pred"], main="log-copper predictions")
pl3 <- spplot(z["cd.pred"], main="log-cadmium predictions")
pl4 <- spplot(z["pb.pred"], main="log-lead predictions")
print(pl1, split = c(1,1,2,2), more=TRUE)
print(pl2, split = c(1,2,2,2), more=TRUE)
print(pl3, split = c(2,1,2,2), more=TRUE)
print(pl4, split = c(2,2,2,2))
z$zn.se = sqrt(z$zn.var)
z$cu.se = sqrt(z$cu.var)
z$pb.se = sqrt(z$pb.var)
z$cd.se = sqrt(z$cd.var)
pl1 <- spplot(z["zn.se"], main="log-zinc std.err.")
pl2 <- spplot(z["cu.se"], main="log-copper std.err.")
pl3 <- spplot(z["cd.se"], main="log-cadmium std.err.")
pl4 <- spplot(z["pb.se"], main="log-lead st.err.")
print(pl1, split = c(1,1,2,2), more=TRUE)
print(pl2, split = c(1,2,2,2), more=TRUE)
print(pl3, split = c(2,1,2,2), more=TRUE)
print(pl4, split = c(2,2,2,2))

rm(meuse.g, x, meuse.fit, z)

# indicator cokriging for the 9 percentiles of zinc:
q <- quantile(meuse$zinc, seq(.1,.9,.1))
meuse.i <- gstat(id = "zn1", formula = I(zinc < q[1])~1, 
	data = meuse, nmax = 7, beta = .1, set = list(order = 4, zero = 1e-5))
meuse.i <- gstat(meuse.i, "zn2", I(zinc < q[2])~1, meuse, nmax = 7, beta=.2)
meuse.i <- gstat(meuse.i, "zn3", I(zinc < q[3])~1, meuse, nmax = 7, beta=.3)
meuse.i <- gstat(meuse.i, "zn4", I(zinc < q[4])~1, meuse, nmax = 7, beta=.4)
meuse.i <- gstat(meuse.i, "zn5", I(zinc < q[5])~1, meuse, nmax = 7, beta=.5)
meuse.i <- gstat(meuse.i, "zn6", I(zinc < q[6])~1, meuse, nmax = 7, beta=.6)
meuse.i <- gstat(meuse.i, "zn7", I(zinc < q[7])~1, meuse, nmax = 7, beta=.7)
meuse.i <- gstat(meuse.i, "zn8", I(zinc < q[8])~1, meuse, nmax = 7, beta=.8)
meuse.i <- gstat(meuse.i, "zn9", I(zinc < q[9])~1, meuse, nmax = 7, beta=.9)
meuse.i <- gstat(meuse.i, model=vgm(1, "Sph", 900, 1), fill.all=T)
x <- variogram(meuse.i, cutoff=1000)
meuse.fit = fit.lmc(x, meuse.i)
plot(x, model = meuse.fit)
z <- predict(meuse.fit, newdata = meuse.grid)
spplot(z, c(3,5,7,9,11,13,15,17,19), 
	names.attr = paste("est.Pr(Zn < ", q, ")", sep = ""))
