library(sp)
library(gstat)
data(meuse.grid)
gridded(meuse.grid) = ~x+y
data(meuse)
coordinates(meuse) = ~x+y

# choose arbitrary line over the grid:
image(meuse.grid["dist"],axes=T)
pp = rbind(c(180000,331000),c(180000,332000),c(181000,333500))
Sl = SpatialLines(list(Lines(list(Line(pp)), "a")))
plot(Sl,add=T,col='green')

# use the default spsample arguments of predict.gstat:
pts=spsample(Sl,n=500,'regular',offset=c(.5,.5))
plot(pts, pch=3, cex=.2, add=T)

v = vgm(.6, "Sph", 900, .06)
out1 = krige(log(zinc)~1, meuse, Sl, v)
out1

points(180333,332167,pch=3,cex=2)

# use the same line as block discretization, and predict for (0,0)
# (because the block discretizing points are not centered)
out2 = krige(log(zinc)~1, meuse, SpatialPoints(matrix(0,1,2)), v, block=coordinates(pts))
out2

compare.krigingLines = function(formula, data, newdata, model) {
	out1 = krige(formula, data, newdata, model)
	pts = spsample(newdata, n=500, 'regular', offset=.5)
	out2 = krige(formula, data, SpatialPoints(matrix(0,1,2)), model, block = coordinates(pts))
	print("difference:")
	as.data.frame(out1)[3:4]- as.data.frame(out2)[3:4]
}

compare.krigingLines(log(zinc)~1, meuse, Sl, v)

# one line, consisting of two line segments:
pp2 = rbind(c(181000,333500),c(181000,332500))
Sl2 = SpatialLines(list(Lines(list(Line(pp),Line(pp2)), "a")))
krige(log(zinc)~1, meuse, Sl2, v)
compare.krigingLines(log(zinc)~1, meuse, Sl2, v)

# two seperate line segments:
Sl3 = SpatialLines(list(Lines(list(Line(pp)), "a"),Lines(list(Line(pp2)),"b")))
krige(log(zinc)~1, meuse, Sl3, v)
