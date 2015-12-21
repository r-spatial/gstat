# $Id: gstat3D.R,v 1.5 2007-02-23 13:34:07 edzer Exp $
# simple demo of 3D interpolation of 50 points with random normal values,
# randomly located in the unit cube
library(sp)
library(gstat)
n <- 50

data3D <- data.frame(x = runif(n), y = runif(n), z = runif(n), v = rnorm(n))
coordinates(data3D) = ~x+y+z

range1D <- seq(from = 0, to = 1, length = 20)
grid3D <- expand.grid(x = range1D, y = range1D, z = range1D)
gridded(grid3D) = ~x+y+z

res3D <- krige(formula = v ~ 1, data3D, grid3D, model = vgm(1, "Exp", .2))

library(lattice)

levelplot(var1.pred ~ x + y | z, as.data.frame(res3D))
rm(n, data3D, range1D, grid3D, res3D)
