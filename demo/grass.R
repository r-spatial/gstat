# $Id: grass.R,v 1.4 2006-02-10 19:05:02 edzer Exp $
# this demo assumes quite a lot:
#  a. it assumes GRASS gis is running
#  b. it assumes that the meuse data zinc variable is available as a site list
#  c. it assumes that mask_map is present, and contains the mask map values
#     (i.e., the study area)

library(sp)
library(GRASS)           # load R GRASS interface

G = gmeta()              # retrieves active data base locations and topology
d = sites.get(G, "zinc") # retrieve zinc observations 
plot(d$east, d$north, asp=1)
names(d)[4] = "zinc"     # rename attribute
mask = rast.get

hist(d$zinc)
hist(log(d$zinc))

mask = rast.get(G, "mask_map")
plot(G, mask$mask.map)
points(d$east,d$north, pch="+")

library(gstat)           # load gstat library

bubble(d, zcol = "zinc", col=c(4,5), maxsize=2)

# explain S formulae: ~
v = variogram(log(zinc)~1, ~east+north, d)
plot(v)

v.mod = vgm(.6, "Sph", 900, .1)
plot(v, model = v.mod)

v.fit = fit.variogram(v, v.mod)
plot(v, model = v.fit)

zinc.g = gstat(NULL, "lzinc", log(zinc)~1, ~east+north, d, model = v.fit)
new.data = data.frame(east = east(G), north = north(G))
new.data[is.na(mask$mask.map), ] = c(NA,NA)

zinc.kr = predict(zinc.g, new.data)
image(zinc.kr)

library(lattice)

levelplot(lzinc.pred~east+north, zinc.kr, asp=1.34, col.regions=bpy.colors(100))

# push prediction and variances grids back into GRASS data base:
rast.put(G, "lzinc.pred", zinc.kr$lzinc.pred)
rast.put(G, "lzinc.var",  zinc.kr$lzinc.var)

# push cross validation residuals back to GRASS data base:
xv = krige.cv(log(zinc)~1, ~east+north, d, v.fit, nmax = 40, verb=F)
sites.put2(G, data = xv, dims = c("east", "north", "residual", "zscore"), 
	lname = "lzinc.xv")
