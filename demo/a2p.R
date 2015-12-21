Rprof()
# import NC SIDS data:
library(sp)
library(maptools)
fname = system.file("shapes/sids.shp", package="maptools")[1]
nc = readShapePoly(fname, proj4string = 
	CRS("+proj=longlat +datum=NAD27 +ellps=clrk66"))

# reproject to UTM17, so we can use Euclidian distances:
library(rgdal)
nc = spTransform(nc, CRS("+proj=utm +zone=17 +datum=WGS84 +ellps=WGS84"))

# create a target (newdata) grid, and plot:
grd = spsample(nc, "regular", n = 1000)
class(grd)
plot(nc, axes = TRUE)
points(grd, pch = 3)

library(gstat)

# area-to-point kriging:
kr = krige0(SID74 ~ 1, nc, grd, vgmArea, ndiscr = 9, 
	vgm = vgm(1, "Exp", 1e5, 0), # point variogram,
	verbose = TRUE)
out = SpatialPixelsDataFrame(grd, data.frame(pred = kr))

pl0 = spplot(nc["SID74"], main = "areas")
pl1 = spplot(out, sp.layout = list("sp.polygons", nc, first=F,col='grey'), 
    main = "points on a grid")
print(pl0, split = c(1,1,1,2), more = TRUE)
print(pl1, split = c(1,2,1,2), more = FALSE)

