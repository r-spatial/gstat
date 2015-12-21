#pdf("wind.pdf")
# PLEASE read the vignette of package spacetime for a more
# clever way to do all this!
library(sp)
library(gstat)
library(rgdal)
library(maptools)
# load wind data, run test:
example(wind)

m = map2SpatialLines(
	map("worldHires", xlim = c(-11,-5.4), ylim = c(51,55.5), plot=F))
proj4string(m) = "+proj=longlat +datum=WGS84 +ellps=WGS84"
m = spTransform(m, CRS("+proj=utm +zone=29 +datum=WGS84 +ellps=WGS84"))

# model temporal autocorrelation
acf(wind[7])
tdiscr = 0:40
lines(tdiscr, exp(- tdiscr/1.5))

# set up data, last year
years = 61
months = 1
jday = c(1,6,11,16,21,26)
sel = wind[wind$year %in% years & 
	wind$month %in% months &
	wind$jday %in% jday,]

#stations = 4:15
stations = 4:15

sels = stack(sel[stations])
sels$t = rep(sel$jday, length(stations))
sels$x = coordinates(wind.loc)[match(sels$ind, wind.loc$Code),1]
sels$y = coordinates(wind.loc)[match(sels$ind, wind.loc$Code),2]
summary(sels)

coordinates(sels) = ~x+y
proj4string(sels) = "+proj=longlat +datum=WGS84 +ellps=WGS84"
sels = spTransform(sels, CRS("+proj=utm +zone=29 +datum=WGS84 +ellps=WGS84"))
grd = makegrid(m, n = 1000)
grd$t = rep(1, nrow(grd))
coordinates(grd) = ~x1+x2
gridded(grd)=TRUE
proj4string(grd) = proj4string(sels)

#sels = as(sels, "data.frame")

# setup grid
covfn = function(x, y = x) { 
	u = spDists(coordinates(x), coordinates(y))
	t = abs(outer(x$t,y$t,"-"))
	0.6 * exp(-u/750000) * exp(-t/1.5)
}
for (i in 1:120) {
	grd$t = rep(i/4, nrow(grd))
	n = paste("t", i/4, sep="")
	grd[[n]] = krige0(sqrt(values)~1, sels, grd, covfn)
}
grd$t = NULL
#grd$pr = out$pred
#library(lattice)
#levelplot(pr~x1+x2|t,grd,col.regions=bpy.colors())
spl = list(list("sp.points", sels,first=F, cex=.5),
	list("sp.lines", m, col='grey'))
spplot(grd, sp.layout = spl, col.regions=bpy.colors())
