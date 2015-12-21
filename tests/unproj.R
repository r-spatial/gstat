if (require(rgdal) == FALSE)
	q()

# for validity of covariance functions on the sphere, see also:
# DOI 10.1007/s11004-011-9344-7
# http://mypage.iu.edu/~srobeson/Pubs/variogram_sphere_mathgeo_2011.pdf

library(sp)
data(meuse)
coordinates(meuse) = ~x+y
proj4string(meuse) = CRS("+init=epsg:28992")
#meuse.ll = spTransform(meuse, CRS("+proj=longlat +datum=WGS84 +ellps=WGS84"))
meuse.ll = spTransform(meuse, CRS("+proj=longlat +ellps=WGS84"))
meuse.ll[1:10,]
library(gstat)
variogram(log(zinc)~1, meuse.ll)

cloud1 = variogram(log(zinc)~1, meuse, cloud=T, cutoff=6000)
cloud2 = variogram(log(zinc)~1, meuse.ll, cloud=T, cutoff=6)

plot(cloud1$dist/1000, cloud2$dist, xlab="Amersfoort, km", ylab = "Long/lat")
abline(0,1)

if (require(fields)) {
  data(ozone2)
  oz = SpatialPointsDataFrame(ozone2$lon.lat, 
		  data.frame(t(ozone2$y)), 
		  proj4string=CRS("+proj=longlat +datum=WGS84 +ellps=WGS84"))
  variogram(X870731~1,oz[!is.na(oz$X870731),])
  utm16 = CRS("+proj=utm +zone=16 +ellps=WGS84")
  oz.utm = spTransform(oz, utm16)
  variogram(X870731~1,oz.utm[!is.na(oz$X870731),])
}

# Timothy Hilton, r-sig-geo, Sept 14, 2008:

foo <-
structure(list(z = c(-1.95824831109744, -1.9158901643563, 4.22211761150161,
3.23356929459598, 1.12038389231868, 0.34613850821113, 1.12589932643631,
23.517912251617, 3.0519158690268, 3.20261431141517, -2.10947106854739
), lon = c(-125.29228, -82.1556, -98.524722, -99.948333, -104.691741,
-79.420833, -105.100533, -88.291867, -72.171478, -121.556944,
-89.34765), lat = c(49.87217, 48.2167, 55.905833, 56.635833,
53.916264, 39.063333, 48.307883, 40.0061, 42.537756, 44.448889,
46.242017)), .Names = c("z", "lon", "lat"), row.names = c(NA,
-11L), class = "data.frame")

coordinates(foo) <- ~lon+lat

proj4string(foo) <- CRS('+proj=longlat +datum=WGS84 +ellps=WGS84')

vg.foo <- variogram(z~1, foo, cloud=TRUE, cutoff=1e10)

cat('==========\nvariogram:\n')
print(head(vg.foo))

cat('==========\nspDistsN1 Distances:\n')
print(spDistsN1(coordinates(foo), coordinates(foo)[1,], longlat=TRUE))
