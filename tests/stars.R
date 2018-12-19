# 0. using sp:

library(sp)
demo(meuse, ask = FALSE)
library(gstat)
v = variogram(log(zinc)~1, meuse)
(v.fit = fit.variogram(v, vgm(1, "Sph", 900, 1)))
k_sp = krige(log(zinc)~1, meuse[-(1:5),], meuse[1:5,], v.fit)
k_sp_grd = krige(log(zinc)~1, meuse, meuse.grid, v.fit)

# 1. using sf:
library(sf)
demo(meuse_sf, ask = FALSE, echo = FALSE)
# reloads meuse as data.frame, so
demo(meuse, ask = FALSE)

library(gstat)
v = variogram(log(zinc)~1, meuse_sf)
(v.fit = fit.variogram(v, vgm(1, "Sph", 900, 1)))
k_sf = krige(log(zinc)~1, meuse_sf[-(1:5),], meuse_sf[1:5,], v.fit)

all.equal(k_sp, as(k_sf, "Spatial"), check.attributes = FALSE)
all.equal(k_sp, as(k_sf, "Spatial"), check.attributes = TRUE)

# 2. using stars for grid:

library(rgdal)
writeGDAL(meuse.grid[,"dist"], "meuse.tif", "GTiff")
library(stars)
(st = setNames(read_stars("meuse.tif"), "dist"))

# compare inputs:
sp = as(st, "Spatial")
fullgrid(meuse.grid) = TRUE
all.equal(sp, meuse.grid["dist"], check.attributes = FALSE)
all.equal(sp, meuse.grid["dist"], check.attributes = TRUE, use.names = FALSE)

# kriging:
st_crs(st) = st_crs(meuse_sf) # GDAL roundtrip messes them up!
k_st = krige(log(zinc)~1, meuse_sf, st, v.fit)
k_st

# handle factors, when going to stars?
k_sp_grd$cls = cut(k_sp_grd$var1.pred, c(0, 5, 6, 7, 8, 9))
st_as_stars(k_sp_grd)
st_as_stars(raster::stack(k_sp_grd)) # check

all.equal(st_redimension(st_as_stars(k_sp_grd)), st_as_stars(raster::stack(k_sp_grd)), check.attributes=FALSE)

