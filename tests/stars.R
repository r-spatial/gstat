Sys.setenv(TZ = "UTC")

Sys.unsetenv("KMP_DEVICE_THREAD_LIMIT")
Sys.unsetenv("KMP_ALL_THREADS")
Sys.unsetenv("KMP_TEAMS_THREAD_LIMIT")
Sys.unsetenv("OMP_THREAD_LIMIT")

# 0. using sp:

suppressPackageStartupMessages(library(sp))
demo(meuse, ask = FALSE)
suppressPackageStartupMessages(library(gstat))
v = variogram(log(zinc)~1, meuse)
(v.fit = fit.variogram(v, vgm(1, "Sph", 900, 1)))
k_sp = krige(log(zinc)~1, meuse[-(1:5),], meuse[1:5,], v.fit)
k_sp_grd = krige(log(zinc)~1, meuse, meuse.grid, v.fit)

# 1. using sf:
suppressPackageStartupMessages(library(sf))
demo(meuse_sf, ask = FALSE, echo = FALSE)
# reloads meuse as data.frame, so
demo(meuse, ask = FALSE)

v = variogram(log(zinc)~1, meuse_sf)
(v.fit = fit.variogram(v, vgm(1, "Sph", 900, 1)))
k_sf = krige(log(zinc)~1, meuse_sf[-(1:5),], meuse_sf[1:5,], v.fit)

all.equal(k_sp, as(k_sf, "Spatial"), check.attributes = FALSE)
all.equal(k_sp, as(k_sf, "Spatial"), check.attributes = TRUE)

# 2. using stars for grid:

suppressPackageStartupMessages(library(stars))
st = st_as_stars(meuse.grid)
st_crs(st)

# compare inputs:
sp = as(st, "Spatial")
fullgrid(meuse.grid) = TRUE
all.equal(sp, meuse.grid["dist"], check.attributes = FALSE)
all.equal(sp, meuse.grid["dist"], check.attributes = TRUE, use.names = FALSE)

# kriging:
st_crs(st) = st_crs(meuse_sf) = NA # GDAL roundtrip messes them up!
k_st = if (Sys.getenv("USER") == "travis") {
	try(krige(log(zinc)~1, meuse_sf, st, v.fit))
} else {
	krige(log(zinc)~1, meuse_sf, st, v.fit)
}
k_st

# handle factors, when going to stars?
k_sp_grd$cls = cut(k_sp_grd$var1.pred, c(0, 5, 6, 7, 8, 9))
st_as_stars(k_sp_grd)
if (require(raster, quietly = TRUE)) {
 print(st_as_stars(raster::stack(k_sp_grd))) # check
 print(all.equal(st_redimension(st_as_stars(k_sp_grd)), st_as_stars(raster::stack(k_sp_grd)), check.attributes=FALSE))
}

suppressPackageStartupMessages(library(spacetime))

tm = as.POSIXct("2019-02-25 15:37:24 CET")
n = 4
s = stars:::st_stars(list(foo = array(1:(n^3), rep(n,3))),
stars:::create_dimensions(list(
  x = stars:::create_dimension(from = 1, to = n, offset = 10, delta = 0.5),
  y = stars:::create_dimension(from = 1, to = n, offset = 0, delta = -0.7),
  time = stars:::create_dimension(values = tm + 1:n)),
  raster = stars:::get_raster(dimensions = c("x", "y")))
  )
s

as.data.frame(s)
plot(s, col = sf.colors(), axes = TRUE)
(s.stfdf = as(s, "STFDF"))
stplot(s.stfdf, scales = list(draw = TRUE))

(s2 = st_as_stars(s.stfdf))
plot(s2, col = sf.colors(), axes = TRUE)
all.equal(s, s2, check.attributes = FALSE)

# multiple simulations:
data(meuse, package = "sp")
data(meuse.grid, package = "sp")
coordinates(meuse.grid) <- ~x+y
gridded(meuse.grid) <- TRUE
meuse.grid = st_as_stars(meuse.grid)
meuse_sf = st_as_sf(meuse, coords = c("x", "y"))
g = gstat(NULL, "zinc", zinc~1, meuse_sf, model = vgm(1, "Exp", 300), nmax = 10)
g = gstat(g, "lead", lead~1, meuse_sf, model = vgm(1, "Exp", 300), nmax = 10, fill.cross = TRUE)
set.seed(123)
## IGNORE_RDIFF_BEGIN
(p = predict(g, meuse.grid, nsim = 5))
## IGNORE_RDIFF_END
