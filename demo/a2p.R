# Rprof()
# import NC SIDS data:
library(sf)
demo(nc, ask = FALSE)

# reproject to UTM17, so we can use Euclidian distances:
nc = st_transform(nc, st_crs("+proj=utm +zone=17 +datum=WGS84 +ellps=WGS84"))

# create a target (newdata) grid, and plot:
grd = st_sample(nc, type = "regular", size = 1000)
plot(st_geometry(nc))
plot(grd, pch = 3, add = TRUE)

library(gstat)

# area-to-point kriging:
kr = krige0(SID74 ~ 1, nc, grd, vgmArea, ndiscr = 9, 
	vgm = vgm(1, "Exp", 1e5, 0), # point variogram,
	verbose = TRUE)
out = st_as_stars(grd)

plot(nc["SID74"], main = "areas", key.pos = 0)
plot(out[1], reset = FALSE)
plot(st_geometry(nc), col = NA, border = 'black', add = TRUE)
