library(sp)
data(meuse)
coordinates(meuse) = c("x", "y")
new.locs <- SpatialPoints(data.frame(
	x = c(181170, 180310, 180205, 178673, 178770, 178270),
	y = c(333250, 332189, 331707, 330066, 330675, 331075)))
library(gstat)
krige(zinc ~ 1,  meuse, new.locs, vgm(1.34e5, "Sph", 800, nug = 2.42e4), 
		block = c(40,40), nmax = 40)

new.locs <- SpatialPoints(data.frame(x = c(181170), y = c(333250)))

disc = c(-15,-5,5,15)
block.irreg <- data.frame(expand.grid(x = disc, y = disc))
block.irreg

# first disable default Gaussian quadrature used for block integration, by
# setting nblockdiscr explicitly:
krige(zinc ~ 1, meuse, new.locs, model = vgm(1.34e5, "Sph", 800, nug = 2.42e4), 
		block = c(40,40), nmax = 40, set = list(nblockdiscr=4))
# now pass the same block discretization as block.irreg:
krige(zinc ~ 1, meuse, new.locs, vgm(1.34e5, "Sph", 800, nug = 2.42e4), 
		block = block.irreg, nmax = 40)
# check weights argument:
block.irreg <- data.frame(expand.grid(x = disc, y = disc), weights = rep(1/16, 16))
krige(zinc ~ 1, meuse, new.locs, vgm(1.34e5, "Sph", 800, nug = 2.42e4), 
		block = block.irreg, nmax = 40)
