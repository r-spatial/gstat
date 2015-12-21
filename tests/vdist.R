library(sp)
library(gstat)

data(meuse)
coordinates(meuse) = ~x+y
data(meuse.grid)
gridded(meuse.grid) = ~x+y

mg = meuse.grid
gridded(mg) = FALSE
mg= mg[1500,]
krige(log(zinc)~1,meuse,mg,vgm(1, "Exp", 300, anis=c(0,0.01)),
	vdist=FALSE, maxdist=1000,nmax=10)
krige(log(zinc)~1,meuse,mg,vgm(1, "Exp", 300, anis=c(0,0.01)),
	vdist=TRUE, maxdist=1000,nmax=10)
