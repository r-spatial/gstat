library(sp)
library(gstat)
library(geoR)
xyz = data.frame(x = c(0,0,1), y = c(0, 1, 1), z = c(1,2,3))
coordinates(xyz)=~x+y
x0 = SpatialPoints(data.frame(x=0,y=.5))
kr1 = krige(z~1,xyz,x0,vgm(1, "Exp", 1))
kr2 = krige.conv(as.geodata(xyz), locations=coordinates(x0), 
	krige=list(cov.model="exponential", cov.par=c(1,1)))
kr1
c(kr2$predict, kr2$krige.var)
kr1[[1]] - kr2$predict
kr1[[2]] - kr2$krige.var
