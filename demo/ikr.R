library(sp)
library(gstat)
data(meuse)
data(meuse.grid)
coordinates(meuse)=~x+y
gridded(meuse.grid)=~x+y
v = variogram(I(zinc < 500)~1,meuse)
plot(v)
vm = fit.variogram(v, vgm(1, "Sph", 300, 1))
plot(v,vm)
vm
# possibly adjust sum of sill to be max. 0.25?
ik = krige(I(zinc>500)~1, meuse, meuse.grid, vm)
spplot(ik[1],col.regions=bpy.colors())
summary(ik[[1]])
# adjust values outside [0,1] to nearest limit:
ik[[1]][ik[[1]]<0] = 0
ik[[1]][ik[[1]]>1] = 1
summary(ik[[1]])
spplot(ik[1],col.regions=bpy.colors())
