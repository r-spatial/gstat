# $Id: block.R,v 1.5 2006-02-10 19:05:02 edzer Exp $
library(sp)
data(meuse)
coordinates(meuse) = ~x+y
data(meuse.grid)
gridded(meuse.grid) = ~x+y

vgm.fit = fit.variogram(variogram(zinc~1, meuse), vgm(1, "Sph", 800, 1))
bl0 = krige(zinc~1, meuse, meuse.grid, model = vgm.fit, block = c(0,0))
bl1 = krige(zinc~1, meuse, meuse.grid, model = vgm.fit, block = c(40,40))
bl2 = krige(zinc~1, meuse, meuse.grid, model = vgm.fit,block = c(100,100))
bl3 = krige(zinc~1, meuse, meuse.grid, model = vgm.fit,block = c(400,400))
bl0$"block=0x0" =     bl0$var1.pred
bl0$"block=40x40" =   bl1$var1.pred
bl0$"block=100x100" = bl2$var1.pred
bl0$"block=400x400" = bl3$var1.pred

plt1 = spplot(bl0, 3:6, layout=c(4,1), col.regions=bpy.colors(), main = "kriging predictions")

bl0$"block=0x0" =     bl0$var1.var
bl0$"block=40x40" =   bl1$var1.var
bl0$"block=100x100" = bl2$var1.var
bl0$"block=400x400" = bl3$var1.var

plt2 = spplot(bl0, 3:6, layout=c(4,1), col.regions=bpy.colors(), main = "kriging standard errors")

print(plt1, split = c(1, 1, 1, 2), more = T)
print(plt2, split = c(1, 2, 1, 2), more = F)
# block krige the full area:

bl = krige(zinc~1, meuse, newdata = SpatialPoints(data.frame(x=0,y=0)),
        model = vgm.fit, block = coordinates(meuse.grid))
bl
# block kriging standard error:
sqrt(bl$var1.var)
# classical statistical standard error of mean:
sqrt(var(meuse$zinc)/155)
