# test -- load data:
library(sp)
data(meuse)
coordinates(meuse) = ~x+y
data(meuse.grid)
gridded(meuse.grid) = ~x+y

library(gstat)
# test -- idw
meuse.grid$idw <- idw0(zinc~1, meuse, meuse.grid)[,1]
x <- idw(zinc~1, meuse, meuse.grid)
all.equal(x$var1.pred, meuse.grid$idw)
spplot(meuse.grid["idw"],col.regions=bpy.colors())
v = vgm(1, "Exp", 500)
# test sk:
x0 <- krige0(zinc~1, meuse, meuse.grid, v, beta = 500, computeVar = TRUE)
rownames(x0$pred)=NULL
x <- krige(zinc~1, meuse, meuse.grid, v, beta = 500)
all.equal(x$var1.pred, x0$pred[,1])
all.equal(x$var1.var, x0$var)
# test ok:
x0 <- krige0(zinc~1, meuse, meuse.grid, v, computeVar = TRUE)
rownames(x0$pred)=NULL
names(x0$var)=NULL
x <- krige(zinc~1, meuse, meuse.grid, v)
all.equal(x$var1.pred, x0$pred[,1])
all.equal(x$var1.var, x0$var)
# test uk:
x0 <- krige0(zinc~sqrt(dist), meuse, meuse.grid, v, computeVar = TRUE)
rownames(x0$pred)=NULL
names(x0$var)=NULL
x <- krige(zinc~sqrt(dist), meuse, meuse.grid, v)
all.equal(x$var1.pred, x0$pred[,1])
all.equal(x$var1.var, x0$var)