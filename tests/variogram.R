library(sp)
library(gstat)
data(meuse)
variogram(log(zinc)~1, ~x+y, meuse)

coordinates(meuse) <- ~ x + y                               
variogram(log(zinc)~1, meuse)

ind=seq(1,155,2)
var1= meuse[ind,]
var2= meuse[-ind,]
g <- gstat(NULL, id = "lead", form = lead ~ 1, data=var1)
g <- gstat(g, id = "zinc", form = zinc ~ 1, data=var2)
v.cross <- variogram(g)
plot(v.cross)

