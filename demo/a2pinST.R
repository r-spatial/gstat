## area2point in space and time
library(sp)
library(gstat)
data("meuse")
coordinates(meuse) <- ~x+y

data("meuse.grid")
coordinates(meuse.grid) <- ~x+y
gridded(meuse.grid) <- TRUE

meuse.coarse.grid <- SpatialGrid(GridTopology(meuse.grid@grid@cellcentre.offset, c(600,600), c(5, 7)))

separableModel <- vgmST("separable",
                        space=vgm(0.85,"Exp", 831, 0.15),
                        time =vgm(0.9,"Exp", 3.25, 0.1),
                        sill=135000)
attr(separableModel,"temporal unit") <- "days"

library(spacetime)
stf <- STF(meuse[sample(155, 25),], Sys.time()-2:0*24*3600)

stf_grid <- STF(geometry(meuse.coarse.grid), stf@time)

krigedSim <- krigeSTSimTB(newdata = stf_grid, modelList = separableModel, nsim = 1, nLyrs = 100)

# area-to-point kriging:
a2pST = krigeST(sim1 ~ 1, krigedSim, stf, modelList = vgmAreaST, ndiscr = 9, 
            model = separableModel, # point variogram,
            verbose = TRUE)

p1 <- stplot(krigedSim, color.key=F)
p2 <- stplot(a2pST, color.key=F)

print(p1, position=c(0,0.5,1,1), more=TRUE)
print(p2, position=c(0,0,1,0.5))
