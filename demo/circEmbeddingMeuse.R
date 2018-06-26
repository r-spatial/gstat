#############################
## Example: Meuse data set ##
#############################
library(sp)
library(gstat)
data("meuse")
coordinates(meuse) <- ~x+y

data("meuse.grid")
coordinates(meuse.grid) <- ~x+y
gridded(meuse.grid) <- TRUE

# variography
empVgm <- variogram(zinc~1, meuse)
modVgm <- fit.variogram(empVgm, vgm(150000, "Sph", 1000, 25000))

# unconditional simulation
unconSim <- krigeSimCE(zinc~1, newdata = meuse.grid, model = modVgm, n=100)
unconSim@data$zinc.simMean <- apply(unconSim@data[,-c(1:5)], 1, mean)

spplot(unconSim[,6:20], main="15 out of 100 unconditional simulations")
spplot(unconSim, "zinc.simMean", main="mean of 100 unconditional simulations")

# conditional simulation
conSim <- krigeSimCE(zinc~1, meuse, meuse.grid, modVgm, n=100)
conSim@data$zinc.simMean <- apply(conSim@data[,-c(1:5)], 1, mean)

spplot(conSim[,6:20], main="15 out of 100 conditional simulations")

# compare with kriging predictor
simKrige <- krige(zinc~1, meuse, meuse.grid, modVgm)
spplot(simKrige, "var1.pred", main="interpolated zinc concentrations")
spplot(conSim, "zinc.simMean", main="mean of 100 conditional simulations")



################################################
## turning bands simulation in space and time ##
################################################

separableModel <- vgmST("separable",
                        space=vgm(0.9,"Exp", 147, 0.1),
                        time =vgm(0.9,"Exp", 3.5, 0.1),
                        sill=40)


stf <- STF(meuse, Sys.time()-20:0*24*3600)

sTime <- Sys.time()
krigedSim <- krigeSTUncSimTB(stf, separableModel, 100)
Sys.time() - sTime

# 27 secs for 100 simulated ST fields of 155 locations and 21 time steps: 325500 values

# plot one simulation along time
stplot(krigedSim[,1:12])

# plot one simulation along time
stplot(krigedSim[1:12,,"sim1"], mode="ts")

# plot the ten simulations of the first day
spplot(krigedSim[,1], paste0("sim",1:10), as.table=TRUE)