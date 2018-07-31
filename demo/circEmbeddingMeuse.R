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
                        space=vgm(0.85,"Exp", 831, 0.15),
                        time =vgm(0.9,"Exp", 3.25, 0.1),
                        sill=135000)
attr(separableModel,"temporal unit") <- "days"

library(spacetime)
stf <- STF(meuse, Sys.time()-20:0*24*3600)

stf_grid <- STF(geometry(meuse.grid), stf@time)

###################
## unconditional ##
###################

sTime <- Sys.time()
krigedSim <- krigeSTSimTB(newdata = stf_grid, modelList = separableModel, nsim = 100, nLyrs = 100)
Sys.time() - sTime

# plot one simulation along time
stplot(krigedSim[,,"sim1"], main="unconditional siumulation")

# plot one simulation along time as time series
stplot(krigedSim[1:21,,"sim1"], mode="ts", main="unconditional siumulation")
stplot(krigedSim[1:21,,"sim2"], mode="ts", main="unconditional siumulation")
stplot(krigedSim[1:21,,"sim3"], mode="ts", main="unconditional siumulation")
stplot(krigedSim[1:21,,"sim4"], mode="ts", main="unconditional siumulation")
stplot(krigedSim[1:21,,"sim5"], mode="ts", main="unconditional siumulation")
stplot(krigedSim[1:21,,"sim6"], mode="ts", main="unconditional siumulation")

# plot the ten simulations of the first day
spplot(krigedSim[,1], paste0("sim",1:10), as.table=TRUE, main="unconditional siumulation")

#################
## conditional ##
#################

sTime <- Sys.time()
krigedSim <- krigeSTSimTB(formula= zinc ~ 1, data = STFDF(geometry(meuse), stf@time, data.frame(zinc=rep(meuse$zinc, 21))),
                          newdata = stf_grid[1:500,], modelList = separableModel, nsim = 10, nLyrs = 500)
Sys.time() - sTime

# plot one simulation along time
stplot(krigedSim[,1:12], main="conditinal simulation")

# plot one simulation along time as time series
stplot(krigedSim[1:12,,"sim1"], mode="ts", main="conditinal simulation")

# plot the ten simulations of the first day
spplot(krigedSim[,1], paste0("sim",1:10), as.table=TRUE, main="conditinal simulation")