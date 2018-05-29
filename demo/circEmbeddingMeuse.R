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