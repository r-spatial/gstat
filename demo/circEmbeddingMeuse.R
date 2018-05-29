#############################
## Example: Meuse data set ##
#############################

# unconditional simulation
unconSim <- krigeSimCE(zinc~1, newdata = meuse.grid, model = modVgm, n=100)
unconSim@data$zinc.simMean <- apply(unconSim@data[,-c(1:5)], 1, mean)

spplot(unconSim[,6:20])
spplot(unconSim, "zinc.simMean")

# conditional simulation
conSim <- krigeSimCE(zinc~1, meuse, meuse.grid, modVgm, n=100)
conSim@data$zinc.simMean <- apply(conSim@data[,-c(1:5)], 1, mean)

spplot(conSim[,6:20])

# compare with kriging predictor
simKrige <- krige(zinc~1, meuse, meuse.grid, modVgm)
spplot(simKrige, "var1.pred")
spplot(simRes, "zinc.simMean")