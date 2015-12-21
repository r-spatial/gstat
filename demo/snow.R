library(sp)
library(gstat)
data(meuse.grid)
gridded(meuse.grid) = ~x+y
data(meuse)
coordinates(meuse) = ~x+y

ncell = 1000000
# sample 1000000 points over meuse.grid:
newd = spsample(meuse.grid, ncell, type="regular")
ncell = dim(coordinates(newd))[1]
v = vgm(0.6, "Sph", 900, 0.05)

library(parallel)

nclus = detectCores()
clus <- c(rep("localhost", nclus))

# set up cluster and data
cl <- makeCluster(clus, type = "SOCK")
clusterEvalQ(cl, library(gstat))
clusterExport(cl, list("meuse", "meuse.grid", "v"))

# split prediction locations:

# either at random:
splt = sample(1:nclus, nrow(coordinates(newd)), replace = TRUE)

# or regular:
splt = rep(1:nclus, each = ceiling(ncell/nclus), length.out = ncell)
newdlst = lapply(as.list(1:nclus), function(w) newd[splt == w,])

# no cluster:
system.time(out.noclus <- krige(log(zinc)~1, meuse, newd, v))

# with cluster:
system.time(out.clus <- do.call("rbind", parLapply(cl, newdlst, function(lst) 
	krige(log(zinc)~1, meuse, lst, v)
)))
stopCluster(cl)
all.equal(out.clus, out.noclus)
gridded(out.clus) = TRUE
image(out.clus)
