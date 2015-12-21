# how to get the BLUE trend coefficients out of a predict.gstat call?
# prepare data
library(sp)
library(gstat)
data(meuse)
coordinates(meuse) = ~x+y
data(meuse.grid)
gridded(meuse.grid) = ~x+y

# create a manual, non-automatic intercept
meuse$Int = rep(1, length(meuse$zinc))
meuse.grid$Int = rep(1, length(meuse.grid$dist))

# create a gstat object without an automatic intercept:
g = gstat(formula = log(zinc)~-1+Int+sqrt(dist), data=meuse, model = vgm(1, "Exp", 300))
newdat = meuse.grid[1:2,c("Int", "dist")]
gridded(newdat) = FALSE
newdat$dist = c(0,1)
newdat$Int = c(1,0)
# (in the more general case of n predictors, make sure that the matrix with
# predictors, or design matrix, equals the identity matrix.)
newdat
out = as.data.frame(predict(g, newdat, BLUE = TRUE))[,3:4]
rownames(out) = c("Intercept", "slope")
colnames(out) = c("BLUE", "Variance")
out
