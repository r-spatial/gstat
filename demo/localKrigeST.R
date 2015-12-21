## FNN local prediction
########################
library(sp)
library(spacetime)
library(gstat)
library(lattice)

# create n space-time points over [0,1] x [0,1] x [Now, Now+some days]
t0 = Sys.time() # now
n = 1000
set.seed(13131) # fix outcomes
x = runif(n)
y = runif(n)
t = t0 + 1e6 * runif(n)
z = rnorm(n)
stidf = STIDF(SpatialPoints(cbind(x,y)), sort(t), data.frame(z=z))

stplot(stidf, number=21, main="random spatio-temporal noise")

# create a regular 20 x 20 x 10 grid of prediction locations:
grd = as(SpatialGrid(GridTopology(c(0.025,0.025), c(.05, .05), c(20,20))), "SpatialPixels")
tgrd = seq(min(t)+10000, max(t)-10000, length.out = 10)
stf = STF(grd, tgrd)

# define a variogram model
sumMetricModel <- vgmST("sumMetric",
                        space=vgm(1/6, "Sph", 0.25, 1/60),
                        time =vgm(2/6, "Exp",  1e5, 1/60),
                        joint=vgm(0.4, "Exp", 0.3, 0.1),
                        stAni=1/1e6)
attr(sumMetricModel, "temporal unit") <- "secs"

dg <- data.frame(spacelag=rep(c(0.001,1:10)/10,6), 
                 timelag=rep(0:5*50e3, each=11))
wireframe(model~spacelag+timelag,
          variogramSurface(sumMetricModel, dist_grid = dg),
          scales=list(arrows=F),
          drape=T, col.regions=bpy.colors(),
          zlim=c(0,1.2),
          main="imposed sum-metric model")

locKrig <- krigeST(z~1, stidf, stf, sumMetricModel, nmax=50, computeVar = T)
stplot(locKrig[,,"var1.pred"], col.regions=bpy.colors(), scales=list(draw=T))
stplot(locKrig[,,"var1.var"], col.regions=bpy.colors(), scales=list(draw=T))