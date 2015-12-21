library(gstat)
d=expand.grid(x=c(-.5,.5), y=c(-.5,.5))
d$z=1:4
vv=vgm(model = "Tab",  covtable = 
	variogramLine(vgm(1, "Sph", 1), 1, n=1e4,min = 0, covariance = TRUE))
vv
krige(z~1,~x+y,d,data.frame(x=0,y=0),vgm(1, "Sph", 1))
krige(z~1,~x+y,d,data.frame(x=0,y=0),vv)
krige(z~1,~x+y,d[1:2,],data.frame(x=0,y=0),vgm(1, "Sph", 1))
krige(z~1,~x+y,d[1:2,],data.frame(x=0,y=0),vv)
