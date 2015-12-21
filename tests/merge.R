options(digits=6)
# illustrates the use of merge, for merging parameters accross variables:
# Z1=m+e1(s)
# Z2=m+e2(s)
# Z1 and Z2 each have a different variogram, but share the parameter m
# see documentation of gstat() function
library(gstat)
d1 = data.frame(x=c(0,2),y=c(0,0),z=c(0,1))
d2 = data.frame(x=c(0,2),y=c(2,2),z=c(4,5))
g = gstat(NULL,"d1", z~1,~x+y,d1,model=vgm(1, "Exp", 1))
g = gstat(g,"d2", z~1,~x+y,d2,model=vgm(1, "Exp", 1), merge=c("d1","d2"))
g = gstat(g, c("d1", "d2"), model = vgm(0.5, "Exp", 1))
predict(g, data.frame(x=1,y=1), debug = 32)

# Z1 and Z2 share a regression slope:
g = gstat(NULL,"d1", z~x,~x+y,d1,model=vgm(1, "Exp", 1))
g = gstat(g,"d2", z~x,~x+y,d2,model=vgm(1, "Exp", 1), 
	merge=list(c("d1",2,"d2",2)))
g = gstat(g, c("d1", "d2"), model = vgm(0.5, "Exp", 1))
predict(g, data.frame(x=1,y=1), debug = 32)
