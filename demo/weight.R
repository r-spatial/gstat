kriging.weights = function(x, formula, newdata, model) {
	weighti = function(x, i, formula,...) {
		ret =rep(0,nrow(x))
		ret[i]=1
		x[[1]]=ret
		krige(formula = formula,locations = x,...)
	}
	ret = sapply(1:nrow(x), weighti, x=x, newdata=newdata[1,], model=model,formula=formula)
	ret = t(sapply(ret, as.data.frame))
	unlist(ret[,3])
}
# example, at first cell of meuse.grid:
require(sp)
require(gstat)
data(meuse)
data(meuse.grid)
coordinates(meuse) = ~x+y
coordinates(meuse.grid) = ~x+y
meuse$wts = kriging.weights(meuse["zinc"], zinc~1, meuse.grid[1,], vgm(1, "Exp", 300))
summary(meuse$wts)
spplot(meuse["wts"], col.regions=bpy.colors(), cuts=(0:10)/20)
