# $Id: sic2004.R,v 1.4 2006-11-22 12:54:16 edzer Exp $
# compared to the original submission, two errors
# were found; one was corrected (dayx~x -> dayx~1), the other not.

# read csv files:
#sic.pred = read.csv("SIC2004_out.csv", header=F)
#names(sic.pred) = c("record", "x", "y")
#input = read.csv("SIC2004_input.csv", header=F)
#names(input) = c("record", "x", "y", "dayx")
#joker = read.csv("SIC2004_joker.csv", header=F)
#names(joker) = c("record", "x", "y", "dayx")

library(sp)
library(gstat)
library(lattice)
data(sic2004) # load directly from R data base
input = sic.test[,1:4]
names(input) = c("record", "x", "y", "dayx")
joker = sic.test[,c(1,2,3,5)]
names(joker) = c("record", "x", "y", "dayx")

do.sic = function(sic.input, sic.output) {
    require(gstat)
    # calculate omnidirectional sample variogram, cutoff 500000, intervals 20000
    sic.input = sic.input[sic.input$record != 838, ]
    v.sample = variogram(dayx~1, ~x+y, input, width=2e4, cutoff=5e5)
    #                                  ^^^^^ the error: this should be sic.input
	#                         ^ this error has been corrected
    # initial spherical model for fit: nug=100, p.sill=400,range=4.5e5
    initial.model = vgm(400, "Sph", 450000, 100)
    # fit the nugget, sill and range to the sample variogram:
    v.fitted = fit.variogram(v.sample, initial.model, fit.method = 2)
    # function returns the output from the call to krige():
    krige(dayx~1, ~x+y, sic.input, sic.output, model = v.fitted, nmax = 125)
}

# do the spatial interpolations:
output = do.sic(input, sic.pred)
joker.output = do.sic(joker, sic.pred)

# write csv files:
#write.table(output, "sic2004_output.csv", 
#    sep=",", col.names=FALSE, row.names=FALSE)
#write.table(joker.output, "sic2004_joker_output.csv",
#    sep=",", col.names=FALSE, row.names=FALSE)

# calculate output:

# read csv files:
sic.stats = function(obs, pred) {
	print("observed:")
	x = obs
	print(c(min=min(x), max=max(x), mean=mean(x), median=median(x), stddev=sqrt(var(x))))
	print("predicted:")
	x = pred
	print(c(min=min(x), max=max(x), mean=mean(x), median=median(x), stddev=sqrt(var(x))))
	print("error:")
	x = pred - obs
	print(c(mae=mean(abs(x)), me=mean(x), corr=cor(pred,obs), rmse=sqrt(mean(x^2))))

}

options(digits=4)
print("input data first")
sic.stats(sic.test[,4], output$var1.pred)
print("joker with error variogram")
sic.stats(sic.test[,4], joker.output$var1.pred)

##################

do.sic = function(sic.input, sic.output) {
    require(gstat)
    v.sample = variogram(dayx~1, ~x+y, input, width=2e4, cutoff=5e5)
    initial.model = vgm(400, "Sph", 450000, 100)
    v.fitted = fit.variogram(v.sample, initial.model, fit.method = 2)
    krige(dayx~1, ~x+y, sic.input, sic.output, model = v.fitted, nmax = 125)
}

# do the spatial interpolations on a grid:
grid1 = do.sic(input, sic.grid)
grid2 = do.sic(joker, sic.grid)

grid1$v2 = grid2$var1.pred
grid1$v2.var = grid2$var1.var

grid1$var1.pred[grid1$var1.pred>200] = 205 # mask larger values!
grid1$v2[grid1$v2>200] = 205               # mask larger values!

greys = c(255, 247, 240, 228, 217, 203, 189, 169, 150, 132, 115, 99, 82, 60, 37, 0)
panel.sic = function(...){
	panel.levelplot(...)
	lpoints(sic.train$x, sic.train$y, pch=22, cex=.5, col=1)
	lpoints(sic.pred$x, sic.pred$y, pch="+", cex=1, col=1)
}
cp1 = contourplot(z~x+y|name, map.to.lev(grid1,z=c(3,5),
    ns=c("data set 1","data set 2")), asp="iso",
    col.regions=rgb(greys/255,greys/255,greys/255),
    at=50+(0:16)*10,region=TRUE,labels=FALSE,
    xlab = "", ylab = "",scales=list(draw=F),panel=panel.sic)
cp2 = contourplot(sqrt(z)~x+y|name, map.to.lev(grid1,z=4,
    ns="standard error"), asp="iso",scales=list(draw=F),
    col.regions=rgb(greys/255,greys/255,greys/255),
    region=TRUE,labels=FALSE, xlab = "", ylab = "", panel=panel.sic)
print(cp1, c(0,0,0.625,1), more=T)
print(cp2, c(0.625,0,1,1), more=F)

wireframe(var1.pred~x+y,grid2,asp=c(diff(range(grid2$y))/diff(range(grid2$x)),0.5), 
    xlab="x",ylab="y",zlab="z",drape=T,
    col.regions=gray(sqrt(seq(from=1.0, to=0.0, length=100))))
