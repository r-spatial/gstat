# $Id: pcb.R,v 1.9 2008-02-01 22:39:44 edzer Exp $
# FIGURE 1:
library(sp)
library(gstat)
library(maptools)
library(lattice)
data(pcb)
coordinates(pcb) = ~x+y
data(ncp.grid)
gridded(ncp.grid) = ~x+y
wsv = readShapePoly(system.file("external/ncp.shp", package="gstat"))
classes = c(.2,1,2,5,10,20)
print(xyplot(y ~ x | as.factor(year), groups = sqrt(PCB138)/3, 
	data = as.data.frame(pcb),
    panel = function(x, y, groups, subscripts, ...) {
				 sp.polygons(wsv, col = grey(.5))
                 panel.xyplot(x, y, cex = groups[subscripts], ...)
			},
    xlab = "x-coordinate", ylab = "y-coordinate",
    key = list(corner = c(0,0), x=0.8, y=0.25, points = list(pch = 1, col = 1,
        cex = sqrt(classes)/3), text = list(as.character(classes))),
	aspect = "iso", as.table = T,
	xlim = c(464000, 739000), ylim = c(5696500, 6131500),
	scales = list(draw = F)))

# FIGURE 2:
pcb$yf = as.factor(pcb$year)
pcb$int = rep(NA, dim(pcb)[1])
c = lm(log(PCB138)~-1+depth+yf,pcb)$coefficients
i = 2
for (f in levels(pcb$yf)) {
	pcb$int[pcb$yf == f] = c[i]
	i = i+1
}
print(xyplot(log(PCB138)~depth | as.factor(year), groups = int, 
	data = as.data.frame(pcb), 
    panel = function(x, y, groups, subscripts, ...) {
                # panel.grid(h=-1, v= 2)
                panel.xyplot(x, y)
				y = c(groups[subscripts][1], groups[subscripts][1] + -0.0641*45)
				llines(c(0, 45), y, lty = 2)
            },
    scales = list(
		y = list(at=log(c(.2, .5, 1, 2, 5, 10, 20)), 
		labels=c(".2",".5","1","2","5","10","20"), alternating = F)
	), xlab = "water depth", ylab = "PCB138", as.table = T)
)

# FIGURE 3:
# ps.options(width=2, height=2)
pcb$res=residuals(lm(log(PCB138)~year+depth, pcb))
v3 = variogram(res ~ year, pcb, dX=.1, bound=c(0,1000,3000,5000,(1:16)*10000))
print(plot(v3, model = vgm(.224,"Exp",17247,.08), plot.numbers = TRUE))

# FIGURE 4:
# next
g.pcb = NULL
merge = list(c("P1986", 2, "P1991", 2), c("P1986", 2, "P1996", 2),
		c("P1986", 2, "P2000", 2))
for (f in levels(pcb$yf)[c(1,4,6,7)])
    g.pcb= gstat(g.pcb, 
		paste("P", as.character(f), sep = ""), 
		log(PCB138)~depth, pcb[pcb$yf == f,],
		merge = merge)
g.pcb = gstat(g.pcb, model = vgm(.224,"Exp",17247,.08), fill.all=T)
v = variogram(g.pcb, cutoff=1e5)
#plot(v, model = fit.lmc(v, g))
#plot(v, model = g,plot.numbers = TRUE)
PCB.cor = matrix(NA, 4,4)
i = 1
for (x in levels(pcb$yf)[c(1,4,6,7)]) {
    j = 1
    for (y in levels(pcb$yf)[c(1,4,6,7)]) {
        A = log(pcb[pcb$yf == x,]$PCB138)
        B.tmp = krige(PCB138~1, pcb[pcb$yf == y,], 
			pcb[pcb$yf == x,], nmax = 1, set=list(debug=0))
        B = log(B.tmp$var1.pred)
        # print(paste(x, y, cor(A,B)))
        PCB.cor[i,j] = cor(A,B)
        j = j + 1
    }
    i = i + 1
}
years=c(1986,1991,1996,2000)
dimnames(PCB.cor)=list(years,years)
PCB.cor.ns = PCB.cor
PCB.cor = 0.5 * (PCB.cor + t(PCB.cor))
i = 1
for (x in levels(pcb$yf)[c(1,4,6,7)]) {
    j = 1
    for (y in levels(pcb$yf)[c(1,4,6,7)]) {
        if (j > i) {
            name = paste(paste("P", x, sep=""), paste("P", y, sep=""),sep = ".")
			print(name)
            g.pcb$model[[name]]["psill"] = 
				g.pcb$model[[name]]["psill"]  * PCB.cor[i,j]
        }
        j = j + 1
    }
    i = i + 1
}
print(plot(v, model = g.pcb, plot.numbers = FALSE))

print(PCB.cor.ns, digits=3)

print(PCB.cor, digits=3)

# FIGURE 5:
pcb.cok = predict(g.pcb, newdata = ncp.grid, debug.level = 0)
levs = c(.1,.2,.5,1,2,5,10,20)

spplot(pcb.cok[c(1,3,5,7)],
	as.table=T, col.regions = bpy.colors(7),
	at = log(levs),
	colorkey = list(at = log(levs), labels = as.character(levs),
		col = bpy.colors(7)), 
	layout = c(4,1)
)

spplot(pcb.cok[c(4,11,6,12,13,8,14,15,16,10)-2],
	skip=c(F,T,T,T,F,F,T,T,F,F,F,T,F,F,F,F),
	as.table=T, layout=c(4,4),
	asp="iso", col.regions = bpy.colors())

X = cbind(rep(1,7), c(1986, 1987, 1989, 1991, 1993, 1996, 2000))
X2=X[c(1,4,6,7),]
lambda = solve(t(X2) %*% X2) %*% t(X2)
dimnames(lambda) = list(NULL, c(1986, 1991, 1996, 2000))
print(lambda[2, ], digits=3)

# FIGURE 7:
pcb.contr = get.contr(pcb.cok, g.pcb, X=lambda[2, ])
# copy coordinates
#pcb.contr$x = pcb.cok$x
#pcb.contr$y = pcb.cok$y

pl1 = spplot(pcb.contr["beta.1"],
	main = "log-PCB138: change estimate",
	col.regions = bpy.colors(100))
pcb.contr$sig = pcb.contr$beta.1 / sqrt(pcb.contr$var.beta.1)
pl2 = spplot(pcb.contr["sig"],
	main = "log-PCB138 change/SE",
	col.regions = bpy.colors(100))
print(pl1, position=c(0,0,0.5,1), more = TRUE)
print(pl2, position=c(0.5,0,1,1), more = FALSE)

cat("source:\n\nEdzer J. Pebesma, Richard N.M. Duin (2005) Spatio-temporal mapping of\nsea floor sediment pollution in the North Sea.  In: Ph. Renard, and\nR. Froidevaux, eds. Proceedings GeoENV 2004 -- Fifth European Conference\non Geostatistics for Environmental Applications; Springer.\n")
