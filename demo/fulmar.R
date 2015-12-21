# $Id: fulmar.R,v 1.3 2006-02-10 19:05:02 edzer Exp $
library(sp)
library(gstat)
data(fulmar)
data(ncp.grid)

glm98 <- glm(formula = fulmar ~ depth + coast, family = quasipoisson,
    data = fulmar[fulmar$year == 1998, ])
glm99 <- glm(formula = fulmar ~ depth + coast, family = quasipoisson,
    data = fulmar[fulmar$year == 1999, ])
fulmar98 = data.frame(fulmar[fulmar$year == 1998,], 
	pr98 = predict(glm98, type = "response"))
fulmar99 <- data.frame(fulmar[fulmar$year == 1999,], 
	pr99 = predict(glm99, type = "response"))
pr98.grd <- predict(glm98, newdata = ncp.grid, type = "response", se.fit=TRUE)
pr99.grd <- predict(glm99, newdata = ncp.grid, type = "response", se.fit=TRUE)
pr <- data.frame(ncp.grid, pr98=pr98.grd$fit, pr99=pr99.grd$fit,
	se98 = pr98.grd$se.fit, se99 = pr99.grd$se.fit)

# B.3 create gstat object
g <- gstat(id = "fulmar98", formula = fulmar~pr98, locations = ~x+y, 
	data = fulmar98, model = vgm(1.89629, "Exp", 50000, 0.852478), 
	beta = c(0,1), variance = "mu")
g <- gstat(g, id = "fulmar99", formula = fulmar~pr99, locations = ~x+y, 
	data = fulmar99, model = vgm(2.52259, "Exp", 50000, 1.76474), 
	beta = c(0,1), variance = "mu")
h <- g
h <- gstat(h, id = c("fulmar98","fulmar99"), 
	model = vgm(2.18, "Exp", 50000, 1.22))

# predict block means for blocks in ncp.grid$area (table 2; cokriging)
library(maptools)
areas.r = readShapePoly(system.file("external/ncp.shp", package="gstat"))
#areas.r <- as.SpatialRings.Shapes(areas.shp$Shapes, areas.shp$att.data$WSVGEB_)
coordinates(pr) = ~x+y
#pr.df = overlay(pr, areas.r, fn = mean)
pr.df = na.omit(as(aggregate(pr, areas.r), "data.frame"))
# match non-empty (and relevant) areas:
#areas = SpatialPolygonsDataFrame(areas.r[c(2,3,4,16),"WSVGEB_"], pr.df[c(1,2,3,5),])#,match.ID=F)
areas = SpatialPolygonsDataFrame(areas.r[c(1,2,12,7),"WSVGEB_"], 
	pr.df[c(1,2,3,5),], match.ID=F)
# areas ID's 0 1 2 14
sk = predict(g, areas)
cok = predict(h, areas)
spplot(cok, c(3,5), names.attr = c("1998", "1999"), 
	main = "Fulmaris glacialis, density estimates\n(by irregular block cokriging)")
sk = as.data.frame(sk)
cok = as.data.frame(cok)
print(data.frame(area = c(1,2,3,16), 
		SK98 = sk$fulmar98.pred, SE98 = sqrt(sk$fulmar98.var),
		SK99 = sk$fulmar99.pred, SE99 = sqrt(sk$fulmar99.var),
		CK98 = cok$fulmar98.pred, SE98 = sqrt(cok$fulmar98.var),
		CK99 = cok$fulmar99.pred, SE99 = sqrt(cok$fulmar99.var)),
	digits=3)
print(data.frame(area = c(1,2,3,16), 
		dSK = sk$fulmar99.pred - sk$fulmar98.pred, 
		SEdSK = sqrt(sk$fulmar98.var+sk$fulmar99.var),
		dCOK = cok$fulmar99.pred - cok$fulmar98.pred, 
		SEdCOK = sqrt(cok$fulmar98.var+cok$fulmar99.var
			 - 2*cok$cov.fulmar98.fulmar99)),
	digits=3)
