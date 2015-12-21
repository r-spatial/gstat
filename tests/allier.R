# Sytze de Bruin's post to r-sig-geo, Nov 16, 2015:
library(sp)
library(gstat)

# some data
x <- c(215, 330, 410, 470, 545)
y <- c(230, 310, 330, 340, 365)
fc <- c(0.211, 0.251, 0.281, 0.262, 0.242)
por <- c(0.438, 0.457, 0.419, 0.430, 0.468)
Allier <- data.frame(x, y, fc, por)
coordinates(Allier) = ~x+y

# gstat object for co-kriging using linear model of co-regionalization
g <- gstat(id=c("fc"), formula=fc~1, data=Allier,
           model=vgm(0.00247, "Sph", 480, 0.00166))
g <- gstat(g, id="por", formula=por~1, data=Allier,
           model=vgm(0.00239, "Sph", 480, 0.00118))
g <- gstat(g, id=c("fc", "por"),
           model=vgm(0.00151, "Sph", 480, -0.00124))

# predict at single point
g$set = list(choleski = 0) # LDL'
A <- predict(g, SpatialPoints(data.frame(x=450, y=350)), debug = 32)
g$set = list(choleski = 1) # Choleski
B <- predict(g, SpatialPoints(data.frame(x=450, y=350)), debug = 32)
identical(A,B)

# TRUE
