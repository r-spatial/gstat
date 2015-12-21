# $Id: krige.R,v 1.5 2007-02-27 22:09:31 edzer Exp $
library(sp)
data(meuse)
coordinates(meuse) = ~x+y
data(meuse.grid)
gridded(meuse.grid) = ~x+y

# ordinary kriging
v <- variogram(log(zinc)~1, meuse)
m <- fit.variogram(v, vgm(1, "Sph", 300, 1))
plot(v, model = m)
lzn.kr <- krige(formula = log(zinc)~1, meuse, meuse.grid, model = m)


pl1 <- spplot(lzn.kr[1], main = "ordinary kriging prediction of log-zinc")
lzn.kr$se = sqrt(lzn.kr$var1.var)
pl2 <- spplot(lzn.kr["se"], main = "ordinary kriging prediction error")

# universal kriging
v <- variogram(log(zinc)~sqrt(dist), meuse)
m <- fit.variogram(v, vgm(1, "Exp", 300, 1))
plot(v, model = m)
lzn.kr <- krige(log(zinc)~sqrt(dist), meuse, meuse.grid, model = m)
pl3 <- spplot(lzn.kr[1], main = "universal kriging prediction of log-zinc")
lzn.kr$se = sqrt(lzn.kr$var1.var)
pl4 <- spplot(lzn.kr["se"], main = "universal kriging prediction error")
print(pl1, split = c(1,1,2,2), more = T)
print(pl2, split = c(1,2,2,2), more = T)
print(pl3, split = c(2,1,2,2), more = T)
print(pl4, split = c(2,2,2,2))
