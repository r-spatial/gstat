\name{wind}
\alias{wind}
\alias{wind.loc}
\title{Ireland wind data, 1961-1978}
\description{
Daily average wind speeds for 1961-1978 at 12 synoptic meteorological 
stations in the Republic of Ireland (Haslett and raftery 1989).
Wind speeds are in knots (1 knot = 0.5418 m/s), at each of the 
stations in the order given in Fig.4 of Haslett and Raftery (1989, see below)
}
\format{
data.frame \code{wind} contains the following columns:
\describe{
\item{year}{year, minus 1900}
\item{month}{month (number) of the year}
\item{day }{day}
\item{RPT}{average wind speed in knots at station RPT}
\item{VAL}{average wind speed in knots at station VAL}
\item{ROS}{average wind speed in knots at station ROS}
\item{KIL}{average wind speed in knots at station KIL}
\item{SHA}{average wind speed in knots at station SHA}
\item{BIR}{average wind speed in knots at station BIR}
\item{DUB}{average wind speed in knots at station DUB}
\item{CLA}{average wind speed in knots at station CLA}
\item{MUL}{average wind speed in knots at station MUL}
\item{CLO}{average wind speed in knots at station CLO}
\item{BEL}{average wind speed in knots at station BEL}
\item{MAL}{average wind speed in knots at station MAL}
}
data.frame \code{wind.loc} contains the following columns:
\describe{
\item{Station}{Station name}
\item{Code }{Station code}
\item{Latitude}{Latitude, in DMS, see examples below}
\item{Longitude}{Longitude, in DMS, see examples below}
\item{MeanWind}{mean wind for each station, metres per second }
}
}
\usage{
data(wind)
}
\author{ Adrian Raftery; imported to R by Edzer Pebesma }
\references{ 
These data were analyzed in detail in the following article:

   Haslett, J. and Raftery, A. E. (1989). Space-time Modelling with
   Long-memory Dependence: Assessing Ireland's Wind Power Resource
   (with Discussion). Applied Statistics 38, 1-50.

and in many later papers on space-time analysis, for example:

   Tilmann Gneiting, Marc G. Genton, Peter Guttorp: Geostatistical
   Space-Time Models, Stationarity, Separability and Full symmetry.
   Ch. 4 in: B. Finkenstaedt, L. Held, V. Isham, Statistical Methods
   for Spatio-Temporal Systems.
}
\note{ This data set comes with the following message:
``Be aware that the dataset is 532494 bytes long (thats over half a
Megabyte).  Please be sure you want the data before you request it.''

The data were obtained on Oct 12, 2008, from:
http://www.stat.washington.edu/raftery/software.html
The data are also available from statlib.

Locations of 11 of the stations (ROS, Rosslare has
been thrown out because it fits poorly the spatial
correlations of the other stations) were obtained from:
http://www.stat.washington.edu/research/reports/2005/tr475.pdf

Roslare lat/lon was obtained from google maps, location Roslare. The mean
wind value for Roslare comes from Fig. 1 in the original paper.

Haslett and Raftery proposed to use a sqrt-transform to stabilize the variance.
}
\keyword{datasets}
\examples{
data(wind)
summary(wind)
wind.loc
library(sp) # char2dms
wind.loc$y = as.numeric(char2dms(as.character(wind.loc[["Latitude"]])))
wind.loc$x = as.numeric(char2dms(as.character(wind.loc[["Longitude"]])))
coordinates(wind.loc) = ~x+y

\dontrun{
# fig 1:
library(maps)
library(mapdata)
map("worldHires", xlim = c(-11,-5.4), ylim = c(51,55.5))
points(wind.loc, pch=16)
text(coordinates(wind.loc), pos=1, label=wind.loc$Station)
}

wind$time = ISOdate(wind$year+1900, wind$month, wind$day)
# time series of e.g. Dublin data:
plot(DUB~time, wind, type= 'l', ylab = "windspeed (knots)", main = "Dublin")

# fig 2:
#wind = wind[!(wind$month == 2 & wind$day == 29),]
wind$jday = as.numeric(format(wind$time, '\%j'))
windsqrt = sqrt(0.5148 * as.matrix(wind[4:15]))
Jday = 1:366
windsqrt = windsqrt - mean(windsqrt)
daymeans = sapply(split(windsqrt, wind$jday), mean)
plot(daymeans ~ Jday)
lines(lowess(daymeans ~ Jday, f = 0.1))

# subtract the trend:
meanwind = lowess(daymeans ~ Jday, f = 0.1)$y[wind$jday]
velocity = apply(windsqrt, 2, function(x) { x - meanwind })

# match order of columns in wind to Code in wind.loc:
pts = coordinates(wind.loc[match(names(wind[4:15]), wind.loc$Code),])

# fig 3, but not really yet...
dists = spDists(pts, longlat=TRUE)
corv = cor(velocity)
sel = !(as.vector(dists) == 0)
plot(as.vector(corv[sel]) ~ as.vector(dists[sel]),
	xlim = c(0,500), ylim = c(.4, 1), xlab = "distance (km.)", 
	ylab = "correlation") 
# plots all points twice, ignores zero distance 
# now really get fig 3:
ros = rownames(corv) == "ROS"
dists.nr = dists[!ros,!ros]
corv.nr = corv[!ros,!ros]
sel = !(as.vector(dists.nr) == 0)
plot(as.vector(corv.nr[sel]) ~ as.vector(dists.nr[sel]), pch = 3,
	xlim = c(0,500), ylim = c(.4, 1), xlab = "distance (km.)", 
	ylab = "correlation") 
# add outlier:
points(corv[ros,!ros] ~ dists[ros,!ros], pch=16, cex=.5)
xdiscr = 1:500
# add correlation model:
lines(xdiscr, .968 * exp(- .00134 * xdiscr))
}
