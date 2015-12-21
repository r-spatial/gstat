\name{tull}
\alias{tull}
\alias{TULLNREG}
\alias{tull36}
\alias{Chlorid92}
\encoding{utf-8}
\title{Südliche Tullnerfeld data set}
\description{
The Südliche Tullnerfeld is a part of the Danube river basin in
central Lower Austria and due to its homogeneous aquifer well suited
for a model-oriented geostatistical analysis. It contains 36 official
water quality measurement stations, which are irregularly spread over the
region. 
}
\format{
  The data frames contain the following columns:
  \describe{
  \item{x}{X location in meter}
  \item{y}{Y location in meter}
\item{S411}{Station name}
\item{S429}{Station name}
\item{S849}{Station name}
\item{S854}{Station name}
\item{S1502}{Station name}
\item{S1584}{Station name}
\item{S1591}{Station name}
\item{S2046}{Station name}
\item{S2047}{Station name}
\item{S2048}{Station name}
\item{S2049}{Station name}
\item{S2051}{Station name}
\item{S2052}{Station name}
\item{S2053}{Station name}
\item{S2054}{Station name}
\item{S2055}{Station name}
\item{S2057}{Station name}
\item{S2058}{Station name}
\item{S2059}{Station name}
\item{S2060}{Station name}
\item{S2061}{Station name}
\item{S2062}{Station name}
\item{S2063}{Station name}
\item{S2064}{Station name}
\item{S2065}{Station name}
\item{S2066}{Station name}
\item{S2067}{Station name}
\item{S2070}{Station name}
\item{S2071}{Station name}
\item{S2072}{Station name}
\item{S2128}{Station name}
\item{S5319}{Station name}
\item{S5320}{Station name}
\item{S5321}{Station name}
\item{S5322}{Station name}
\item{S5323}{Station name}
  }
}
\usage{
data(tull)
}
\references{ 
Werner G. Müller, Collecting Spatial Data, 3rd edition.
Springer Verlag, Heidelberg, 2007
}
\note{
This data set was obtained on May 6, 2008 from
\url{ http://www.ifas.jku.at/e5361/index_ger.html }.
The author of the book that uses it is found at:
\url{ http://www.ifas.jku.at/e2571/e2604/index_ger.html }
}
\keyword{datasets}
\examples{
data(tull)

# TULLNREG = read.csv("TULLNREG.csv")

# I modified tulln36des.csv, such that the first line only contained: x,y
# resulting in row.names that reflect the station ID, as in
# tull36 = read.csv("tulln36des.csv")

# Chlorid92 was read & converted by:
#Chlorid92=read.csv("Chlorid92.csv")
#Chlorid92$Datum = as.POSIXct(strptime(Chlorid92$Datum, "%d.%m.%y"))

summary(tull36)
summary(TULLNREG)
summary(Chlorid92)

# stack & join data to x,y,Date,Chloride form:
cl.st = stack(Chlorid92[-1])
names(cl.st) = c("Chloride", "Station")
cl.st$Date = rep(Chlorid92$Datum, length(names(Chlorid92))-1)
cl.st$x = tull36[match(cl.st[,"Station"], row.names(tull36)), "x"]
cl.st$y = tull36[match(cl.st[,"Station"], row.names(tull36)), "y"]
# library(lattice)
# xyplot(Chloride~Date|Station, cl.st)
# xyplot(y~x|Date, cl.st, asp="iso", layout=c(16,11))
summary(cl.st)
plot(TULLNREG, pch=3, asp=1)
points(y~x, cl.st, add=TRUE, pch=16)
}
