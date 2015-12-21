% $Id: variogramLine.Rd,v 1.3 2008-08-19 07:27:02 edzer Exp $
\name{variogramLine}
\alias{variogramLine}
\alias{getGammas}
\title{ Semivariance Values For a Given Variogram Model }
\description{ Generates a semivariance values given a variogram model }
\usage{
variogramLine(object, maxdist, n = 200, min = 1.0e-6 * maxdist, 
	dir = c(1,0,0), covariance = FALSE, ..., dist_vector, debug.level = 0)
}
\arguments{
\item{object}{ variogram model for which we want semivariance function values }
\item{maxdist}{ maximum distance for which we want semivariance values }
\item{n}{ number of points }
\item{min}{ minimum distance; a value slightly larger than zero is usually
used to avoid the discontinuity at distance zero if a nugget component is
present }
\item{dir}{ direction vector: unit length vector pointing the direction in 
x (East-West), y (North-South) and z (Up-Down) }
\item{covariance}{logical; if TRUE return covariance values, otherwise
return semivariance values }
\item{...}{ignored}
\item{dist_vector}{numeric vector or matrix with distance values}
\item{debug.level}{gstat internal debug level}
}
\value{
a data frame of dimension (\code{n} x 2), with columns distance and gamma
(semivariances or covariances), or in case \code{dist_vector} is a matrix, a
conforming matrix with semivariance/covariance values is returned.
}
\note{variogramLine is used to generate data for plotting a variogram model.}
\author{ Edzer Pebesma }
\seealso{ \link{plot.gstatVariogram}}

\examples{
variogramLine(vgm(5, "Exp", 10, 5), 10, 10)
# anisotropic variogram, plotted in E-W direction:
variogramLine(vgm(1, "Sph", 10, anis=c(0,0.5)), 10, 10)
# anisotropic variogram, plotted in N-S direction:
variogramLine(vgm(1, "Sph", 10, anis=c(0,0.5)), 10, 10, dir=c(0,1,0))
variogramLine(vgm(1, "Sph", 10, anis=c(0,0.5)), dir=c(0,1,0), dist_vector = 0.5)
variogramLine(vgm(1, "Sph", 10, anis=c(0,0.5)), dir=c(0,1,0), dist_vector = c(0, 0.5, 0.75))
}

\keyword{models}
