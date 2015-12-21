% $Id: variogramLine.Rd,v 1.3 2008-08-19 07:27:02 edzer Exp $
\name{vgmArea}
\alias{vgmArea}
\title{ point-point, point-area or area-area semivariance }
\description{ Compute point-point, point-area or area-area variogram values from point model }
\usage{
vgmArea(x, y = x, vgm, ndiscr = 16, verbose = FALSE, covariance = TRUE)
}
\arguments{
\item{x}{ object of class \link[sp]{SpatialPoints} or \link[sp]{SpatialPolygons}}
\item{y}{ object of class \link[sp]{SpatialPoints} or \link[sp]{SpatialPolygons}}
\item{vgm}{ variogram model, see \link{vgm}}
\item{ndiscr}{ number of points to discretize an area, using \link[sp]{spsample}}
\item{verbose}{ give progress bar }
\item{covariance}{ logical; compute covariances, rather than semivariances? }
}
\value{ semivariance or covariance matrix of dimension \code{length(x)} x \code{lenght(y)}}

\author{ Edzer Pebesma }

\examples{
library(sp)
demo(meuse, ask = FALSE, echo = FALSE)
vgmArea(meuse[1:5,], vgm = vgm(1, "Exp", 1000)) # point-point
vgmArea(meuse[1:5,], meuse.area, vgm = vgm(1, "Exp", 1000)) # point-area
}

\keyword{models}
