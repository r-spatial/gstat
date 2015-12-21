% $Id: spplot.vcov.Rd,v 1.2 2007-11-16 12:59:47 edzer Exp $
\name{spplot.vcov}
\alias{spplot.vcov}
\title{ Plot map matrix of prediction error variances and covariances }
\description{
Plot map matrix of prediction error variances and covariances
}
\usage{
spplot.vcov(x, ...)
}
\arguments{
\item{x}{ Object of class SpatialPixelsDataFrame or SpatialGridDataFrame,
resulting from a krige call with multiple variables (cokriging }
\item{...}{remaining arguments passed to spplot }
}
\value{ The plotted object, of class trellis; see \code{spplot} in
package \pkg{sp}. }
\author{ Edzer Pebesma }
\keyword{dplot}

