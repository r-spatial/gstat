% $Id: hscat.Rd,v 1.3 2008-02-04 10:06:44 edzer Exp $
\name{hscat}
\alias{hscat}
\title{
Produce h-scatterplot
}
\description{
Produces h-scatterplots, where point pairs having specific separation distances
are plotted. This function is a wrapper around xyplot.
}
\usage{
hscat(formula, data, breaks, pch = 3, cex = .6, mirror = FALSE, 
	variogram.alpha = 0, as.table = TRUE,...)
}
\arguments{
\item{formula}{ specifies the dependent variable }
\item{data}{ data where the variable in formula is resolved }
\item{breaks}{ distance class boundaries }
\item{pch}{ plotting symbol }
\item{cex}{ plotting symbol size }
\item{mirror}{ logical; duplicate all points mirrored along x=y?
(note that correlations are those of the points plotted) }
\item{variogram.alpha}{ parameter to be passed as alpha parameter to
\link{variogram}; if alpha is specified it will only affect xyplot by
being passed through ...}
\item{as.table}{logical; if \code{TRUE}, panels plot top-to-bottom }
\item{...}{ parameters, passed to variogram and xyplot }
}
\value{
an object of class trellis; normally the h scatter plot
}
\author{ Edzer Pebesma }
\note{
Data pairs are plotted once, so the h-scatterplot are not symmetric.
}

\references{ \url{http://www.gstat.org/}

Pebesma, E.J., 2004. Multivariable geostatistics in S: the gstat package.
Computers and Geosciences, 30: 683-691.

}
\examples{
library(sp)
data(meuse)
coordinates(meuse) = ~x+y
hscat(log(zinc)~1, meuse, c(0, 80, 120, 250, 500, 1000))
}

\keyword{models}
