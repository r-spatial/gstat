% $Id: plot.pointPairs.Rd,v 1.5 2006-12-12 20:44:07 edzer Exp $
\name{plot.pointPairs}
\alias{plot.pointPairs}
\title{
Plot a point pairs, identified from a variogram cloud
}
\description{
Plot a point pairs, identified from a variogram cloud
}
\usage{
\method{plot}{pointPairs}(x, data, xcol = data$x, ycol = data$y, xlab = "x coordinate", 
ylab = "y coordinate", col.line = 2, line.pch = 0, main = "selected point pairs", ...) 
}
\arguments{
\item{x}{ object of class "pointPairs", obtained from the function 
\link{plot.variogramCloud}, containing point pair indices }
\item{data}{ data frame to which the indices refer (from which the
variogram cloud was calculated) }
\item{xcol}{ numeric vector with x-coordinates of data }
\item{ycol}{ numeric vector with y-coordinates of data }
\item{xlab}{ x-axis label }
\item{ylab}{ y-axis label }
\item{col.line}{ color for lines connecting points }
\item{line.pch}{ if non-zero, symbols are also plotted at the middle of
line segments, to mark lines too short to be visible on the plot;
the color used is \code{col.line}; the value passed to this argument
will be used as plotting symbol (pch) }
\item{main}{ title of plot }
\item{...}{ arguments, further passed to \code{xyplot}}
}
\value{
plots the data locations, with lines connecting the point pairs identified
(and refered to by indices in) x
}
\references{ \url{http://www.gstat.org}}
\author{ Edzer Pebesma }
\seealso{
\link{plot.variogramCloud}
}
\examples{
### The following requires interaction, and is therefore outcommented
#data(meuse)
#coordinates(meuse) = ~x+y
#vgm1 <- variogram(log(zinc)~1, meuse, cloud = TRUE)
#pp <- plot(vgm1, id = TRUE)
### Identify the point pairs
#plot(pp, data = meuse) # meuse has x and y as coordinates
}
\keyword{dplot}
