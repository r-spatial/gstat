% $Id: vgm.panel.Rd,v 1.3 2008-10-30 13:47:05 edzer Exp $
\name{vgm.panel.xyplot}
\alias{vgm.panel.xyplot}
\alias{panel.pointPairs}
\title{ panel functions for most of the variogram plots through lattice }
\description{
Variogram plots contain symbols and lines; more control over them can
be gained by writing your own panel functions, or extending the ones
described here; see examples.
}
\usage{
vgm.panel.xyplot(x, y, subscripts, type = "p", pch = plot.symbol$pch, 
    col, col.line = plot.line$col, col.symbol = plot.symbol$col, 
    lty = plot.line$lty, cex = plot.symbol$cex, ids, lwd = plot.line$lwd, 
    model = model, direction = direction, labels, shift = shift, mode = mode, ...) 
panel.pointPairs(x, y, type = "p", pch = plot.symbol$pch, col, col.line = 
	plot.line$col, col.symbol = plot.symbol$col, lty = plot.line$lty, 
	cex = plot.symbol$cex, lwd = plot.line$lwd, pairs = pairs, 
	line.pch = line.pch, ...) 
}
\arguments{
\item{x}{ x coordinates of points in this panel}
\item{y}{ y coordinates of points in this panel}
\item{subscripts }{ subscripts of points in this panel}
\item{type}{ plot type: "l" for connected lines }
\item{pch}{ plotting symbol }
\item{col}{ symbol and line color (if set) }
\item{col.line}{ line color }
\item{col.symbol}{ symbol color }
\item{lty}{ line type for variogram model }
\item{cex}{ symbol size }
\item{ids}{ gstat model ids }
\item{lwd}{ line width }
\item{model}{ variogram model }
\item{direction}{ direction vector \code{c(dir.horizontal, dir.ver)}}
\item{labels}{ labels to plot next to points }
\item{shift}{ amount to shift the label right of the symbol }
\item{mode}{ to be set by calling function only }
\item{line.pch}{ symbol type to be used for point of selected point pairs, e.g.
to highlight point pairs with distance close to zero }
\item{pairs}{ two-column matrix with pair indexes to be highlighted }
\item{...}{ parameters that get passed to \link{lpoints} }
}
\value{
ignored; the enclosing function returns a plot of class \code{trellis}
}
\references{ \url{http://www.gstat.org/} }
\author{ Edzer Pebesma }
\seealso{ \link{plot.gstatVariogram}, \link{vgm}}
\examples{ 
library(sp)
data(meuse)
coordinates(meuse) <- c("x", "y")
library(lattice)
mypanel = function(x,y,...) { 
	vgm.panel.xyplot(x,y,...)
	panel.abline(h=var(log(meuse$zinc)), color = 'red')
}
plot(variogram(log(zinc)~1,meuse), panel = mypanel)
}
\keyword{models}
