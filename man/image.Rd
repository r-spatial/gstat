% $Id: image.Rd,v 1.10 2007-11-16 12:59:47 edzer Exp $
\name{image}
\alias{image.data.frame}
\alias{image}
\alias{xyz2img}
\title{
Image Gridded Coordinates in Data Frame
}
\description{
Image gridded data, held in a data frame, keeping the
right aspect ratio for axes, and the right cell shape
}
\usage{
\method{image}{data.frame}(x, zcol = 3, xcol = 1, ycol = 2, asp = 1, ...)
xyz2img(xyz, zcol = 3, xcol = 1, ycol = 2, tolerance = 10 * .Machine$double.eps)
}
\arguments{
\item{x}{ data frame (or matrix) with x-coordinate,
y-coordinate, and z-coordinate in its columns }
\item{zcol}{ column number or name of z-variable }
\item{xcol}{ column number or name of x-coordinate }
\item{ycol}{ column number or name of y-coordinate }
\item{asp}{ aspect ratio for the x and y axes }
\item{...}{ arguments, passed to image.default }
\item{xyz}{data frame (same as \code{x})}
\item{tolerance}{ maximum allowed deviation for coordinats from being
exactly on a regularly spaced grid }
}
\value{
\link{image.data.frame} plots an image from gridded data, organized
in arbritrary order, in a data frame. It uses \link{xyz2img} and
\link{image.default} for this. In the S-Plus version, \link{xyz2img}
tries to make an image object with a size such that it will plot with
an equal aspect ratio; for the R version, image.data.frame uses the
\code{asp=1} argument to guarantee this.

\link{xyz2img} returns a list with components: \code{z}, a matrix
containing the z-values; \code{x}, the increasing coordinates of the
rows of \code{z}; \code{y}, the increasing coordinates of the columns
of \code{z}. This list is suitable input to \link{image.default}.
}
\note{
I wrote this function before I found out about \code{levelplot},
a Lattice/Trellis function that lets you control the aspect ratio by
the \code{aspect} argument, and that automatically draws a legend, and
therefore I now prefer levelplot over \code{image}. Plotting points
on a levelplots is probably done with providing a panel function and
using \code{lpoints}.

(for S-Plus only -- ) it is hard (if not impossible) to get exactly right
cell shapes (e.g., square for a square grid) without altering the size of
the plotting region, but this function tries hard to do so by extending
the image to plot in either x- or y-direction.  The larger the grid, the
better the approximation. Geographically correct images can be obtained
by modifiying \code{par("pin")}. Read the examples, image a 2 x 2 grid,
and play with \code{par("pin")} if you want to learn more about this.
}
\author{ Edzer Pebesma }
\examples{
library(sp)
data(meuse)
data(meuse.grid)
g <- gstat(formula=log(zinc)~1,locations=~x+y,data=meuse,model=vgm(1,"Exp",300))
x <- predict(g, meuse.grid)
image(x, 4, main="kriging variance and data points")
points(meuse$x, meuse$y, pch = "+")
# non-square cell test:
image(x[((x$y - 20) \%\% 80) == 0,], main = "40 x 80 cells")
image(x[((x$x - 20) \%\% 80) == 0,], main = "80 x 40 cells")
# the following works for square cells only:
oldpin <- par("pin")
ratio <- length(unique(x$x))/length(unique(x$y))
par(pin = c(oldpin[2]*ratio,oldpin[2]))
image(x, main="Exactly square cells, using par(pin)")
par(pin = oldpin)
library(lattice)
levelplot(var1.var~x+y, x, aspect = "iso", main = "kriging variance")
}

\keyword{dplot}

