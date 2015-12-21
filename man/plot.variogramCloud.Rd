% $Id: plot.variogramCloud.Rd,v 1.6 2006-12-12 20:44:07 edzer Exp $
\name{plot.variogramCloud}
\alias{plot.variogramCloud}
\title{
Plot and Identify Data Pairs on Sample Variogram Cloud
}
\description{
Plot a sample variogram cloud, possibly with identification of
individual point pairs
}
\usage{
\method{plot}{variogramCloud}(x, identify = FALSE, digitize = FALSE, xlim, ylim, xlab, ylab, 
	keep = FALSE, ...)
}
\arguments{
\item{x}{ object of class \code{variogramCloud}}
\item{identify}{ logical; if TRUE, the plot allows identification of
a series of individual point pairs that correspond to individual variogram
cloud points (use left mouse button to select; right mouse button ends) }
\item{digitize}{ logical; if TRUE, select point pairs by digitizing a
region with the mouse (left mouse button adds a point, right mouse button
ends) }
\item{xlim}{ limits of x-axis }
\item{ylim}{ limits of y-axis }
\item{xlab}{ x axis label }
\item{ylab}{ y axis label }
\item{keep}{ logical; if TRUE and \code{identify} is TRUE, the labels
identified and their position are kept and glued to object x, which is
returned. Subsequent calls to plot this object will now have the labels
shown, e.g. to plot to hardcopy }
\item{...}{ parameters that are passed through to \link{plot.gstatVariogram}
(in case of identify = FALSE) or to plot (in case of identify = TRUE) }
}
\value{
If \code{identify} or \code{digitize} is TRUE, a data frame of class
\code{pointPairs} with in its rows the point pairs identified (pairs of
row numbers in the original data set); if identify is F, a plot of the
variogram cloud, which uses \link{plot.gstatVariogram}

If in addition to \code{identify}, \code{keep} is also TRUE, an object
of class \code{variogramCloud} is returned, having attached to it attributes
"sel" and "text", which will be used in subsequent calls to plot.variogramCloud
with \code{identify} set to FALSE, to plot the text previously identified.

If in addition to \code{digitize}, \code{keep} is also TRUE, an object of
class \code{variogramCloud} is returned, having attached to it attribute
"poly", which will be used in subsequent calls to plot.variogramCloud
with \code{digitize} set to FALSE, to plot the digitized line.

In both of the \code{keep = TRUE} cases, the attribute \code{ppairs} of
class \code{pointPairs} is present, containing the point pairs identified.
}
\references{ 
\url{http://www.gstat.org/}
}
\author{ Edzer Pebesma }
\seealso{
\link{variogram},
\link{plot.gstatVariogram},
\link{plot.pointPairs},
\link{identify},
\link{locator}
}
\examples{
library(sp)
data(meuse)
coordinates(meuse) = ~x+y
plot(variogram(log(zinc)~1, meuse, cloud=TRUE))
## commands that require interaction:
# x <- variogram(log(zinc)~1, loc=~x+y, data=meuse, cloud=TRUE)
# plot(plot(x, identify = TRUE), meuse)
# plot(plot(x, digitize = TRUE), meuse)
}

\keyword{dplot}
