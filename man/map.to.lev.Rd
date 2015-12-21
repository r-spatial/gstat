% $Id: map.to.lev.Rd,v 1.4 2006-02-10 19:03:27 edzer Exp $
\name{map.to.lev}
\alias{map.to.lev}
\title{ rearrange data frame for plotting with levelplot }
\description{ rearrange data frame for plotting with levelplot }
\usage{
map.to.lev(data, xcol = 1, ycol = 2, zcol = c(3, 4), ns = names(data)[zcol])
}
\arguments{
\item{data}{ data frame, e.g. output from \link{krige} or \link[gstat]{predict} }
\item{xcol}{ x-coordinate column number }
\item{ycol}{ y-coordinate column number }
\item{zcol}{ z-coordinate column number range }
\item{ns}{names of the set of z-columns to be viewed}
}
\value{
data frame with the following elements:
\item{x}{ x-coordinate for each row}
\item{y}{ y-coordinate for each row} 
\item{z}{ column vector with each of the elements in columns \code{zcol}
of \code{data} stacked } 
\item{name}{ factor; name of each of the stacked \code{z} columns } 
}
\seealso{ 
\link{image.data.frame}, \link{krige};
for examples see \link[gstat]{predict};
\code{levelplot} in package lattice.
}
\keyword{dplot}
