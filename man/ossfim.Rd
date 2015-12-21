% $Id: ossfim.Rd,v 1.3 2006-02-10 19:03:27 edzer Exp $
\name{ossfim}
\alias{ossfim}
\title{ Kriging standard errors as function of grid spacing and block size}
\description{ Calculate, for a given variogram model, ordinary block
kriging standard errors as a function of sampling spaces and block sizes }
\usage{
ossfim(spacings = 1:5, block.sizes = 1:5, model, nmax = 25, debug = 0)
}
\arguments{
\item{spacings}{range of grid (data) spacings to be used}
\item{block.sizes}{ range of block sizes to be used}
\item{model}{variogram model, output of \code{vgm}}
\item{nmax}{set the kriging neighbourhood size}
\item{debug}{debug level; set to 32 to see a lot of output}
}

\value{ data frame with columns \code{spacing} (the grid spacing),
\code{block.size} (the block size), and \code{kriging.se} (block kriging
standard error) }

\references{ 
Burrough, P.A., R.A. McDonnell (1999) Principles of Geographical
Information Systems. Oxford University Press (e.g., figure 10.11 on
page 261)

Burgess, T.M., R. Webster, A.B. McBratney (1981) Optimal interpolation and
isarithmic mapping of soil properties. IV Sampling strategy.  The journal
of soil science 32(4), 643-660.

McBratney, A.B., R. Webster (1981) The design of optimal sampling schemes
for local estimation and mapping of regionalized variables: 2 program
and examples. Computers and Geosciences 7: 335-365.
}

\author{ Edzer Pebesma }
\note{ The idea is old, simple, but still of value. If you want to map
a variable with a given accuracy, you will have to sample it. Suppose
the variogram of the variable is known. Given a regular sampling scheme,
the kriging standard error decreases when either (i) the data spacing
is smaller, or (ii) predictions are made for larger blocks. This function
helps quantifying this relationship. Ossfim probably refers to ``optimal
sampling scheme for isarithmic mapping''.
}

\seealso{
\link{krige}
}
\examples{
\dontrun{
x <- ossfim(1:15,1:15, model = vgm(1,"Exp",15))
library(lattice)
levelplot(kriging.se~spacing+block.size, x, 
  main = "Ossfim results, variogram 1 Exp(15)")
}
# if you wonder about the decrease in the upper left corner of the graph,
# try the above with nmax set to 100, or perhaps 200.
}
\keyword{models}
