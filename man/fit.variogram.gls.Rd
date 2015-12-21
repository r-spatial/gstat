% $Id: fit.variogram.gls,v 1.4 2009-02-20 13:53:38 edzer Exp $
\name{fit.variogram.gls}
\alias{fit.variogram.gls}
\title{ GLS fitting of variogram parameters }
\description{
Fits variogram parameters (nugget, sill, range) to variogram cloud,
using GLS (generalized least squares) fitting. Only for direct variograms.
}
\usage{
fit.variogram.gls(formula, data, model, maxiter = 30, 
		eps = .01, trace = TRUE, ignoreInitial = TRUE, cutoff = Inf,
		plot = FALSE)
}
\arguments{
\item{formula}{formula defining the response vector and (possible) 
regressors; in case of absence of regressors, use e.g. \code{z~1}}
\item{data}{object of class Spatial}
\item{model}{variogram model to be fitted, output of \code{vgm}}
\item{maxiter}{maximum number of iterations}
\item{eps}{ convergence criterium }
\item{trace}{ logical; if TRUE, prints parameter trace}
\item{ignoreInitial}{ logical; 
if FALSE, initial parameter are taken from model; 
if TRUE, initial values of model are
ignored and taken from variogram cloud: 
nugget: \code{mean(y)/2}, sill: \code{mean(y)/2}, range \code{median(h0)/4}
with \code{y} the semivariance cloud value and \code{h0} the distances }
\item{cutoff}{maximum distance up to which point pairs are taken into
consideration}
\item{plot}{logical; if TRUE, a plot is returned with variogram cloud and
fitted model; else, the fitted model is returned.}
}

\value{ an object of class "variogramModel"; see \link{fit.variogram}; if
\code{plot} is TRUE, a plot is returned instead. }

\references{ 
Mueller, W.G., 1999: Least-squares fitting from the variogram
cloud. Statistics \& Probability Letters, 43, 93-98.

Mueller, W.G., 2007: Collecting Spatial Data. Springer, Heidelberg.
}

\author{ Edzer Pebesma }
\note{
Inspired by the code of Mihael Drinovac, which was again inspired by
code from Ernst Glatzer, author of package vardiag.
}
\seealso{
\link{fit.variogram},
}
\examples{
library(sp)
data(meuse)
coordinates(meuse) = ~x+y
\dontrun{
fit.variogram.gls(log(zinc)~1, meuse[1:40,], vgm(1, "Sph", 900,1))
}
}
\keyword{models}
