% $Id: fit.variogram.reml.Rd,v 1.4 2009-02-20 13:53:38 edzer Exp $
\name{fit.variogram.reml}
\alias{fit.variogram.reml}
\title{ REML Fit Direct Variogram Partial Sills to Data }
\description{
Fit Variogram Sills to Data, using REML (only for direct variograms;
not for cross variograms)
}
\usage{
fit.variogram.reml(formula, locations, data, model, debug.level = 1, set, degree = 0)
}
\arguments{
\item{formula}{formula defining the response vector and (possible) 
regressors; in case of absence of regressors, use e.g. \code{z~1}}
\item{locations}{ spatial data locations; a formula with the
coordinate variables in the right hand (dependent variable) side. }
\item{data}{data frame where the names in formula and locations 
are to be found}
\item{model}{variogram model to be fitted, output of \code{vgm}}
\item{debug.level}{debug level; set to 65 to see the iteration trace and
log likelihood}
\item{set}{additional options that can be set; use \code{set=list(iter=100)}
to set the max. number of iterations to 100.  }
\item{degree}{order of trend surface in the location, between 0 and 3}
}

\value{ an object of class "variogramModel"; see \link{fit.variogram} }

\references{ 
Christensen, R. Linear models for multivariate, Time Series,
and Spatial Data, Springer, NY, 1991. 

Kitanidis, P., Minimum-Variance Quadratic Estimation of Covariances of
Regionalized Variables, Mathematical Geology 17 (2), 195--208, 1985 }

\author{ Edzer Pebesma }
\note{
This implementation only uses REML fitting of sill parameters. For each
iteration, an \eqn{n \times n}{n x n} matrix is inverted, with $n$ the number of
observations, so for large data sets this method becomes 
demanding. I guess there is much more to likelihood variogram fitting in 
package \code{geoR}, and probably also in \code{nlme}.
}
\seealso{
\link{fit.variogram},
}
\examples{
library(sp)
data(meuse)
fit.variogram.reml(log(zinc)~1, ~x+y, meuse, model = vgm(1, "Sph", 900,1))
}
\keyword{models}
