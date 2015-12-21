\name{variogramSurface}
\alias{variogramSurface}
\title{Semivariance values for a given spatio-temporal variogram model}
\description{
Generates a surface of semivariance values given a spatio-temporal variogram model (one of separable, productSum, sumMetric, simpleSumMetric or metric)
}
\usage{
variogramSurface(model, dist_grid, ...)
}

\arguments{
  \item{model}{
  A spatio-temporal variogram model generated through \code{\link{vgmST}} or \code{\link{fit.StVariogram}}.
}
  \item{dist_grid}{
  A data.frame with two columns: \code{spacelag} and \code{timelag}.
}
  \item{\dots}{
Additional arguments passed on to the underlying variogram functions.
}
}
\value{
A data.frame with columns \code{spacelag}, \code{timelag} and \code{gamma}.
}
\author{
Benedikt Graeler
}
\seealso{
See \code{\link{variogramLine}} for the spatial version and \code{\link{fit.StVariogram}} for the estimation of spatio-temporal variograms.
}
\examples{
separableModel <- vgmST("separable", 
                        space=vgm(0.86, "Exp", 476, 0.14),
                        time =vgm(   1, "Exp",   3, 0),
                        sill=113)

data(vv)

if(require(lattice)) {
plot(vv, separableModel, wireframe=TRUE, all=TRUE)
}

# plotting of sample and model variogram
plot(vv, separableModel)

}
\keyword{models}
\keyword{spatio-temporal}
