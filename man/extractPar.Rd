\name{extractPar}
\alias{extractPar}
\alias{extractParNames}

\title{
Extracting parameters and their names from a spatio-temporal variogram model
}
\description{
All spatio-temporal variogram models have a different set of parameters. These functions extract the parameters and their names from the spatio-temporal variogram model. Note, this function is as well used to pass the parameters to the optim function. The arguments lower and upper passed to optim should follow the same structure. 
}
\usage{
extractPar(model)
extractParNames(model)
}

\arguments{
  \item{model}{a spatio-temporal variogram model from \code{\link{vgmST}}}
}

\value{
A named numeric vector of parameters or a vector of characters holding the parameters' names.
}

\author{
Benedikt Graeler
}


\seealso{
\code{\link{fit.StVariogram}} and \code{\link{vgmST}}
}

\examples{
sumMetricModel <- vgmST("sumMetric",
                        space=vgm(30, "Sph", 200,  6),
                        time =vgm(30, "Sph",  15,  7),
                        joint=vgm(60, "Exp",  84, 22),
                        stAni=100)

extractPar(sumMetricModel)
extractParNames(sumMetricModel)
}