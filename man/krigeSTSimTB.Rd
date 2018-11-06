\name{krigeSTSimTB}
\alias{krigeSTSimTB}

\title{conditional/unconditional spatio-temporal simulation}
\description{
conditional/unconditional spatio-temporal simulation based on turning bands
}
\usage{
krigeSTSimTB(formula, data, newdata, modelList, nsim, progress = TRUE, 
             nLyrs = 500, tGrid = NULL, sGrid = NULL, ceExt = 2, nmax = Inf)
}
\arguments{
  \item{formula}{the formula of the kriging predictor}
  \item{data}{conditioning data}
  \item{newdata}{locations in space and time where the simulation is carried out}
  \item{modelList}{the spatio-temporal variogram (from \code{\link{vgmST}}) defining the spatio-temporal covariance structure of the simulated Gaussian random field}
  \item{nsim}{number of simulations}
  \item{progress}{boolean; whether the progress should be shown in progress bar}
  \item{nLyrs}{number of layers used in the turning bands approach (default = 500) }
  \item{tGrid}{optional explicit temporal griding that shall be used}
  \item{sGrid}{optional explicit spatial griding that shall be used}
  \item{ceExt}{expansion in the circulant embedding, defaults to 2}
  \item{nmax}{number of nearest neighbours that shall e used, defaults to 'Inf' meaning all available points are used}
}
\value{a spatio-temporal data frame with \code{nSim} simulations}
\references{
Turning bands

Lantuejoul, C. (2002) Geostatistical Simulation: Models and Algorithms. Springer.

Matheron, G. (1973). The intrinsic random functions and their applications. Adv. Appl. Probab., 5, 439-468.

Strokorb, K., Ballani, F., and Schlather, M. (2014) Tail correlation functions of max-stable processes: Construction principles, recovery and diversity of some mixing max-stable processes with identical TCF. Extremes, Submitted.

Turning layers

Schlather, M. (2011) Construction of covariance functions and unconditional simulation of random fields. In Porcu, E., Montero, J.M. and Schlather, M., Space-Time Processes and Challenges Related to Environmental Problems. New York: Springer.
}

\author{Benedikt Graeler}

\seealso{\code{\link{krigeSimCE}}}
\examples{
# see demo('circEmbeddingMeuse')
}