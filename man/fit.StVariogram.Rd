\name{fit.StVariogram}
\alias{fit.StVariogram}
\title{Fit a spatio-temporal sample variogram to a sample variogram}
\description{
Fits a spatio-temporal variogram of a given type to spatio-temporal sample variogram.
}
\usage{
fit.StVariogram(object, model, ...,  method = "L-BFGS-B",
	lower, upper, fit.method = 6, stAni=NA, wles)
}
\arguments{
  \item{object}{The spatio-temporal sample variogram. Typically output from \code{\link{variogramST}}}
  \item{model}{The desired spatio-temporal model defined through \code{\link{vgmST}}.}
  \item{\dots}{further arguments passed to \code{\link{optim}}. \code{\link{extractParNames}} provides the parameter structure of spatio-temporal variogram models that help to provide sensible \code{upper} and \code{lower} limits.}
  \item{lower}{Lower limits used by optim. If missing, the smallest well defined values are used (mostly near 0).}
  \item{upper}{Upper limits used by optim. If missing, the largest well defined values are used (mostly \code{Inf}).}
  \item{method}{fit method, pass to \code{\link{optim}}}
  \item{fit.method}{an integer between 0 and 13 determine the fitting routine (i.e. weighting of the squared residuals in the LSE). Values 0 to 6 correspond with the pure spatial version (see \code{\link{fit.variogram}}). See the details section for the meaning of the other values (partly experimental).}
  \item{stAni}{The spatio-temporal anisotropy that is used in the weighting. Might be missing if the desired spatio-temporal variogram model does contain a spatio-temporal anisotropy parameter (this might cause bad convergence behaviour). The default is \code{NA} and will be understood as identity (1 temporal unit = 1 spatial unit). As this only in very few cases a valid assumption, a warning is issued.}
  \item{wles}{Should be missing; only for backwards compatibility, \code{wles = TRUE} corresponds to \code{fit.method = 1} and \code{wles = FALSE} corresponds to \code{fit.method = 6}.}
}

\details{
The following list summarizes the meaning of the \code{fit.method} argument which is essential a weighting of the squared residuals in the least-squares estimation. Please note, that weights based on the models gamma value might fail to converge properly due to the dependence of weights on the variogram estimate:
\describe{
  \item{\code{fit.method = 0}}{no fitting, however the MSE between the provided variogram model and sample variogram surface is calculated.}
  \item{\code{fit.method = 1}}{Number of pairs in the spatio-temporal bin: \eqn{N_j}{N[j]}}
  \item{\code{fit.method = 2}}{Number of pairs in the spatio-temporal bin divided by the square of the current variogram model's value: \eqn{N_j/\gamma(h_j, u_j)^2}{N[j]/gamma(h[j],u[j])^2}}
  \item{\code{fit.method = 3}}{Same as \code{fit.method = 1} for compatibility with \code{\link{fit.variogram}} but as well evaluated in R.}
  \item{\code{fit.method = 4}}{Same as \code{fit.method = 2} for compatibility with \code{\link{fit.variogram}} but as well evaluated in R.}
  \item{\code{fit.method = 5}}{Reserved for REML for compatibility with \code{\link{fit.variogram}}, not yet implemented.}
  \item{\code{fit.method = 6}}{No weights.}
  \item{\code{fit.method = 7}}{Number of pairs in the spatio-temporal bin divided by the square of the bin's metric distance. If \code{stAni} is not specified, the model's  parameter is used to calculate the metric distance across space and time: \eqn{N_j/(h_j^2 + {\rm stAni}^2\cdot u_j^2)}{N[j]/(h[j]^2+ stAni^2*u[j]^2)}}
  \item{\code{fit.method = 8}}{Number of pairs in the spatio-temporal bin divided by the square of the bin's spatial distance. \eqn{N_j/h_j^2}{N[j]/h[j]^2}. Note that the 0 distances are replaced by the smallest non-zero distances to avoid division by zero.}
  \item{\code{fit.method = 9}}{Number of pairs in the spatio-temporal bin divided by the square of the bin's temporal distance. \eqn{N_j/u_j^2}{N[j]/u[j]^2}. Note that the 0 distances are replaced by the smallest non-zero distances to avoid division by zero.}
  \item{\code{fit.method = 10}}{Reciprocal of the square of the current variogram model's value: \eqn{1/\gamma(h_j,u_j)^2}{1/gamma(h[j],u[j])^2}}
  \item{\code{fit.method = 11}}{Reciprocal of the square of the bin's metric distance. If \code{stAni} is not specified, the model's  parameter is used to calculate the metric distance across space and time: \eqn{1/(h_j^2 + {\rm stAni}^2\cdot u_j^2)}{1/(h[j]^2+ stAni^2*u[j]^2)}}
  \item{\code{fit.method = 12}}{Reciprocal of the square of the bin's spatial distance. \eqn{1/h_j^2}{1/h[j]^2}. Note that the 0 distances are replaced by the smallest non-zero distances to avoid division by zero.}
  \item{\code{fit.method = 13}}{Reciprocal of the square of the bin's temporal distance. \eqn{1/u_j^2}{1/u[j]^2}. Note that the 0 distances are replaced by the smallest non-zero distances to avoid division by zero.}
}
See also Table 4.2 in the gstat manual for the original spatial version.
}

\value{
Returns a spatio-temporal variogram model, as S3 class StVariogramModel. It carries the temporal and spatial unit as attributes \code{"temporal unit"} and \code{"spatial unit"} in order to allow \code{\link{krigeST}} to adjust for different units. The units are obtained from the provided empirical variogram. Further attributes are the optim output \code{"optim.output"} and the always not weighted mean squared error \code{"MSE"}.
}

\author{
Benedikt Graeler
}

\seealso{
\code{\link{fit.variogram}} for the pure spatial case. \code{\link{extractParNames}} helps to understand the parameter structure of spatio-temporal variogram models.
}
\examples{
# separable model: spatial and temporal sill will be ignored
# and kept constant at 1-nugget respectively. A joint sill is used.
\dontrun{
separableModel <- vgmST("separable", 
                        method = "Nelder-Mead", # no lower & upper needed
                        space=vgm(0.9,"Exp", 123, 0.1),
                        time =vgm(0.9,"Exp", 2.9, 0.1),
                        sill=100)

data(vv)
separableModel <- fit.StVariogram(vv, separableModel,
                                  method="L-BFGS-B",
                                  lower=c(10,0,0.01,0,1),
                                  upper=c(500,1,20,1,200))
plot(vv, separableModel)
} # dontrun
}

\keyword{models}
