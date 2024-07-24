\name{variogramST}
\alias{variogramST}
\title{Calculate Spatio-Temporal Sample Variogram}
\description{
Calculates the sample variogram from spatio-temporal data.
}
\usage{
variogramST(formula, locations, data, ..., tlags = 0:15,  cutoff,
	width = cutoff/15, boundaries = seq(0, cutoff, width),
	progress = interactive(), pseudo = TRUE, assumeRegular = FALSE, 
	na.omit = FALSE, cores = 1)
}
\arguments{
\item{formula}{formula, specifying the dependent variable.}
\item{locations}{A STFDF or STSDF containing the variable; kept for 
  compatibility reasons with variogram, either \code{locations} or \code{data} 
  must be provided.}
\item{data}{A \code{\link[spacetime:STFDF-class]{STFDF}}, \code{\link[spacetime:STSDF-class]{STSDF}} or \code{\link[spacetime:STIDF]{STIDF}} containing the variable.}
\item{...}{any other arguments that will be passed to the underlying 
  \code{\link{variogram}} function. In case of using data of type \code{\link[spacetime:STIDF-class]{STIDF}}, the argument \code{tunit} is recommended (and only used in the case of STIDF) to set the temporal unit of the \code{tlags}. Additionally, \code{twindow} can be passed to control the temporal window used for temporal distance calculations. This builds on the property of xts being ordered and only the next \code{twindow} instances are considered. This avoids the need of huge temporal distance matrices. The default uses twice the number as the average difference goes into the temporal cutoff.}
\item{tlags}{ integer; time lags to consider or in case \code{data} is of class \code{\link[spacetime:STIDF-class]{STIDF}} the actual temporal boundaries with time unit given by \code{tunit} otherwise the same unit as \code{\link{diff}} on the index of the time slot will generate is assumed.}
\item{cutoff}{ spatial separation distance up to which point pairs
are included in semivariance estimates; as a default, the
length of the diagonal of the box spanning the data is divided by three. }
\item{width}{ the width of subsequent distance intervals into which
data point pairs are grouped for semivariance estimates, by default the 
  \code{cutoff} is divided into 15 equal lags.}
\item{boundaries}{ numerical vector with distance interval upper boundaries; 
values should be strictly increasing }
\item{progress}{ logical; if TRUE, show text progress bar }
\item{pseudo}{ integer; use pseudo cross variogram for computing
time-lagged spatial variograms? -1: find out from coordinates -- if they
are equal then yes, else no; 0: no; 1: yes. }
\item{assumeRegular}{logical; whether the time series should be assumed regular. 
  The first time step is assumed to be representative for the whole series. Note,
  that temporal lags are considered by index, and no check is made whether pairs
  actually have the desired separating distance.}
\item{na.omit}{shall all \code{NA} values in the spatio-temporal variogram be dropped? In case where complete rows or columns in the variogram consists of \code{NA} only, \code{plot} might produce a distorted picture.}
\item{cores}{number of cores to use in parallel}
}

\value{The spatio-temporal sample variogram contains besides the fields 
  \code{np}, \code{dist} and \code{gamma} the spatio-temporal fields, 
  \code{timelag}, \code{spacelag} and \code{avgDist}, the first of which indicates the time lag 
  used, the second and third different spatial lags. \code{spacelag} is the midpoint in the spatial
  lag intervals as passed by the parameter \code{boundaries}, whereas \code{avgDist} is the average 
  distance between the point pairs found in a distance interval over all temporal lags (i.e. the 
  averages of the values \code{dist} per temporal lag.) To compute variograms for space lag \eqn{h} and
  time lag \eqn{t}, the pseudo cross variogram \eqn{(Z_i(s)-Z_{i+t}(s+h))^2} is averaged over all time 
  lagged observation sets \eqn{Z_i} and \eqn{Z_{i+t}} available (weighted by the number of pairs involved).
}

\references{ 
Cressie, N.A.C., 1993, Statistics for Spatial Data, Wiley.

Cressie, N., C. Wikle, 2011, Statistics for Spatio-temporal Data, Wiley.

\url{http://www.gstat.org/}

Pebesma, E.J., 2004. Multivariable geostatistics in S: the gstat package.
Computers and Geosciences, 30: 683-691.
}
\author{ Edzer Pebesma, Benedikt Graeler }

\seealso{
\code{\link{plot.StVariogram}},
for variogram models: \code{\link{vgmST}},
to fit a spatio-temporal variogram model to a spatio-temporal sample variogram: 
\code{\link{fit.StVariogram}}
}
\examples{
# The following spatio-temporal variogram has been calcualted through
# vv = variogram(PM10~1, r5to10, width=20, cutoff = 200, tlags=0:5)
# in the vignette "st".

data(vv)
str(vv)
plot(vv)
}

\keyword{models}
