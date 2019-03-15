% $Id: krige.cv.Rd,v 1.18 2009-10-30 16:11:21 edzer Exp $
\name{krige.cv}
\docType{methods}
\alias{krige.cv}
\alias{krige.cv.spatial}
\alias{krige.cv.locations}
\alias{gstat.cv}
\alias{krige.cv,formula,formula-method}
\alias{krige.cv,formula,Spatial-method}
\alias{krige.cv,formula,sf-method}
\title{ (co)kriging cross validation, n-fold or leave-one-out }
\description{ Cross validation functions for simple, ordinary or universal
point (co)kriging, kriging in a local neighbourhood.  }
\usage{
gstat.cv(object, nfold, remove.all = FALSE, verbose = interactive(), 
	all.residuals = FALSE, ...)
krige.cv(formula, locations, ...)
krige.cv.locations(formula, locations, data, model = NULL, ..., beta = NULL, 
	nmax = Inf, nmin = 0, maxdist = Inf, nfold = nrow(data), 
	verbose = interactive(), debug.level = 0)
krige.cv.spatial(formula, locations, model = NULL, ..., beta = NULL, 
	nmax = Inf, nmin = 0, maxdist = Inf, nfold = nrow(locations), 
	verbose = interactive(), debug.level = 0)
}
\arguments{
  \item{object}{ object of class gstat; see function \link{gstat}}
  \item{nfold}{ integer; if larger than 1, then apply n-fold cross validation; 
  if \code{nfold} equals \code{nrow(data)} (the default), apply leave-one-out cross 
  validation; if set to e.g. 5, five-fold cross validation is done. To specify the
  folds, pass an integer vector of length \code{nrow(data)} with fold indexes. }
  \item{remove.all}{ logical; if TRUE, remove observations at cross validation
  locations not only for the first, but for all subsequent variables as well }
  \item{verbose}{ logical; if FALSE, progress bar is suppressed }
  \item{all.residuals}{ logical; if TRUE, residuals for all variables are
  returned instead of for the first variable only}
  \item{\dots}{ other arguments that will be passed to \link[gstat]{predict}
  in case of \code{gstat.cv}, or to \link{gstat} in case of \code{krige.cv}}
  \item{formula}{ formula that defines the dependent variable as a linear
  model of independent variables; suppose the dependent variable has name
  \code{z}, for ordinary and simple kriging use the formula \code{z~1};
  for simple kriging also define \code{beta} (see below); for universal
  kriging, suppose \code{z} is linearly dependent on \code{x} and \code{y},
  use the formula \code{z~x+y}}
  \item{locations}{ data object deriving from class \code{Spatial} or \code{sf} }
  \item{data}{ data frame (deprecated); should contain the dependent variable, independent
  variables, and coordinates; only to be provided if \code{locations} is a formula}
  \item{model}{ variogram model of dependent variable (or its residuals), 
  defined by a call to \link{vgm} or \link{fit.variogram}}
  \item{beta}{ only for simple kriging (and simulation based on simple
  kriging); vector with the trend coefficients (including intercept);
  if no independent variables are defined the model only contains an
  intercept and this should be the simple kriging mean }
  \item{nmax}{ for local kriging: the number of nearest observations that
  should be used for a kriging prediction or simulation, where nearest
  is defined in terms of the space of the spatial locations. By default,
  all observations are used }
  \item{nmin}{ for local kriging: if the number of nearest observations
  within distance \code{maxdist} is less than \code{nmin}, a missing 
  value will be generated; see maxdist }
  \item{maxdist}{ for local kriging: only observations within a distance
  of \code{maxdist} from the prediction location are used for prediction
  or simulation; if combined with \code{nmax}, both criteria apply }
  \item{debug.level}{ print debugging information; 0 suppresses
  debug information }
}

\section{Methods}{
\describe{
\item{formula = "formula", locations = "formula"}{ 
locations specifies which coordinates in \code{data} refer to spatial coordinates
}
\item{formula = "formula", locations = "Spatial"}{ 
Object locations knows about its own spatial locations
}
}}

\details{ Leave-one-out cross validation (LOOCV) visits a data point, and
predicts the value at that location by leaving out the observed value,
and proceeds with the next data point. (The observed value is left
out because kriging would otherwise predict the value itself.) N-fold
cross validation makes a partitions the data set in N parts. For all
observation in a part, predictions are made based on the remaining N-1
parts; this is repeated for each of the N parts. N-fold cross validation
may be faster than LOOCV. }

\value{ data frame containing the coordinates of \code{data} or those
of the first variable in \code{object}, and columns of prediction and
prediction variance of cross validated data points, observed values,
residuals, zscore (residual divided by kriging standard error), and fold.

If \code{all.residuals} is true, a data frame with residuals for all
variables is returned, without coordinates. }

\references{ \url{http://www.gstat.org/} }
\author{ Edzer Pebesma }
\note{ Leave-one-out cross validation seems to be much faster in plain
(stand-alone) gstat, apparently quite a bit of the effort is spent moving 
data around from R to gstat. }

\seealso{ \link{krige}, \link{gstat}, \link[gstat]{predict} }

\examples{
library(sp)
data(meuse)
coordinates(meuse) <- ~x+y
m <- vgm(.59, "Sph", 874, .04)
# five-fold cross validation:
x <- krige.cv(log(zinc)~1, meuse, m, nmax = 40, nfold=5)
bubble(x, "residual", main = "log(zinc): 5-fold CV residuals")

# multivariable; thanks to M. Rufino:
meuse.g <- gstat(id = "zn", formula = log(zinc) ~ 1, data = meuse)
meuse.g <- gstat(meuse.g, "cu", log(copper) ~ 1, meuse)
meuse.g <- gstat(meuse.g, model = vgm(1, "Sph", 900, 1), fill.all = TRUE)
x <- variogram(meuse.g, cutoff = 1000)
meuse.fit = fit.lmc(x, meuse.g)
out = gstat.cv(meuse.fit, nmax = 40, nfold = 5) 
summary(out)
out = gstat.cv(meuse.fit, nmax = 40, nfold = c(rep(1,100), rep(2,55))) 
summary(out)
# mean error, ideally 0:
mean(out$residual)
# MSPE, ideally small
mean(out$residual^2)
# Mean square normalized error, ideally close to 1
mean(out$zscore^2)
# correlation observed and predicted, ideally 1
cor(out$observed, out$observed - out$residual)
# correlation predicted and residual, ideally 0
cor(out$observed - out$residual, out$residual)
}
\keyword{ models }
