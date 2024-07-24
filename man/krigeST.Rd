% $Id: krige.Rd,v 1.25 2010-01-12 12:15:48 edzer Exp $
\name{krigeST}
\docType{methods}
\alias{krigeST}
\alias{krige,formula,ST-method} % calling krigeST
\alias{krigeSTTg}
% \alias{idw,formula,ST-method}

\title{ Ordinary global Spatio-Temporal Kriging}
\description{
Function for ordinary global and local and trans Gaussian spatio-temporal kriging on point support
}
\usage{
krigeST(formula, data, newdata, modelList, beta, y, ...,
        nmax = Inf, stAni = NULL,
        computeVar = FALSE,	fullCovariance = FALSE,
        bufferNmax=2, progress=TRUE)
krigeSTTg(formula, data, newdata, modelList, y, nmax=Inf, stAni=NULL,
                      bufferNmax=2, progress=TRUE, lambda = 0)
}
\arguments{
 \item{formula}{ formula that defines the dependent variable as a linear
  model of independent variables; suppose the dependent variable has name
  \code{z}, for ordinary and simple kriging use the formula \code{z~1};
  for simple kriging also define \code{beta} (see below); for universal
  kriging, suppose \code{z} is linearly dependent on \code{x} and \code{y},
  use the formula \code{z~x+y}}
 \item{data}{ ST object: should contain the dependent variable and independent
  variables.}
 \item{newdata}{ ST object with prediction/simulation locations in
 space and time; should contain attribute columns with the independent variables
 (if present).}
 \item{modelList}{ object of class \code{StVariogramModel}, created by
 \code{\link{vgmST}} - see below or the function \code{\link{vgmAreaST}} for area-to-point kriging. For the general kriging case: a list with named elements: \code{space}, \code{time}
 and/or \code{joint} depending 
  on the spatio-temporal covariance family, and an entry \code{stModel}. 
  Currently 
  implemented families that may be used for \code{stModel} 
  are \code{separable}, \code{productSum}, \code{metric}, 
  \code{sumMetric} and \code{simpleSumMetric}. 
  See the examples section in \code{\link{fit.StVariogram}} or 
  \code{\link{variogramSurface}} for details on how to define spatio-temporal 
  covariance models. \code{krigeST} will look for a "temporal unit" attribute
  in the provided modelList in order to adjust the temporal scales.}
 \item{y}{matrix; to krige multiple fields in a single step, pass data 
  as columns of matrix \code{y}. This will ignore the value of the 
  response in \code{formula}.}
 \item{beta}{The (known) mean for simple kriging.}
 \item{nmax}{The maximum number of neighbouring locations for a spatio-temporal local neighbourhood}
 \item{stAni}{a spatio-temporal anisotropy scaling assuming a metric spatio-temporal space. Used only for the selection of the closest neighbours. This scaling needs only to be provided in case the model does not have a stAni parameter, or if a different one should be used for the neighbourhood selection. Mind the correct spatial unit. Currently, no coordinate conversion is made for the neighbourhood selection (i.e. Lat and Lon require a spatio-temporal anisotropy scaling in degrees per second).}
 \item{\dots}{ further arguments used for instance to pass the model into vgmAreaST for area-to-point kriging}
 \item{computeVar}{logical; if TRUE, prediction variances will be returned}
 \item{fullCovariance}{logical; if FALSE a vector with prediction 
  variances will be returned, if TRUE the full covariance matrix of 
  all predictions will be returned}
  \item{bufferNmax}{factor with which nmax is multiplied for an extended search radius (default=2). Set to 1 for no extension of the search radius.}
  \item{progress}{whether a progress bar shall be printed for local spatio-temporal kriging; default=TRUE}
  \item{lambda}{The value of lambda used in the box-cox transformation.}
}
\details{
 Function \code{krigeST} is a R implementation of the kriging function from 
 \link{gstat} using spatio-temporal covariance models following the 
 implementation of \code{\link{krige0}}. Function \code{krigeST} offers some 
 particular methods for ordinary spatio-temporal (ST) kriging. In particular,
 it does not support block kriging or kriging in a distance-based
 neighbourhood, and does not provide simulation.

 If \code{data} is of class \code{sftime}, then \code{newdata} MUST be
 of class \code{stars} or \code{sftime}, i.e. mixing form old-style
 classes (package spacetime) and new-style classes (sf, stars, sftime)
 is not supported.
}

\value{
 An object of the same class as \code{newdata} (deriving from 
 \code{\link[spacetime:ST-class]{ST}}). Attributes columns contain prediction and prediction 
 variance.
}

\references{ 
Benedikt Graeler, Edzer Pebesma, Gerard Heuvelink.  Spatio-Temporal
Geostatistics using gstat.  The R Journal 8(1), 204--218.
\url{https://journal.r-project.org/archive/2016/RJ-2016-014/index.html}

N.A.C. Cressie, 1993, Statistics for Spatial Data,
Wiley. 

\url{http://www.gstat.org/}

Pebesma, E.J., 2004. Multivariable geostatistics in S: the gstat package.
Computers and Geosciences, 30: 683-691.  }
\author{ Edzer Pebesma, Benedikt Graeler }

\seealso{ \code{\link{krige0}}, \code{\link{gstat}}, \code{\link[gstat]{predict}}, \code{\link{krigeTg}} }

\examples{
library(sp)
library(spacetime)
sumMetricVgm <- vgmST("sumMetric",
                      space = vgm( 4.4, "Lin", 196.6,  3),
                      time  = vgm( 2.2, "Lin",   1.1,  2),
                      joint = vgm(34.6, "Exp", 136.6, 12),
                      stAni = 51.7)

data(air)
suppressWarnings(proj4string(stations) <- CRS(proj4string(stations)))
rural = STFDF(stations, dates, data.frame(PM10 = as.vector(air)))

rr <- rural[,"2005-06-01/2005-06-03"]
rr <- as(rr,"STSDF")

x1 <- seq(from=6,to=15,by=1)
x2 <- seq(from=48,to=55,by=1)

DE_gridded <- SpatialPoints(cbind(rep(x1,length(x2)), rep(x2,each=length(x1))), 
                            proj4string=CRS(proj4string(rr@sp)))
gridded(DE_gridded) <- TRUE
DE_pred <- STF(sp=as(DE_gridded,"SpatialPoints"), time=rr@time)
DE_kriged <- krigeST(PM10~1, data=rr, newdata=DE_pred,
                     modelList=sumMetricVgm)
gridded(DE_kriged@sp) <- TRUE
stplot(DE_kriged)
}
\keyword{ models }
