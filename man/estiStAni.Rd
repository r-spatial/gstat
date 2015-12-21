\name{estiStAni}
\alias{estiStAni}

\title{
Estimation of the spatio-temporal anisotropy
}
\description{
Estimation of the spatio-temporal anisotropy without an underlying spatio-temporal model. Different methods are implemented using a linear model to predict the temporal gamma values or the ratio of the ranges of a spatial and temporal variogram model or a spatial variogram model to predict the temporal gamma values or the spatio-temporal anisotropy value as used in a metric spatio-temporal variogram. 
}
\usage{
estiStAni(empVgm, interval, method = "linear", spatialVgm, 
          temporalVgm, s.range=NA, t.range=NA)
}

\arguments{
  \item{empVgm}{
  An empirical spatio-temporal variogram.
}
  \item{interval}{
  A search interval for the optimisation of the spatio-temporal anisotropy parameter
}
  \item{method}{
  A character string determining the method to be used (one of \code{linear}, \code{range}, \code{vgm} or \code{metric}, see below for details)
}
  \item{spatialVgm}{
  A spatial variogram definition from the call to \code{\link{vgm}}. The model is optimised based on the pure spatial values in \code{empVgm}.
}
  \item{temporalVgm}{
  A temporal variogram definition from the call to \code{\link{vgm}}. The model is optimised based on the pure temporal values in \code{empVgm}.}
  \item{s.range}{
  A spatial cutoff value applied to the empirical variogram \code{empVgm}.  
  }
  \item{t.range}{
  A temporal cutoff value applied to the empirical variogram \code{empVgm}.  
  }
}
\details{
\describe{
  \item{linear}{
  A linear model is fitted to the pure spatial gamma values based on the spatial distances. An optimal scaling is searched to stretch the temporal distances such that the linear model explains best the pure temporal gamma values. This assumes (on average) a linear relationship between distance and gamma, hence it is advisable to use only those pairs of pure spatial (pure temporal) distance and gamma value that show a considerable increase (i.e. drop all values beyond the range by setting values for \code{s.range} and \code{t.range}).
  }
  \item{range}{
  A spatial and temporal variogram model is fitted to the pure spatial and temporal gamma values respectively. The spatio-temporal anisotropy estimate is the ratio of the spatial range over the temporal range.
  }
  \item{vgm}{
  A spatial variogram model is fitted to the pure spatial gamma values. An optimal scaling is used to stretch the temporal distances such that the spatial variogram model explains best the pure temporal gamma values.
  }
  \item{metric}{
  A metric spatio-temporal variogram model is fitted with \code{joint} component according to the defined spatial variogram \code{spatialVgm}. The starting value of \code{stAni} is the mean of the \code{interval} parameter (see \code{\link{vgmST}} for the metric variogram definition). The spatio-temporal anisotropy as estimated in the spatio-temporal variogram is returned. Note that the parameter \code{interval} is only used to set the starting value. Hence, the estimate might exceed the given interval.
  }
}
}
\value{A scalar representing the spatio-temporal anisotropy estimate.}

\note{Different methods might lead to very different estimates. All but the \code{linear} approach are sensitive to the variogram model selection.}

\author{Benedikt Graeler}

\examples{
data(vv)

estiStAni(vv, c(10, 150))
estiStAni(vv, c(10, 150), "vgm", vgm(80, "Sph", 120, 20))
}
