% $Id: krigeTg.Rd,v 1.10 2009-08-17 14:38:28 edzer Exp $
\name{krigeTg}
\docType{methods}
\alias{krigeTg}
\title{ TransGaussian kriging using Box-Cox transforms }
\description{
TransGaussian (ordinary) kriging function using Box-Cox transforms  
}
\usage{
krigeTg(formula, locations, newdata, model = NULL, ...,
	nmax = Inf, nmin = 0, maxdist = Inf, block = numeric(0),
	nsim = 0, na.action = na.pass, debug.level = 1,
	lambda = 1.0)
}
\arguments{
 \item{formula}{ formula that defines the dependent variable as a linear
  model of independent variables; suppose the dependent variable has name
  \code{z}, for ordinary and use a formula like \code{z~1}; the dependent
  variable should be NOT transformed.
  }
 \item{locations}{ object of class \code{Spatial}, with observations }
 \item{newdata}{ Spatial object with prediction/simulation locations; 
  the coordinates should have names as defined in \code{locations} }
 \item{model}{ variogram model of the TRANSFORMED dependent variable, 
  see \link{vgm}, or \link{fit.variogram}}
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
 \item{block}{ does not function correctly, afaik }
 \item{nsim}{ does not function correctly, afaik }
 \item{na.action}{ function determining what should be done with missing
  values in 'newdata'.  The default is to predict 'NA'.  Missing values 
  in coordinates and predictors are both dealt with. }
 \item{lambda}{value for the Box-Cox transform}
 \item{debug.level}{debug level, passed to \link[gstat]{predict}; use -1 
 to see progress in percentage, and 0 to suppress all printed information }
 \item{\dots}{ other arguments that will be passed to \link{gstat}}
}
\details{
Function \code{krigeTg} uses transGaussian kriging as explained in
\url{https://www.math.umd.edu/~bnk/bak/Splus/kriging.html}. 

As it uses the R/gstat krige function to derive everything, it needs in
addition to ordinary kriging on the transformed scale a simple kriging
step to find m from the difference between the OK and SK prediction
variance, and a kriging/BLUE estimation step to obtain the estimate
of \eqn{\mu}{mu}.

For further details, see \link{krige} and \link[gstat]{predict}.
}

\value{
an SpatialPointsDataFrame object containing the fields:
\code{m} for the m (Lagrange) parameter for each location;
\code{var1SK.pred} the \eqn{c_0 C^{-1}}{c0 Cinv} correction obtained by
\code{muhat} for the mean estimate at each location;
\code{var1SK.var} the simple kriging variance;
\code{var1.pred} the OK prediction on the transformed scale;
\code{var1.var} the OK kriging variance on the transformed scale;
\code{var1TG.pred} the transGaussian kriging predictor;
\code{var1TG.var} the transGaussian kriging variance, obtained by
\eqn{\phi'(\hat{\mu},\lambda)^2 \sigma^2_{OK}}{phi'(muhat, lambda)^2 * var1.var}

}

\references{ N.A.C. Cressie, 1993, Statistics for Spatial Data,
Wiley. 

\url{http://www.gstat.org/}

}

\author{ Edzer Pebesma }

\seealso{ \link{gstat}, \link[gstat]{predict} }

\examples{
\donttest{
library(sp)
data(meuse)
coordinates(meuse) = ~x+y
data(meuse.grid)
gridded(meuse.grid) = ~x+y
v = vgm(1, "Exp", 300)
x1 = krigeTg(zinc~1,meuse,meuse.grid,v, lambda=1) # no transform
x2 = krige(zinc~1,meuse,meuse.grid,v)
summary(x2$var1.var-x1$var1TG.var)
summary(x2$var1.pred-x1$var1TG.pred)
lambda = -0.25
m = fit.variogram(variogram((zinc^lambda-1)/lambda ~ 1,meuse), vgm(1, "Exp", 300))
x = krigeTg(zinc~1,meuse,meuse.grid,m,lambda=-.25)
spplot(x["var1TG.pred"], col.regions=bpy.colors())
summary(meuse$zinc)
summary(x$var1TG.pred)
}
}
\keyword{ models }
