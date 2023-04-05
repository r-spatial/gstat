% $Id: variogram.Rd,v 1.23 2009-11-02 21:33:17 edzer Exp $
\name{variogram}
\alias{variogram}
\alias{variogram.gstat}
\alias{variogram.formula}
\alias{variogram.default}
\alias{print.gstatVariogram}
\alias{print.variogramCloud}
\alias{as.data.frame.variogramCloud}
\title{ Calculate Sample or Residual Variogram or Variogram Cloud }
\description{
Calculates the sample variogram from data, or in case of a linear model
is given, for the residuals, with options for directional, robust,
and pooled variogram, and for irregular distance intervals.

In case spatio-temporal data is provided, the function \code{\link{variogramST}}
is called with a different set of parameters.
}
\usage{
\method{variogram}{gstat}(object, ...)
\method{variogram}{formula}(object, locations = coordinates(data), data, ...)
\method{variogram}{default}(object, locations, X, cutoff, width = cutoff/15,
	alpha = 0, beta = 0, tol.hor = 90/length(alpha), tol.ver =
	90/length(beta), cressie = FALSE, dX = numeric(0), boundaries =
	numeric(0), cloud = FALSE, trend.beta = NULL, debug.level = 1,
	cross = TRUE, grid, map = FALSE, g = NULL, ..., projected = TRUE, 
	lambda = 1.0, verbose = FALSE, covariogram = FALSE, PR = FALSE, 
	pseudo = -1)
\method{print}{gstatVariogram}(x, ...)
\method{print}{variogramCloud}(x, ...)
}
\arguments{
\item{object}{object of class \code{gstat}; in this form, direct
and cross (residual) variograms are calculated for all variables and
variable pairs defined in \code{object}; in case of \code{variogram.formula},
formula defining the response vector and (possible) 
regressors, in case of absence of regressors, use e.g. \code{z~1};
in case of \code{variogram.default}: list with for each variable 
the vector with responses (should not be called directly)
}
\item{data}{data frame where the names in formula are to be found}
\item{locations}{ spatial data locations.  For variogram.formula: a
formula with only the coordinate variables in the right hand (explanatory
variable) side e.g. \code{~x+y}; see examples.

For variogram.default: list with coordinate matrices, each with the number
of rows matching that of corresponding vectors in y; the number of columns
should match the number of spatial dimensions spanned by the data (1
(x), 2 (x,y) or 3 (x,y,z)).  }
\item{...}{any other arguments that will be passed to
\link{variogram.default} (ignored)}
\item{X}{ (optional) list with for each variable the matrix with
regressors/covariates; the number of rows should match that of the
correspoding element in y, the number of columns equals the number of
regressors (including intercept) }
\item{cutoff}{ spatial separation distance up to which point pairs
are included in semivariance estimates; as a default, the
length of the diagonal of the box spanning the data is divided by three. }
\item{width}{ the width of subsequent distance intervals into which
data point pairs are grouped for semivariance estimates }
\item{alpha}{ direction in  plane (x,y), in positive degrees clockwise
from positive y (North): alpha=0 for direction North (increasing y),
alpha=90 for direction East (increasing x); optional a vector of
directions in (x,y) }
\item{beta}{ direction in z, in positive degrees up from the (x,y) plane; }
optional a vector of directions
\item{tol.hor}{ horizontal tolerance angle in degrees }
\item{tol.ver}{ vertical tolerance angle in degrees }
\item{cressie}{ logical; if TRUE, use Cressie''s robust variogram estimate;
if FALSE use the classical method of moments variogram estimate }
\item{dX}{ include a pair of data points $y(s_1),y(s_2)$ taken at
locations $s_1$ and $s_2$ for sample variogram calculation only when
$||x(s_1)-x(s_2)|| < dX$ with and $x(s_i)$ the vector with regressors at
location $s_i$, and $||.||$ the 2-norm.  This allows pooled estimation of
within-strata variograms (use a factor variable as regressor, and dX=0.5),
or variograms of (near-)replicates in a linear model (addressing point
pairs having similar values for regressors variables) }
\item{boundaries}{ numerical vector with distance interval upper boundaries; 
values should be strictly increasing }
\item{cloud}{ logical; if TRUE, calculate the semivariogram cloud }
\item{trend.beta}{vector with trend coefficients, in case they are
known. By default, trend coefficients are estimated from the data.}
\item{debug.level}{ integer; set gstat internal debug level }
\item{cross}{ logical or character; if FALSE, no cross variograms are computed
when object is of class \code{gstat} and has more than one variable; if
TRUE, all direct and cross variograms are computed; if
equal to "ST", direct and cross variograms are computed for all pairs 
involving the first (non-time lagged) variable; if equal to "ONLY",
only cross variograms are computed (no direct variograms). }
\item{formula}{formula, specifying the dependent variable and possible covariates}
\item{x}{ object of class \code{variogram} or \code{variogramCloud}
to be printed}
\item{grid}{ grid parameters, if data are gridded (not to be called
directly; this is filled automatically) }
\item{map}{ logical; if TRUE, and \code{cutoff} and \code{width}
are given, a variogram map is returned. This requires package
sp. Alternatively, a map can be passed, of class SpatialDataFrameGrid
(see sp docs) }
\item{g}{ NULL or object of class gstat; may be used to pass settable
parameters and/or variograms; see example }
\item{projected}{logical; if FALSE, data are assumed to be unprojected,
meaning decimal longitude/latitude. For projected data, Euclidian
distances are computed, for unprojected great circle distances
(km). In \code{variogram.formula} or \code{variogram.gstat}, for data
deriving from class Spatial, projection is detected automatically using
\code{is.projected}}
\item{lambda}{test feature; not working (yet)}
\item{verbose}{logical; print some progress indication}
\item{pseudo}{ integer; use pseudo cross variogram for computing
time-lagged spatial variograms? -1: find out from coordinates -- if they
are equal then yes, else no; 0: no; 1: yes. }
\item{covariogram}{logical; compute covariogram instead of variogram?}
\item{PR}{logical; compute pairwise relative variogram (does NOT check
whether variable is strictly positive)}
}

\value{ If map is TRUE (or a map is passed), a grid map is returned
containing the (cross) variogram map(s). See package sp.

In other cases, an object of class "gstatVariogram" with the 
following fields:
\item{np}{the number of point pairs for this estimate; 
in case of a \code{variogramCloud} see below}
\item{dist}{the average distance of all point pairs considered
for this estimate}
\item{gamma}{the actual sample variogram estimate}
\item{dir.hor}{the horizontal direction}
\item{dir.ver}{the vertical direction}
\item{id}{the combined id pair}

If cloud is TRUE: an object of class \code{variogramCloud}, with the field
\code{np} encoding the numbers of the point pair that contributed to a
variogram cloud estimate, as follows. The first point is found by 1 + the
integer division of np by the \code{.BigInt} attribute of the returned
object, the second point by 1 + the remainder of that division. 
\link{as.data.frame.variogramCloud} returns no \code{np} field,
but does the decoding into:
\item{left}{for variogramCloud: data id (row number) of one of 
the data pair}
\item{right}{for variogramCloud: data id (row number) of the other 
data in the pair}

In case of a spatio-temporal variogram is sought see \code{\link{variogramST}} 
for details.
}
\note{ \code{variogram.default} should not be called by users directly,
as it makes many assumptions about the organization of the data, that
are not fully documented (but of course, can be understood from reading
the source code of the other \code{variogram} methods) 

Successfully setting \code{gridded() <- TRUE} may trigger a branch that
will fail unless dx and dy are identical, and not merely similar
to within machine epsilon.  }

\references{ 
Cressie, N.A.C., 1993, Statistics for Spatial Data, Wiley.

Cressie, N., C. Wikle, 2011, Statistics for Spatio-temporal Data, Wiley.

\url{http://www.gstat.org/}

Pebesma, E.J., 2004. Multivariable geostatistics in S: the gstat package.
Computers and Geosciences, 30: 683-691.
}
\author{ Edzer Pebesma }
\note{
\code{variogram.line} is DEPRECATED; it is and was never meant as a variogram
method, but works automatically as such by the R dispatch system. Use
variogramLine instead.
}
\seealso{
\link{print.gstatVariogram},
\link{plot.gstatVariogram},
\link{plot.variogramCloud};
for variogram models: \link{vgm},
to fit a variogram model to a sample variogram: 
\link{fit.variogram}
\code{\link{variogramST}} for details on the spatio-temporal sample variogram.
}
\examples{
library(sp)
data(meuse)
# no trend:
coordinates(meuse) = ~x+y
variogram(log(zinc)~1, meuse)
# residual variogram w.r.t. a linear trend:
variogram(log(zinc)~x+y, meuse)
# directional variogram:
variogram(log(zinc)~x+y, meuse, alpha=c(0,45,90,135))
variogram(log(zinc)~1, meuse, width=90, cutoff=1300)

# GLS residual variogram:
v = variogram(log(zinc)~x+y, meuse)
v.fit = fit.variogram(v, vgm(1, "Sph", 700, 1))
v.fit
set = list(gls=1)
v
g = gstat(NULL, "log-zinc", log(zinc)~x+y, meuse, model=v.fit, set = set)
variogram(g)

if (require(sf)) {
  proj4string(meuse) = CRS("+init=epsg:28992")
  meuse.ll = sf::st_transform(sf::st_as_sf(meuse), sf::st_crs("+proj=longlat +datum=WGS84"))
# variogram of unprojected data, using great-circle distances, returning km as units
  print(variogram(log(zinc) ~ 1, meuse.ll))
}

}
\keyword{models}
