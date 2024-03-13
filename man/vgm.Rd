% $Id: vgm.Rd,v 1.16 2008-10-09 14:25:20 edzer Exp $
\name{vgm}
\alias{vgm}
\alias{print.variogramModel}
\alias{plot.variogramModel}
\alias{as.vgm.variomodel}
\title{
Generate, or Add to Variogram Model
}
\description{
Generates a variogram model, or adds to an existing model. 
\code{print.variogramModel} prints the essence of a variogram  model.
}
\usage{
vgm(psill = NA, model, range = NA, nugget, add.to, anis, kappa = 0.5, ..., covtable,
	Err = 0)
\method{print}{variogramModel}(x, ...)
\method{plot}{variogramModel}(x, cutoff, ..., type = 'l')
as.vgm.variomodel(m)
}
\arguments{
\item{psill}{ (partial) sill of the variogram model component, or model: see Details }
\item{model}{ model type, e.g. "Exp", "Sph", "Gau", or "Mat". Can be a character vector of model 
types combined with c(), e.g. c("Exp", "Sph"), in which case the best fitting is returned. 
Calling vgm() without a model argument returns a data.frame with available models. }
\item{range}{ range parameter of the variogram model component; in case of anisotropy: major range }
\item{kappa}{ smoothness parameter for the Matern class of variogram 
models }
\item{nugget}{ nugget component of the variogram (this basically adds a nugget
compontent to the model); if missing, nugget component is omitted }
\item{add.to}{ the variogram model to which we want to add a component (structure) }
\item{anis}{ anisotropy parameters: see notes below }
\item{x}{ a variogram model to print or plot}
\item{...}{ arguments that will be passed to \code{print}, e.g. \code{digits}
(see examples), or to \code{variogramLine} for the plot method }
\item{covtable}{ if model is \code{Tab}, instead of model parameters a one-dimensional 
covariance table can be passed here. See covtable.R in tests directory, and example below. }
\item{Err}{ numeric; if larger than zero, the measurement error variance
component that will not be included to the kriging equations, i.e. kriging
will now smooth the process Y instead of predict the measured Z, 
where Z=Y+e, and Err is the variance of e}
\item{m}{ object of class \code{variomodel}, see \pkg{geoR}}
\item{cutoff}{maximum distance up to which variogram values are computed}
\item{type}{plot type}
}
\value{
If a single model is passed, an object of class \code{variogramModel}
extending \code{data.frame}.

In case a vector ofmodels is passed, an object of class
\code{variogramModelList} which is a list of \code{variogramModel}
objects.

When called without a model argument, a data.frame with available models
is returned, having two columns: short (abbreviated names, to be used
as model argument: "Exp", "Sph" etc) and long (with some description).

as.vgm.variomodel tries to convert an object of class variomodel (geoR)
to vgm.
}
\author{ Edzer Pebesma }
\details{ If only the first argument (\code{psill}) is given a
\code{character} value/vector indicating one or more models, as in \code{vgm("Sph")},
then this taken as a shorthand form of \code{vgm(NA,"Sph",NA,NA)},
i.e. a spherical variogram with nugget and unknown parameter values;
see examples below. Read \link{fit.variogram} to find out how
\code{NA} variogram parameters are given initial values for a fitting
a model, based on the sample variogram. Package \code{automap}
gives further options for automated variogram modelling. }

\note{ 
Geometric anisotropy can be modelled for each individual simple model
by giving two or five anisotropy parameters, two for two-dimensional
and five for three-dimensional data. In any case, the range defined
is the range in the direction of the strongest correlation, or the
major range. Anisotropy parameters define which direction this is (the
main axis), and how much shorter the range is in (the) direction(s)
perpendicular to this main axis.

In two dimensions, two parameters define an anisotropy ellipse, say
\code{anis = c(30, 0.5)}. The first parameter, \code{30}, refers to
the main axis direction: it is the angle for the principal direction
of continuity (measured in degrees, clockwise from positive Y, i.e. North).
The second parameter, \code{0.5}, is the anisotropy ratio, the ratio
of the minor range to the major range (a value between 0 and 1).  So,
in our example, if the range in the major direction (North-East) is 100,
the range in the minor direction (South-East) is 0.5 x 100 = 50.

In three dimensions, five values should be given in the form \code{anis
= c(p,q,r,s,t)}.  Now, $p$ is the angle for the principal direction of
continuity (measured in degrees, clockwise from Y, in direction of X),
$q$ is the dip angle for the principal direction of continuity (measured
in positive degrees up from horizontal), $r$ is the third rotation angle
to rotate the two minor directions around the principal direction defined
by $p$ and $q$. A positive angle acts counter-clockwise while looking
in the principal direction. Anisotropy ratios $s$ and $t$ are the ratios
between the major range and each of the two minor ranges. The anisotropy code
was taken from GSLIB. Note that in \url{http://www.gslib.com/sec_gb.html}
it is reported that this code has a bug. Quoting from this
site: ``The third angle in all GSLIB programs operates in the opposite
direction than specified in the GSLIB book.  Explanation - The books
says (pp27) the angle is measured clockwise when looking toward
the origin (from the postive principal direction), but it should be
counter-clockwise. This is a documentation error. Although rarely used,
the correct specification of the third angle is critical if used.''

(Note that \code{anis = c(p,s)} is equivalent to \code{anis = c(p,0,0,s,1)}.)

The implementation in gstat for 2D and 3D anisotropy was taken from
the gslib (probably 1992) code.  I have seen a paper where it is argued
that the 3D anisotropy code implemented in gslib (and so in gstat)
is in error, but I have not corrected anything afterwards.
}
\references{ \url{http://www.gstat.org/} 

Pebesma, E.J., 2004. Multivariable geostatistics in S: the gstat package.
Computers and Geosciences, 30: 683-691.

Deutsch, C.V. and Journel, A.G., 1998.  GSLIB: Geostatistical software
library and user's guide, second edition, Oxford University Press.

For the validity of variogram models on the sphere, see Huang,
Chunfeng, Haimeng Zhang, and Scott M. Robeson. On the validity
of commonly used covariance and variogram functions on the
sphere. Mathematical Geosciences 43.6 (2011): 721-733.

}
\seealso{
\link{show.vgms} to view the available models,
\link{fit.variogram},
\link{variogramLine},
\link{variogram} for the sample variogram.
}
\examples{
vgm()
vgm("Sph")
vgm(NA, "Sph", NA, NA)
vgm(, "Sph") # "Sph" is second argument: NO nugget in this case
vgm(10, "Exp", 300)
x <- vgm(10, "Exp", 300)
vgm(10, "Nug", 0)
vgm(10, "Exp", 300, 4.5)
vgm(10, "Mat", 300, 4.5, kappa = 0.7)
vgm( 5, "Exp", 300, add.to = vgm(5, "Exp", 60, nugget = 2.5))
vgm(10, "Exp", 300, anis = c(30, 0.5))
vgm(10, "Exp", 300, anis = c(30, 10, 0, 0.5, 0.3))
# Matern variogram model:
vgm(1, "Mat", 1, kappa=.3)
x <- vgm(0.39527463, "Sph", 953.8942, nugget = 0.06105141)
x
print(x, digits = 3);
# to see all components, do
print.data.frame(x)
vv=vgm(model = "Tab",  covtable = 
	variogramLine(vgm(1, "Sph", 1), 1, n=1e4, min = 0, covariance = TRUE))
vgm(c("Mat", "Sph"))
vgm(, c("Mat", "Sph")) # no nugget
} 
\keyword{models}
