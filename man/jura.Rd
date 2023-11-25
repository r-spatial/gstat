% $Id: jura.Rd,v 1.3 2007-03-08 09:35:49 edzer Exp $
\name{jura}
\alias{jura}
\alias{prediction.dat}
\alias{validation.dat}
\alias{transect.dat}
\alias{juragrid.dat}
\alias{jura.grid}
\alias{jura.pred}
\alias{jura.val}
\title{Jura data set}
\description{
The jura data set from Pierre Goovaerts' book (see references
below). It contains four \code{data.frame}s: prediction.dat, validation.dat
and transect.dat and juragrid.dat, and three \code{data.frame}s with
consistently coded land use and rock type factors, as well as geographic
coordinates. The examples below show how to transform these into
spatial (sp) objects in a local coordinate system and in geographic
coordinates, and how to transform to metric coordinate reference
systems.
}
\format{
  The \code{data.frames} prediction.dat and validation.dat contain the following fields:
  \describe{
	\item{Xloc}{ X coordinate, local grid km }
	\item{Yloc}{ Y coordinate, local grid km }
	\item{Landuse}{ see book and below }
	\item{Rock}{ see book and below }
	\item{Cd}{ mg cadmium \eqn{\mbox{kg}^{-1}}{kg^-1} topsoil }
	\item{Co}{ mg cobalt \eqn{\mbox{kg}^{-1}}{kg^-1} topsoil }
	\item{Cr}{ mg chromium \eqn{\mbox{kg}^{-1}}{kg^-1} topsoil }
	\item{Cu}{ mg copper \eqn{\mbox{kg}^{-1}}{kg^-1} topsoil }
	\item{Ni}{ mg nickel \eqn{\mbox{kg}^{-1}}{kg^-1} topsoil }
	\item{Pb}{ mg lead \eqn{\mbox{kg}^{-1}}{kg^-1} topsoil }
	\item{Zn}{ mg zinc \eqn{\mbox{kg}^{-1}}{kg^-1} topsoil }
  }
  The \code{data.frame} juragrid.dat only has the first four fields.
  In addition the \code{data.frame}s jura.pred, jura.val and jura.grid also
  have inserted third and fourth fields giving geographic coordinates:
  \describe{
  	\item{long}{ Longitude, WGS84 datum }
	\item{lat}{ Latitude, WGS84 datum }
  } 
} 
\usage{
data(jura)
}
\author{ Data preparation by David Rossiter (dgr2@cornell.edu) 
and Edzer Pebesma; georeferencing by David Rossiter }
\references{ 
Goovaerts, P. 1997. Geostatistics for Natural Resources Evaluation. Oxford
Univ. Press, New-York, 483 p. Appendix C describes (and gives) the Jura
data set.

Atteia, O., Dubois, J.-P., Webster, R., 1994, Geostatistical analysis of
soil contamination in the Swiss Jura: Environmental Pollution 86, 315-327

Webster, R., Atteia, O., Dubois, J.-P., 1994, Coregionalization of trace
metals in the soil in the Swiss Jura: European Journal of Soil Science
45, 205-218

}
\note{
The points data sets were obtained from
http://home.comcast.net/~pgoovaerts/book.html, which seems to be
no longer available; the grid data were kindly provided by Pierre
Goovaerts.

The following codes were used to convert \code{prediction.dat}
and \code{validation.dat} to \code{jura.pred} and \code{jura.val}
(see examples below):

Rock Types: 1: Argovian, 2: Kimmeridgian, 3: Sequanian, 4: Portlandian,
5: Quaternary.

Land uses: 1: Forest, 2: Pasture (Weide(land), Wiese, Grasland),
3: Meadow (Wiese, Flur, Matte, Anger), 4: Tillage (Ackerland,
bestelltes Land)

Points 22 and 100 in the validation set
(\code{validation.dat[c(22,100),]}) seem not to lie exactly on the
grid originally intended, but are kept as such to be consistent with
the book.

Georeferencing was based on two control points in the Swiss grid system
shown as Figure 1 of Atteia et al. (see above) and further points digitized
on the tentatively georeferenced scanned map. RMSE 2.4 m. Location of points
in the field was less precise.
}

\keyword{datasets}
\examples{
data(jura)
summary(prediction.dat)
summary(validation.dat)
summary(transect.dat)
summary(juragrid.dat)

# the following commands were used to create objects with factors instead
# of the integer codes for Landuse and Rock:
\dontrun{
  jura.pred = prediction.dat
  jura.val = validation.dat
  jura.grid = juragrid.dat

  jura.pred$Landuse = factor(prediction.dat$Landuse, 
	labels=levels(juragrid.dat$Landuse))
  jura.pred$Rock = factor(prediction.dat$Rock, 
	labels=levels(juragrid.dat$Rock))
  jura.val$Landuse = factor(validation.dat$Landuse, 
	labels=levels(juragrid.dat$Landuse))
  jura.val$Rock = factor(validation.dat$Rock, 
	labels=levels(juragrid.dat$Rock))
}

# the following commands convert data.frame objects into spatial (sp) objects
#   in the local grid:
require(sp)
coordinates(jura.pred) = ~Xloc+Yloc
coordinates(jura.val) = ~Xloc+Yloc
coordinates(jura.grid) = ~Xloc+Yloc
gridded(jura.grid) = TRUE

# the following commands convert the data.frame objects into spatial (sp) objects
#   in WGS84 geographic coordinates
# example is given only for jura.pred, do the same for jura.val and jura.grid
# EPSG codes can be found by searching make_EPSG()
jura.pred <- as.data.frame(jura.pred)
coordinates(jura.pred) = ~ long + lat
proj4string(jura.pred) = CRS("+init=epsg:4326")
}
