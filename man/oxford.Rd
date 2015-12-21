% $Id: oxford.Rd,v 1.5 2006-02-10 19:03:27 edzer Exp $
\name{oxford}
\alias{oxford}
\title{Oxford soil samples}
\description{
Data: 126 soil augerings on a 100 x 100m square grid, with 6 columns
and 21 rows. Grid is oriented with long axis North-north-west to
South-south-east Origin of grid is South-south-east point, 100m outside
grid.

Original data are part of a soil survey carried out by P.A. Burrough in
1967. The survey area is located on the chalk downlands on the Berkshire
Downs in Oxfordshire, UK. Three soil profile units were recognised on
the shallow Rendzina soils; these are Ia - very shallow, grey calcareous
soils less than 40cm deep over chalk;  Ct - shallow to moderately deep,
grey-brown calcareous soils  on calcareous colluvium, and Cr: deep,
moderately acid, red-brown clayey soils.  These soil profile classes
were registered at every augering.

In addition, an independent landscape soil map was made by interpolating
soil boundaries between these soil types, using information from the
changes in landform.  Because the soil varies over short distances,
this field mapping caused some soil borings to receive a different
classification from the classification based on the point data.

Also registered at each auger point were the site elevation (m), the
depth to solid chalk rock (in cm) and the depth to lime in cm. Also,
the percent clay content, the Munsell colour components of VALUE and
CHROMA , and the lime content of the soil (as tested using HCl)  were
recorded for the top two soil layers (0-20cm and 20-40cm).

Samples of topsoil taken as a bulk sample within a circle of radius 2.5m
around each sample point were used for the laboratory determination of Mg
(ppm), OM1 \%, CEC as mequ/100g air dry soil, pH, P as ppm and K (ppm).

}
\format{
This data frame contains the following columns:
\describe{
\item{PROFILE}{profile number}
\item{XCOORD}{x-coordinate, field, non-projected}
\item{YCOORD}{y-coordinate, field, non-projected}
\item{ELEV}{elevation, m.}
\item{PROFCLASS}{soil class, obtained by classifying the soil profile at the
sample site }
\item{MAPCLASS}{soil class, obtained by looking up the site location in the soil map}
\item{VAL1}{Munsell colour component VALUE, 0-20 cm}
\item{CHR1}{Munsell colour component CHROMA, 20-40 cm}
\item{LIME1}{Lime content (tested using HCl), 0-20 cm}
\item{VAL2}{Munsell colour component VALUE, 0-20 cm}
\item{CHR2}{Munsell colour component CHROMA, 20-40 cm}
\item{LIME2}{Lime content (tested using HCl), 20-40 cm}
\item{DEPTHCM}{soil depth, cm }
\item{DEP2LIME}{depth to lime, cm }
\item{PCLAY1}{percentage clay, 0-20 cm}
\item{PCLAY2}{percentage clay, 20-40 cm}
\item{MG1}{Magnesium content (ppm), 0-20 cm}
\item{OM1}{organic matter (\%), 0-20 cm}
\item{CEC1}{CES as mequ/100g air dry soil, 0-20 cm}
\item{PH1}{pH, 0-20 cm}
\item{PHOS1}{Phosphorous, 0-20 cm, ppm}
\item{POT1}{K (potassium), 0-20 cm, ppm}
}
}
\usage{
data(oxford)
}
\author{ P.A. Burrough; compiled for R by Edzer Pebesma }
\references{ 
P.A. Burrough, R.A. McDonnell, 1998. Principles of Geographical Information
Systems. Oxford University Press.
}
\note{ \code{oxford.jpg}, in the gstat package external directory (see
example below), shows an image of the soil map for the region }
\keyword{datasets}
\examples{
data(oxford)
summary(oxford)
# open the following file with a jpg viewer:
system.file("external/oxford.jpg", package="gstat")
}
