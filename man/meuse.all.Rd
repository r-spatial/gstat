% $Id: meuse.all.Rd,v 1.5 2006-02-10 19:03:27 edzer Exp $
\name{meuse.all}
\alias{meuse.all}
\title{Meuse river data set -- original, full data set}
\description{
This data set gives locations and top soil heavy metal concentrations
(ppm), along with a number of soil and landscape variables, collected in
a flood plain of the river Meuse, near the village Stein. Heavy metal
concentrations are bulk sampled from an area of approximately 15 m
x 15 m.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{sample}{sample number}
    \item{x}{a numeric vector; x-coordinate (m)
	in RDM (Dutch topographical map coordinates) }
    \item{y}{a numeric vector; y-coordinate (m)
	in RDM (Dutch topographical map coordinates)}
	\item{cadmium}{topsoil cadmium concentration, ppm.; note that zero cadmium
	values in the original data set have been shifted to 0.2 (half the lowest 
	non-zero value) }
	\item{copper}{topsoil copper concentration, ppm. }
	\item{lead}{topsoil lead concentration, ppm. }
	\item{zinc}{topsoil zinc concentration, ppm. } 
	\item{elev}{relative elevation}
	\item{om}{organic matter, as percentage }
	\item{ffreq}{flooding frequency class}
	\item{soil}{soil type}
	\item{lime}{lime class}
	\item{landuse}{landuse class}
	\item{dist.m}{distance to river Meuse (metres), as obtained during
	the field survey}
	\item{in.pit}{logical; indicates whether this is a sample taken in
	a pit}
	\item{in.meuse155}{logical; indicates whether the sample is part of
	the \code{meuse} (i.e., filtered) data set; in addition to the samples
	in a pit, an sample (139) with outlying zinc content was removed }
	\item{in.BMcD}{logical; indicates whether the sample is used as part
	of the subset of 98 points in the various interpolation examples of
	Burrough and McDonnell}
  } 
} 
\usage{
data(meuse.all)
}
\author{ The actual field data were collected by Ruud van Rijn and
Mathieu Rikken; data compiled for R by Edzer Pebesma }
\references{ 
P.A. Burrough, R.A. McDonnell, 1998. Principles of Geographical Information
Systems. Oxford University Press.

\url{http://www.gstat.org/}
}
\note{ \code{sample} refers to original sample number.  Eight samples
were left out because they were not indicative for the metal content of
the soil. They were taken in an old pit. One sample contains an outlying
zinc value, which was also discarded for the meuse (155) data set.
}
\seealso{\link{meuse.alt}}
\keyword{datasets}
\examples{
data(meuse.all)
summary(meuse.all)
}
