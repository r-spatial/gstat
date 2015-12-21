% $Id: sic97.Rd,v 1.2 2008-10-30 13:47:05 edzer Exp $
\name{sic97}
\alias{sic_obs}
\alias{sic_full}
\alias{demstd}
\title{ Spatial Interpolation Comparison 1997 data set: Swiss Rainfall}
\description{
The text below is copied from the data item at ai-geostats,
\url{https://wiki.52north.org/bin/view/AI_GEOSTATS/WebHome}

}

\format{
  The data frames contain the following columns:
  \describe{
   \item{ID}{this integer value is the number (unique value) of
   the monitoring station}
   \item{rainfall}{ rainfall amount, in 10th of mm }
  }
}
\note{
See the pdf that accompanies the original file for a description of the data.
The .dxf file with the Swiss border is not included here.
}

\usage{
data(sic97) # 
}
\author{ 
Gregoire Dubois and others.
}

\references{ \url{https://wiki.52north.org/bin/view/AI_GEOSTATS/WebHome}}

\keyword{datasets}
\examples{
data(sic97) 
image(demstd)
points(sic_full, pch=1)
points(sic_obs, pch=3)
}
