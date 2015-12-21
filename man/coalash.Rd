% $Id: coalash.Rd,v 1.2 2007-03-21 15:14:25 edzer Exp $
\name{coalash}
\alias{coalash}
\title{Coal ash samples from a mine in Pennsylvania}
\description{
Data obtained from Gomez and Hazen (1970, Tables 19 and 20) on coal ash
for the Robena Mine Property in Greene County Pennsylvania.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{x}{a numeric vector; x-coordinate; reference unknown }
    \item{y}{a numeric vector; x-coordinate; reference unknown }
	\item{coalash}{the target variable}
  } 
} 
\usage{
data(coalash)
}
\author{ unknown; R version prepared by Edzer Pebesma; data obtained from
\url{http://www.stat.uiowa.edu/~dzimmer/spatialstats/}, Dale Zimmerman's
course page }
\references{ N.A.C. Cressie, 1993, Statistics for Spatial Data,
Wiley. 

Gomez, M. and Hazen, K. (1970). Evaluating sulfur and ash
distribution in coal seems by statistical response surface regression
analysis. U.S. Bureau of Mines Report RI 7377.

see also fields manual: 
\url{http://www.image.ucar.edu/GSP/Software/Fields/fields.manual.coalashEX.Krig.shtml}
}
\note{ data are also present in package fields, as coalash.  }
\keyword{datasets}
\examples{
data(coalash)
summary(coalash)
}
