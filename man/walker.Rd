% $Id: walker.Rd,v 1.2 2006-02-10 19:03:27 edzer Exp $
\name{walker}
\alias{walker}
\alias{walker.exh}
\title{Walker Lake sample and exhaustive data sets}
\description{
This is the Walker Lake data sets (sample and exhaustive
data set), used in Isaaks and Srivastava's Applied Geostatistics.
}
\format{
  This data frame contains the following columns:
  \describe{
  \item{Id}{Identification Number}
  \item{X}{Xlocation in meter}
  \item{Y}{Ylocation in meter}
  \item{V}{V variable, concentration in ppm}
  \item{U}{U variable, concentration in ppm}
  \item{T}{T variable, indicator variable}
  }
}
\usage{
data(walker)
}
\references{ 
Applied Geostatistics
by Edward H. Isaaks, R. Mohan Srivastava;
Oxford University Press.
}
\note{ This data sets was obtained from the data sets on
ai-geostats (link no longer functioning) }
\keyword{datasets}
\examples{
library(sp)
data(walker)
summary(walker)
summary(walker.exh)
}
