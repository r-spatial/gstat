% $Id: pcb.Rd,v 1.8 2009-07-03 12:10:55 edzer Exp $
\name{pcb}
\alias{pcb}
\title{PCB138 measurements in sediment at the NCP, the Dutch part of the North Sea}
\description{
PCB138 measurements in sediment at the NCP, which is the 
Dutch part of the North Sea
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{year}{ measurement year }
    \item{x}{ x-coordinate; UTM zone 31 }
    \item{y}{ y-coordinate; UTM zone 31 }
    \item{coast}{ distance to coast of the Netherlands, in km.}
    \item{depth}{ sea water depth, m. }
    \item{PCB138}{ PCB-138, measured on the sediment fraction smaller than
	63 \eqn{\mu}{mu}, in \eqn{\mu g/kg}{mu g/kg} dry matter; BUT SEE NOTE BELOW }
	\item{yf}{ year; as factor}
  }
}

\usage{
data(pcb)
}
\note{ A note of caution: The PCB-138 data are provided only to be able
to re-run the analysis done in Pebesma and Duin (2004; see references
below). If you want to use these data for comparison with PCB measurements
elsewhere, or if you want to compare them to regulation standards, or
want to use these data for any other purpose, you should first contact
\url{mailto:basisinfodesk@rikz.rws.minvenw.nl}.  The reason for this is
that several normalisations were carried out that are not reported here,
nor in the paper below. }

\references{ 
Pebesma, E. J., and Duin, R. N. M. (2005). Spatial patterns of temporal
change in North Sea sediment quality on different spatial scales. In
P. Renard, H. Demougeot-Renard and R. Froidevaux (Eds.), Geostatistics for
Environmental Applications: Proceedings of the Fifth European Conference
on Geostatistics for Environmental Applications (pp. 367-378): Springer.
}

\seealso{\link{ncp.grid}}
\keyword{datasets}
\examples{
data(pcb)
library(lattice)
xyplot(y~x|as.factor(yf), pcb, aspect = "iso")
# demo(pcb)
}
