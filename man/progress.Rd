% $Id: hscat.Rd,v 1.3 2008-02-04 10:06:44 edzer Exp $
\name{progress}
\alias{get_gstat_progress}
\alias{set_gstat_progress}
\title{ Get or set progress indicator }
\description{ Get or set progress indicator }
\usage{
get_gstat_progress()
set_gstat_progress(value)
}
\arguments{
\item{value}{ logical }
}
\value{
return the logical value indicating whether progress bars should be given
}
\author{ Edzer Pebesma }

\examples{
set_gstat_progress(FALSE)
get_gstat_progress()
}
\keyword{models}
