\name{vgmAreaST}
\alias{vgmAreaST}
\title{
Function that returns the covariances for areas
}
\description{
Function that returns the covariances for areas based on spatio-temporal point variograms for use in the spatio-temporal area-to-point kriging
}
\usage{
vgmAreaST(x, y = x, model, ndiscrSpace = 16, verbose = FALSE, covariance = TRUE)
}
\arguments{
  \item{x}{spatio-temporal data frame}
  \item{y}{spatio-temporal data frame}
  \item{model}{spatio-temporal variogram model for point support}
  \item{ndiscrSpace}{number of discretisation in space}
  \item{verbose}{Boolean: default to FALSE, set to TRUE for debugging}
  \item{covariance}{Boolean: whether the covariance shall be evaluated, currently disfunction and set to TRUE}
}
\value{The covariance between 'x'  and 'y'.}
\author{Benedikt Graeler}

\seealso{\code{\link{vgmArea}}}
\examples{
# see demo('a2pinST')
}
