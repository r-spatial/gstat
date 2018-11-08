\name{krigeSimCE}
\alias{krigeSimCE}
\title{Simulation based on circulant embedding}
\description{Simulating a conditional/unconditional Gaussian random field via kriging and circulant embedding}
\usage{
krigeSimCE(formula, data, newdata, model, n = 1, ext = 2)
}

\arguments{
  \item{formula}{ the formula of the kriging predictor}
  \item{data}{ spatial data frame that conditions the simulation}
  \item{newdata}{locations in space where the Gaussian random field shall be simulated}
  \item{model}{a vgm model that defines the spatial covariance structure}
  \item{n}{number of simulations}
  \item{ext}{extension factor of the circulant embedding, defaults to 2}
}
\value{A spatial data frame as defined in \code{newdata} with \code{n} simulations.}

\references{
Davies, Tilman M., and David Bryant: "On circulant embedding for Gaussian random fields in R." 
Journal of Statistical Software 55.9 (2013): 1-21.
See i.e. the supplementary files at (retrieved 2018-05-25): 
https://www.jstatsoft.org/index.php/jss/article/downloadSuppFile/v055i09/v55i09.R
}

\author{Benedikt Graeler}

\seealso{\code{\link{krigeSTSimTB}}}

\examples{# see demo('circEmbeddingMeuse')}