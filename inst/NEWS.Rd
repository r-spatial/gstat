\name{NEWS} 
\title{News for Package 'gstat'}
\section{Changes in version 1.1-1 (2015-11-xx)}{
  \itemize{
	\item Further cleaned up C source code, got rid of lex and yacc files
  }
}

\section{Changes in version 1.1-0 (2015-10-17)}{
  \itemize{
	\item remove meschach matrix library, rewrote interface, link to R's lapack
	\item improve notification in case of singular matrices
	\item remove all C code that was not used by R package gstat
    \item add \code{Makevars}, remove \code{configure}
	\item remove links to \code{ai-geostats.org}
	\item wrap \code{fit.StVariogram} example in \code{dontrun}
  }
}

\section{Changes in version 1.0-26 (2015-08-26)}{
  \itemize{
  \item use ordered spatial index when selecting nearest strongest correlated
     neighbours in local kriging to avoid warning of spacetime
  \item update spatio-temporal geostatitics vignettes; add R Journal 
   \href{http://cran.uni-muenster.de/web/packages/gstat/vignettes/spatio-temporal-kriging.pdf}{draft}
  \item add spatio-temporal PM10 interpolation 
   \href{http://gstat.r-forge.r-project.org/STpred.html}{movie (mp4)}

  }
}
