\name{NEWS} 
\title{News for Package 'gstat'}
\section{Changes in version 1.1-5 (2017-03-10)}{
 \itemize{
  \item fix auto-choosing of variogram parameters if only variogram model type is given,
        https://github.com/edzer/gstat/issues/17
 }
}
\section{Changes in version 1.1-4 (2017-01-10)}{
 \itemize{
   \item Fix https://github.com/edzer/gstat/issues/17
   \item Fix https://github.com/edzer/gstat/issues/14
   \item Fix https://github.com/edzer/gstat/issues/12
   \item Fix great circle distance bug; see https://stat.ethz.ch/pipermail/r-sig-geo/2016-December/025251.html
 }
}
\section{Changes in version 1.1-3 (2016-03-30)}{
  \itemize{
	\item Merge pull request #4 from BenGraeler/master r-journal ms merge
	 updates vignette "spatio-temporal-kriging" to revised paper
	\item Merge pull request #3 from BenGraeler/master update st paper demos
	\item demo/stkrige.R: - remove empty station
	\item demo/00Index, demo/stkrige-crossvalidation.R,
	demo/stkrige-prediction.R, demo/stkrige.R, man/krigeST.Rd: - adds R scripts from vignette spatio-temporal kriging as demos
	\item Merge pull request #2 from BenGraeler/master additions for space-time paper
	\item NAMESPACE, R/krige0.R, man/krigeST.Rd: - adds trans Gaussian kriging for space and time
  }
}

\section{Changes in version 1.1-2 (2016-02-22)}{
  \itemize{
    \item fixed \code{memcpy} for overlapping regions error, found in valgrind checks by Brian Ripley
	\item fixed a few more (small) memory leaks
  }
}

\section{Changes in version 1.1-1 (2016-02-19)}{
  \itemize{
	\item Further cleaned up C source code, got rid of lex and yacc files
	\item Improve \code{fit.variogram} to choose initial values for range, sill and nugget before fitting, and fit over a range of model types
	\item allow \code{NA} values in \code{vgm}
	\item Fixed \url{https://github.com/edzer/gstat/issues/1}
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
  \item update spatio-temporal geostatitics vignettes; add R Journal draft as vignette
  (Spatio-Temporal Geostatistics using gstat)
  \item add spatio-temporal PM10 interpolation 
   \href{http://gstat.r-forge.r-project.org/STpred.html}{movie (mp4)}
  }
}
