# version 2.0-6

* fixes `object 'ret' not found` bug introduced by #63; #65 #66 #70

# version 2.0-5

* use multiple cores in `variogramST`, using pkg future; #63 by @sigmafelix

* fix bug with conditional simulation using `stars` target grid and nsim=1, #58

# version 2.0-4

* fix CRAN warning issue

# version 2.0-3

* fix bug in support for `sf` objects; #46

* fix `krigeTg` for the case when data or newdata are of class `sf` or `sfc`; #51

# version 2.0-1

* try to fix CRS issue found on OSX

* check for temporal unit in krigeST opposed to covariance function and assign the temporal unit found in krigeST to the spatio-temporal variogram for consistency in case none has been provided by the model; add warning when ST variogram doesn't have a `temporal unit` attribute; #42

# version 2.0-0

* add plot.variogramModel; #40

* support `sf` and `stars` objects, for vector and raster data, respectively; #39

#  version 1.1-6 
 
  * address warnings from Tomas Kalibera's static code checking

#  version 1.1-5
 
  * fix auto-choosing of variogram parameters if only variogram model type is given,
        https://github.com/edzer/gstat/issues/17

#  version 1.1-4
 
   * Fix #17
   * Fix #14
   * Fix #12
   * Fix great circle distance bug; see https://stat.ethz.ch/pipermail/r-sig-geo/2016-December/025251.html

#  version 1.1-3
  
* Merge pull request #4 from BenGraeler/master r-journal ms merge updates vignette "spatio-temporal-kriging" to revised paper
* Merge pull request #3 from BenGraeler/master update st paper demos
* demo/stkrige.R: - remove empty station
* demo/00Index, demo/stkrige-crossvalidation.R, demo/stkrige-prediction.R, demo/stkrige.R, man/krigeST.Rd: - adds R scripts from vignette spatio-temporal kriging as demos
* Merge pull request #2 from BenGraeler/master additions for space-time paper
* NAMESPACE, R/krige0.R, man/krigeST.Rd: - adds trans Gaussian kriging for space and time

#  version 1.1-2
  
* fixed `memcpy` for overlapping regions error, found in valgrind checks by Brian Ripley
* fixed a few more (small) memory leaks

#  version 1.1-1
  
* Further cleaned up C source code, got rid of lex and yacc files
* Improve `fit.variogram` to choose initial values for range, sill and nugget before fitting, and fit over a range of model types
* allow `NA` values in `vgm`
* Fixed #1

#  version 1.1-0
  
* remove meschach matrix library, rewrote interface, link to R's lapack
* improve notification in case of singular matrices
* remove all C code that was not used by R package gstat
* add `Makevars`, remove `configure`
* remove links to `ai-geostats.org`
* wrap `fit.StVariogram` example in `dontrun`

#  version 1.0-26
  
* use ordered spatial index when selecting nearest strongest correlated neighbours in local kriging to avoid warning of spacetime
* update spatio-temporal geostatitics vignettes; add R Journal draft as vignette (Spatio-Temporal Geostatistics using gstat)
* add spatio-temporal PM10 interpolation [movie](http://gstat.r-forge.r-project.org/STpred.html)

#  version before 1.0-26

* see [ChangeLog](https://github.com/r-spatial/gstat/blob/master/inst/ChangeLog)
