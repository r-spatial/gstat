# Changelog

## version 2.1-2

CRAN release: 2024-09-05

- [`variogram()`](../reference/variogram.md) supports `stars` (raster)
  objects, benefiting from them being gridded

## version 2.1-1

CRAN release: 2023-04-06

## version 2.1-0

CRAN release: 2022-10-19

- import `sftime`; modify [`krigeST()`](../reference/krigeST.md)
  variogram functions to accept `sftime` objects for `data` (as
  alternative to `STI` or `STIDF`), and `stars` or `sftime` objects for
  `newdata`; [\#108](https://github.com/r-spatial/gstat/issues/108) with
  great help from [@henningte](https://github.com/henningte)

- fix Gaussian cosimulation, probably introduced in 2016,
  [\#106](https://github.com/r-spatial/gstat/issues/106)

## version 2.0-7

CRAN release: 2021-03-19

- return `NA` as estimate when prediction/simulation fails;
  [\#80](https://github.com/r-spatial/gstat/issues/80)

## version 2.0-6

CRAN release: 2020-05-18

- fixes `object 'ret' not found` bug introduced by
  [\#63](https://github.com/r-spatial/gstat/issues/63);
  [\#65](https://github.com/r-spatial/gstat/issues/65)
  [\#66](https://github.com/r-spatial/gstat/issues/66)
  [\#70](https://github.com/r-spatial/gstat/issues/70)

## version 2.0-5

CRAN release: 2020-04-04

- use multiple cores in `variogramST`, using pkg future;
  [\#63](https://github.com/r-spatial/gstat/issues/63) by
  [@sigmafelix](https://github.com/sigmafelix)

- fix bug with conditional simulation using `stars` target grid and
  nsim=1, [\#58](https://github.com/r-spatial/gstat/issues/58)

## version 2.0-4

CRAN release: 2020-01-21

- fix CRAN warning issue

## version 2.0-3

CRAN release: 2019-09-26

- fix bug in support for `sf` objects;
  [\#46](https://github.com/r-spatial/gstat/issues/46)

- fix `krigeTg` for the case when data or newdata are of class `sf` or
  `sfc`; [\#51](https://github.com/r-spatial/gstat/issues/51)

## version 2.0-1

CRAN release: 2019-04-25

- try to fix CRS issue found on OSX

- check for temporal unit in krigeST opposed to covariance function and
  assign the temporal unit found in krigeST to the spatio-temporal
  variogram for consistency in case none has been provided by the model;
  add warning when ST variogram doesn’t have a `temporal unit`
  attribute; [\#42](https://github.com/r-spatial/gstat/issues/42)

## version 2.0-0

CRAN release: 2019-02-28

- add plot.variogramModel;
  [\#40](https://github.com/r-spatial/gstat/issues/40)

- support `sf` and `stars` objects, for vector and raster data,
  respectively; [\#39](https://github.com/r-spatial/gstat/issues/39)

## version 1.1-6

CRAN release: 2018-04-02

- address warnings from Tomas Kalibera’s static code checking

## version 1.1-5

CRAN release: 2017-03-12

- fix auto-choosing of variogram parameters if only variogram model type
  is given, <https://github.com/edzer/gstat/issues/17>

## version 1.1-4

CRAN release: 2017-01-11

- Fix [\#17](https://github.com/r-spatial/gstat/issues/17)
- Fix [\#14](https://github.com/r-spatial/gstat/issues/14)
- Fix [\#12](https://github.com/r-spatial/gstat/issues/12)
- Fix great circle distance bug; see
  <https://stat.ethz.ch/pipermail/r-sig-geo/2016-December/025251.html>

## version 1.1-3

CRAN release: 2016-03-31

- Merge pull request [\#4](https://github.com/r-spatial/gstat/issues/4)
  from BenGraeler/master r-journal ms merge updates vignette
  “spatio-temporal-kriging” to revised paper
- Merge pull request [\#3](https://github.com/r-spatial/gstat/issues/3)
  from BenGraeler/master update st paper demos
- demo/stkrige.R: - remove empty station
- demo/00Index, demo/stkrige-crossvalidation.R,
  demo/stkrige-prediction.R, demo/stkrige.R, man/krigeST.Rd: - adds R
  scripts from vignette spatio-temporal kriging as demos
- Merge pull request [\#2](https://github.com/r-spatial/gstat/issues/2)
  from BenGraeler/master additions for space-time paper
- NAMESPACE, R/krige0.R, man/krigeST.Rd: - adds trans Gaussian kriging
  for space and time

## version 1.1-2

CRAN release: 2016-02-23

- fixed `memcpy` for overlapping regions error, found in valgrind checks
  by Brian Ripley
- fixed a few more (small) memory leaks

## version 1.1-1

CRAN release: 2016-02-21

- Further cleaned up C source code, got rid of lex and yacc files
- Improve `fit.variogram` to choose initial values for range, sill and
  nugget before fitting, and fit over a range of model types
- allow `NA` values in `vgm`
- Fixed [\#1](https://github.com/r-spatial/gstat/issues/1)

## version 1.1-0

CRAN release: 2015-10-18

- remove meschach matrix library, rewrote interface, link to R’s lapack
- improve notification in case of singular matrices
- remove all C code that was not used by R package gstat
- add `Makevars`, remove `configure`
- remove links to `ai-geostats.org`
- wrap `fit.StVariogram` example in `dontrun`

## version 1.0-26

CRAN release: 2015-08-26

- use ordered spatial index when selecting nearest strongest correlated
  neighbours in local kriging to avoid warning of spacetime
- update spatio-temporal geostatitics vignettes; add R Journal draft as
  vignette (Spatio-Temporal Geostatistics using gstat)
- add spatio-temporal PM10 interpolation
  [movie](http://gstat.r-forge.r-project.org/STpred.md)
