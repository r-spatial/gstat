# GLS fitting of variogram parameters

Fits variogram parameters (nugget, sill, range) to variogram cloud,
using GLS (generalized least squares) fitting. Only for direct
variograms.

## Usage

``` r
fit.variogram.gls(formula, data, model, maxiter = 30, 
    eps = .01, trace = TRUE, ignoreInitial = TRUE, cutoff = Inf,
    plot = FALSE)
```

## Arguments

- formula:

  formula defining the response vector and (possible) regressors; in
  case of absence of regressors, use e.g. `z~1`

- data:

  object of class Spatial

- model:

  variogram model to be fitted, output of `vgm`

- maxiter:

  maximum number of iterations

- eps:

  convergence criterium

- trace:

  logical; if TRUE, prints parameter trace

- ignoreInitial:

  logical; if FALSE, initial parameter are taken from model; if TRUE,
  initial values of model are ignored and taken from variogram cloud:
  nugget: `mean(y)/2`, sill: `mean(y)/2`, range `median(h0)/4` with `y`
  the semivariance cloud value and `h0` the distances

- cutoff:

  maximum distance up to which point pairs are taken into consideration

- plot:

  logical; if TRUE, a plot is returned with variogram cloud and fitted
  model; else, the fitted model is returned.

## Value

an object of class "variogramModel"; see
[fit.variogram](fit.variogram.md); if `plot` is TRUE, a plot is returned
instead.

## References

Mueller, W.G., 1999: Least-squares fitting from the variogram cloud.
Statistics and Probability Letters, 43, 93-98.

Mueller, W.G., 2007: Collecting Spatial Data. Springer, Heidelberg.

## Author

Edzer Pebesma

## Note

Inspired by the code of Mihael Drinovac, which was again inspired by
code from Ernst Glatzer, author of package vardiag.

## See also

[fit.variogram](fit.variogram.md),

## Examples

``` r
library(sp)
data(meuse)
coordinates(meuse) = ~x+y
if (FALSE) { # \dontrun{
fit.variogram.gls(log(zinc)~1, meuse[1:40,], vgm(1, "Sph", 900,1))
} # }
```
