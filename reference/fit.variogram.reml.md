# REML Fit Direct Variogram Partial Sills to Data

Fit Variogram Sills to Data, using REML (only for direct variograms; not
for cross variograms)

## Usage

``` r
fit.variogram.reml(formula, locations, data, model, debug.level = 1, set, degree = 0)
```

## Arguments

- formula:

  formula defining the response vector and (possible) regressors; in
  case of absence of regressors, use e.g. `z~1`

- locations:

  spatial data locations; a formula with the coordinate variables in the
  right hand (dependent variable) side.

- data:

  data frame where the names in formula and locations are to be found

- model:

  variogram model to be fitted, output of `vgm`

- debug.level:

  debug level; set to 65 to see the iteration trace and log likelihood

- set:

  additional options that can be set; use `set=list(iter=100)` to set
  the max. number of iterations to 100.

- degree:

  order of trend surface in the location, between 0 and 3

## Value

an object of class "variogramModel"; see
[fit.variogram](fit.variogram.md)

## References

Christensen, R. Linear models for multivariate, Time Series, and Spatial
Data, Springer, NY, 1991.

Kitanidis, P., Minimum-Variance Quadratic Estimation of Covariances of
Regionalized Variables, Mathematical Geology 17 (2), 195–208, 1985

## Author

Edzer Pebesma

## Note

This implementation only uses REML fitting of sill parameters. For each
iteration, an \\n \times n\\ matrix is inverted, with \$n\$ the number
of observations, so for large data sets this method becomes demanding. I
guess there is much more to likelihood variogram fitting in package
`geoR`, and probably also in `nlme`.

## See also

[fit.variogram](fit.variogram.md),

## Examples

``` r
library(sp)
data(meuse)
fit.variogram.reml(log(zinc)~1, ~x+y, meuse, model = vgm(1, "Sph", 900,1))
#>   model      psill range
#> 1   Nug 0.02549524     0
#> 2   Sph 0.61080053   900
```
