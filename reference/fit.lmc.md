# Fit a Linear Model of Coregionalization to a Multivariable Sample Variogram

Fit a Linear Model of Coregionalization to a Multivariable Sample
Variogram; in case of a single variogram model (i.e., no nugget) this is
equivalent to Intrinsic Correlation

## Usage

``` r
fit.lmc(v, g, model, fit.ranges = FALSE, fit.lmc = !fit.ranges, 
correct.diagonal = 1.0, ...)
```

## Arguments

- v:

  multivariable sample variogram, output of [variogram](variogram.md)

- g:

  gstat object, output of [gstat](gstat.md)

- model:

  variogram model, output of [vgm](vgm.md); if supplied this value is
  used as initial value for each fit

- fit.ranges:

  logical; determines whether the range coefficients (excluding that of
  the nugget component) should be fitted; or logical vector: determines
  for each range parameter of the variogram model whether it should be
  fitted or fixed.

- fit.lmc:

  logical; if TRUE, each coefficient matrices of partial sills is
  guaranteed to be positive definite

- correct.diagonal:

  multiplicative correction factor to be applied to partial sills of
  direct variograms only; the default value, 1.0, does not correct. If
  you encounter problems with singular covariance matrices during
  cokriging or cosimulation, you may want to try to increase this to
  e.g. 1.01

- ...:

  parameters that get passed to [fit.variogram](fit.variogram.md)

## Value

returns an object of class `gstat`, with fitted variograms;

## Author

Edzer Pebesma

## Note

This function does not use the iterative procedure proposed by M.
Goulard and M. Voltz (Math. Geol., 24(3): 269-286; reproduced in
Goovaerts' 1997 book) but uses simply two steps: first, each variogram
model is fitted to a direct or cross variogram; next each of the partial
sill coefficient matrices is approached by its in least squares sense
closest positive definite matrices (by setting any negative eigenvalues
to zero).

The argument `correct.diagonal` was introduced by experience: by zeroing
the negative eigenvalues for fitting positive definite partial sill
matrices, apparently still perfect correlation may result, leading to
singular cokriging/cosimulation matrices. If someone knows of a more
elegant way to get around this, please let me know.

## See also

[variogram](variogram.md), [vgm](vgm.md),
[fit.variogram](fit.variogram.md), `demo(cokriging)`
