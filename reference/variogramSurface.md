# Semivariance values for a given spatio-temporal variogram model

Generates a surface of semivariance values given a spatio-temporal
variogram model (one of separable, productSum, sumMetric,
simpleSumMetric or metric)

## Usage

``` r
variogramSurface(model, dist_grid, covariance = FALSE)
```

## Arguments

- model:

  A spatio-temporal variogram model generated through
  [`vgmST`](vgmST.md) or [`fit.StVariogram`](fit.StVariogram.md).

- dist_grid:

  A data.frame with two columns: `spacelag` and `timelag`.

- covariance:

  Whether the covariance should be computed instead of the variogram
  (default: FALSE).

## Value

A data.frame with columns `spacelag`, `timelag` and `gamma`.

## Author

Benedikt Graeler

## See also

See [`variogramLine`](variogramLine.md) for the spatial version and
[`fit.StVariogram`](fit.StVariogram.md) for the estimation of
spatio-temporal variograms.

## Examples

``` r
separableModel <- vgmST("separable", 
                        space=vgm(0.86, "Exp", 476, 0.14),
                        time =vgm(   1, "Exp",   3, 0),
                        sill=113)

data(vv)

if(require(lattice)) {
plot(vv, separableModel, wireframe=TRUE, all=TRUE)
}


# plotting of sample and model variogram
plot(vv, separableModel)

```
