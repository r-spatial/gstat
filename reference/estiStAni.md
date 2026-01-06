# Estimation of the spatio-temporal anisotropy

Estimation of the spatio-temporal anisotropy without an underlying
spatio-temporal model. Different methods are implemented using a linear
model to predict the temporal gamma values or the ratio of the ranges of
a spatial and temporal variogram model or a spatial variogram model to
predict the temporal gamma values or the spatio-temporal anisotropy
value as used in a metric spatio-temporal variogram.

## Usage

``` r
estiStAni(empVgm, interval, method = "linear", spatialVgm, 
          temporalVgm, s.range=NA, t.range=NA)
```

## Arguments

- empVgm:

  An empirical spatio-temporal variogram.

- interval:

  A search interval for the optimisation of the spatio-temporal
  anisotropy parameter

- method:

  A character string determining the method to be used (one of `linear`,
  `range`, `vgm` or `metric`, see below for details)

- spatialVgm:

  A spatial variogram definition from the call to [`vgm`](vgm.md). The
  model is optimised based on the pure spatial values in `empVgm`.

- temporalVgm:

  A temporal variogram definition from the call to [`vgm`](vgm.md). The
  model is optimised based on the pure temporal values in `empVgm`.

- s.range:

  A spatial cutoff value applied to the empirical variogram `empVgm`.

- t.range:

  A temporal cutoff value applied to the empirical variogram `empVgm`.

## Details

- linear:

  A linear model is fitted to the pure spatial gamma values based on the
  spatial distances. An optimal scaling is searched to stretch the
  temporal distances such that the linear model explains best the pure
  temporal gamma values. This assumes (on average) a linear relationship
  between distance and gamma, hence it is advisable to use only those
  pairs of pure spatial (pure temporal) distance and gamma value that
  show a considerable increase (i.e. drop all values beyond the range by
  setting values for `s.range` and `t.range`).

- range:

  A spatial and temporal variogram model is fitted to the pure spatial
  and temporal gamma values respectively. The spatio-temporal anisotropy
  estimate is the ratio of the spatial range over the temporal range.

- vgm:

  A spatial variogram model is fitted to the pure spatial gamma values.
  An optimal scaling is used to stretch the temporal distances such that
  the spatial variogram model explains best the pure temporal gamma
  values.

- metric:

  A metric spatio-temporal variogram model is fitted with `joint`
  component according to the defined spatial variogram `spatialVgm`. The
  starting value of `stAni` is the mean of the `interval` parameter (see
  [`vgmST`](vgmST.md) for the metric variogram definition). The
  spatio-temporal anisotropy as estimated in the spatio-temporal
  variogram is returned. Note that the parameter `interval` is only used
  to set the starting value. Hence, the estimate might exceed the given
  interval.

## Value

A scalar representing the spatio-temporal anisotropy estimate.

## Note

Different methods might lead to very different estimates. All but the
`linear` approach are sensitive to the variogram model selection.

## Author

Benedikt Graeler

## Examples

``` r
data(vv)

estiStAni(vv, c(10, 150))
#> [1] 83.57463
estiStAni(vv, c(10, 150), "vgm", vgm(80, "Sph", 120, 20))
#> [1] 62.83509
```
