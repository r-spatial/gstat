# Constructing a spatio-temporal variogram

Constructs a spatio-temporal variogram of a given type checking for a
minimal set of parameters.

## Usage

``` r
vgmST(stModel, ..., space, time, joint, sill, k, nugget, stAni, temporalUnit)
```

## Arguments

- stModel:

  A string identifying the spatio-temporal variogram model (see details
  below). Only the string before an optional "\_" is used to identify
  the model. This mechanism can be used to identify different fits of
  the same model (`separable_A` and `separable_B` will be interpreted as
  separable models, but carry different names).

- ...:

  unused, but ensure an exact match of the following parameters.

- space:

  A spatial variogram.

- time:

  A temporal variogram.

- joint:

  A joint spatio-temporal variogram.

- sill:

  A joint spatio-temporal sill.

- k:

  The weighting of the product in the product-sum model.

- nugget:

  A joint spatio-temporal nugget.

- stAni:

  A spatio-temporal anisotropy; the number of space units equivalent to
  one time unit.

- temporalUnit:

  length one character vector, indicating the temporal unit (like secs)

## Details

The different implemented spatio-temporal variogram models have the
following required parameters (see as well the example section)

- separable::

  A variogram for `space` and `time` each and a joint spatio-temporal
  `sill` (variograms may have a separate nugget effect, but their joint
  sill will be 1) generating the call

      vgmST("separable", space, time, sill)

- productSum::

  A variogram for `space` and `time` each, and the weighting of product
  `k` generating the call

      vgmST("productSum", space, time, k)

- sumMetric::

  A variogram (potentially including a nugget effect) for `space`,
  `time` and `joint` each and a spatio-temporal anisotropy ratio `stAni`
  generating the call

      vgmST("sumMetric", space, time, joint, stAni)

- simpleSumMetric::

  A variogram (without nugget effect) for `space`, `time` and `joint`
  each, a joint spatio-temporal `nugget` effect and a spatio-temporal
  anisotropy ratio `stAni` generating the call

      vgmST("simpleSumMetric", space, time, joint, nugget, stAni)

- metric::

  A spatio-temporal `joint` variogram (potentially including a nugget
  effect) and `stAni` generating the call

      vgmST("metric", joint, stAni)

## Value

Returns an S3 object of class `StVariogramModel`.

## Author

Benedikt Graeler

## See also

[`fit.StVariogram`](fit.StVariogram.md) for fitting,
[`variogramSurface`](variogramSurface.md) to plot the variogram and
[`extractParNames`](extractPar.md) to better understand the parameter
structure of spatio-temporal variogram models.

## Examples

``` r
# separable model: spatial and temporal sill will be ignored
# and kept constant at 1-nugget respectively. A joint sill is used.
separableModel <- vgmST("separable", 
                        space=vgm(0.9,"Exp", 147, 0.1),
                        time =vgm(0.9,"Exp", 3.5, 0.1),
                        sill=40)

# product sum model: spatial and temporal nugget will be ignored and kept
# constant at 0. Only a joint nugget is used.
prodSumModel <- vgmST("productSum",
                      space=vgm(39, "Sph", 343, 0),
                      time= vgm(36, "Exp",   3, 0), 
                      k=15)

# sum metric model: spatial, temporal and joint nugget will be estimated
sumMetricModel <- vgmST("sumMetric",
                        space=vgm( 6.9, "Lin", 200, 3.0),
                        time =vgm(10.3, "Lin",  15, 3.6),
                        joint=vgm(37.2, "Exp",  84,11.7),
                        stAni=77.7)
                       
# simplified sumMetric model, only a overall nugget is fitted. The spatial, 
# temporal and jont nuggets are set to 0.
simpleSumMetricModel <- vgmST("simpleSumMetric",
                              space=vgm(20,"Lin", 150, 0),
                              time =vgm(20,"Lin", 10,  0),
                              joint=vgm(20,"Exp", 150, 0),
                              nugget=1, stAni=15)

# metric model
metricModel <- vgmST("metric",
                     joint=vgm(60, "Exp", 150, 10),
                     stAni=60)
```
