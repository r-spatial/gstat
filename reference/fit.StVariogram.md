# Fit a spatio-temporal sample variogram to a sample variogram

Fits a spatio-temporal variogram of a given type to spatio-temporal
sample variogram.

## Usage

``` r
fit.StVariogram(object, model, ...,  method = "L-BFGS-B",
  lower, upper, fit.method = 6, stAni=NA, wles)
```

## Arguments

- object:

  The spatio-temporal sample variogram. Typically output from
  [`variogramST`](variogramST.md)

- model:

  The desired spatio-temporal model defined through [`vgmST`](vgmST.md).

- ...:

  further arguments passed to
  [`optim`](https://rdrr.io/r/stats/optim.html).
  [`extractParNames`](extractPar.md) provides the parameter structure of
  spatio-temporal variogram models that help to provide sensible `upper`
  and `lower` limits.

- lower:

  Lower limits used by optim. If missing, the smallest well defined
  values are used (mostly near 0).

- upper:

  Upper limits used by optim. If missing, the largest well defined
  values are used (mostly `Inf`).

- method:

  fit method, pass to [`optim`](https://rdrr.io/r/stats/optim.html)

- fit.method:

  an integer between 0 and 13 determine the fitting routine (i.e.
  weighting of the squared residuals in the LSE). Values 0 to 6
  correspond with the pure spatial version (see
  [`fit.variogram`](fit.variogram.md)). See the details section for the
  meaning of the other values (partly experimental).

- stAni:

  The spatio-temporal anisotropy that is used in the weighting. Might be
  missing if the desired spatio-temporal variogram model does contain a
  spatio-temporal anisotropy parameter (this might cause bad convergence
  behaviour). The default is `NA` and will be understood as identity (1
  temporal unit = 1 spatial unit). As this only in very few cases a
  valid assumption, a warning is issued.

- wles:

  Should be missing; only for backwards compatibility, `wles = TRUE`
  corresponds to `fit.method = 1` and `wles = FALSE` corresponds to
  `fit.method = 6`.

## Details

The following list summarizes the meaning of the `fit.method` argument
which is essential a weighting of the squared residuals in the
least-squares estimation. Please note, that weights based on the models
gamma value might fail to converge properly due to the dependence of
weights on the variogram estimate:

- `fit.method = 0`:

  no fitting, however the MSE between the provided variogram model and
  sample variogram surface is calculated.

- `fit.method = 1`:

  Number of pairs in the spatio-temporal bin: \\N_j\\

- `fit.method = 2`:

  Number of pairs in the spatio-temporal bin divided by the square of
  the current variogram model's value: \\N_j/\gamma(h_j, u_j)^2\\

- `fit.method = 3`:

  Same as `fit.method = 1` for compatibility with
  [`fit.variogram`](fit.variogram.md) but as well evaluated in R.

- `fit.method = 4`:

  Same as `fit.method = 2` for compatibility with
  [`fit.variogram`](fit.variogram.md) but as well evaluated in R.

- `fit.method = 5`:

  Reserved for REML for compatibility with
  [`fit.variogram`](fit.variogram.md), not yet implemented.

- `fit.method = 6`:

  No weights.

- `fit.method = 7`:

  Number of pairs in the spatio-temporal bin divided by the square of
  the bin's metric distance. If `stAni` is not specified, the model's
  parameter is used to calculate the metric distance across space and
  time: \\N_j/(h_j^2 + {\rm stAni}^2\cdot u_j^2)\\

- `fit.method = 8`:

  Number of pairs in the spatio-temporal bin divided by the square of
  the bin's spatial distance. \\N_j/h_j^2\\. Note that the 0 distances
  are replaced by the smallest non-zero distances to avoid division by
  zero.

- `fit.method = 9`:

  Number of pairs in the spatio-temporal bin divided by the square of
  the bin's temporal distance. \\N_j/u_j^2\\. Note that the 0 distances
  are replaced by the smallest non-zero distances to avoid division by
  zero.

- `fit.method = 10`:

  Reciprocal of the square of the current variogram model's value:
  \\1/\gamma(h_j,u_j)^2\\

- `fit.method = 11`:

  Reciprocal of the square of the bin's metric distance. If `stAni` is
  not specified, the model's parameter is used to calculate the metric
  distance across space and time: \\1/(h_j^2 + {\rm stAni}^2\cdot
  u_j^2)\\

- `fit.method = 12`:

  Reciprocal of the square of the bin's spatial distance. \\1/h_j^2\\.
  Note that the 0 distances are replaced by the smallest non-zero
  distances to avoid division by zero.

- `fit.method = 13`:

  Reciprocal of the square of the bin's temporal distance. \\1/u_j^2\\.
  Note that the 0 distances are replaced by the smallest non-zero
  distances to avoid division by zero.

See also Table 4.2 in the gstat manual for the original spatial version.

## Value

Returns a spatio-temporal variogram model, as S3 class StVariogramModel.
It carries the temporal and spatial unit as attributes `"temporal unit"`
and `"spatial unit"` in order to allow [`krigeST`](krigeST.md) to adjust
for different units. The units are obtained from the provided empirical
variogram. Further attributes are the optim output `"optim.output"` and
the always not weighted mean squared error `"MSE"`.

## Author

Benedikt Graeler

## See also

[`fit.variogram`](fit.variogram.md) for the pure spatial case.
[`extractParNames`](extractPar.md) helps to understand the parameter
structure of spatio-temporal variogram models.

## Examples

``` r
# separable model: spatial and temporal sill will be ignored
# and kept constant at 1-nugget respectively. A joint sill is used.
if (FALSE) { # \dontrun{
separableModel <- vgmST("separable", 
                        method = "Nelder-Mead", # no lower & upper needed
                        space=vgm(0.9,"Exp", 123, 0.1),
                        time =vgm(0.9,"Exp", 2.9, 0.1),
                        sill=100)

data(vv)
separableModel <- fit.StVariogram(vv, separableModel,
                                  method="L-BFGS-B",
                                  lower=c(10,0,0.01,0,1),
                                  upper=c(500,1,20,1,200))
plot(vv, separableModel)
} # } # dontrun
```
