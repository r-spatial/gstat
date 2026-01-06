# conditional/unconditional spatio-temporal simulation

conditional/unconditional spatio-temporal simulation based on turning
bands

## Usage

``` r
krigeSTSimTB(formula, data, newdata, modelList, nsim, progress = TRUE, 
             nLyrs = 500, tGrid = NULL, sGrid = NULL, ceExt = 2, nmax = Inf)
```

## Arguments

- formula:

  the formula of the kriging predictor

- data:

  conditioning data

- newdata:

  locations in space and time where the simulation is carried out

- modelList:

  the spatio-temporal variogram (from [`vgmST`](vgmST.md)) defining the
  spatio-temporal covariance structure of the simulated Gaussian random
  field

- nsim:

  number of simulations

- progress:

  boolean; whether the progress should be shown in progress bar

- nLyrs:

  number of layers used in the turning bands approach (default = 500)

- tGrid:

  optional explicit temporal griding that shall be used

- sGrid:

  optional explicit spatial griding that shall be used

- ceExt:

  expansion in the circulant embedding, defaults to 2

- nmax:

  number of nearest neighbours that shall e used, defaults to 'Inf'
  meaning all available points are used

## Value

a spatio-temporal data frame with `nSim` simulations

## References

Turning bands

Lantuejoul, C. (2002) Geostatistical Simulation: Models and Algorithms.
Springer.

Matheron, G. (1973). The intrinsic random functions and their
applications. Adv. Appl. Probab., 5, 439-468.

Strokorb, K., Ballani, F., and Schlather, M. (2014) Tail correlation
functions of max-stable processes: Construction principles, recovery and
diversity of some mixing max-stable processes with identical TCF.
Extremes, Submitted.

Turning layers

Schlather, M. (2011) Construction of covariance functions and
unconditional simulation of random fields. In Porcu, E., Montero, J.M.
and Schlather, M., Space-Time Processes and Challenges Related to
Environmental Problems. New York: Springer.

## Author

Benedikt Graeler

## See also

[`krigeSimCE`](krigeSimCE.md)

## Examples

``` r
# see demo('circEmbeddingMeuse')
```
