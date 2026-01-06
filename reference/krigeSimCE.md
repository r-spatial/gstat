# Simulation based on circulant embedding

Simulating a conditional/unconditional Gaussian random field via kriging
and circulant embedding

## Usage

``` r
krigeSimCE(formula, data, newdata, model, n = 1, ext = 2)
```

## Arguments

- formula:

  the formula of the kriging predictor

- data:

  spatial data frame that conditions the simulation

- newdata:

  locations in space where the Gaussian random field shall be simulated

- model:

  a vgm model that defines the spatial covariance structure

- n:

  number of simulations

- ext:

  extension factor of the circulant embedding, defaults to 2

## Value

A spatial data frame as defined in `newdata` with `n` simulations.

## References

Davies, Tilman M., and David Bryant: "On circulant embedding for
Gaussian random fields in R." Journal of Statistical Software 55.9
(2013): 1-21. See i.e. the supplementary files at (retrieved
2018-05-25):
https://www.jstatsoft.org/index.php/jss/article/downloadSuppFile/v055i09/v55i09.R

## Author

Benedikt Graeler

## See also

[`krigeSTSimTB`](krigeSTSimTB.md)

## Examples

``` r
# see demo('circEmbeddingMeuse')
```
