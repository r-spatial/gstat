# Function that returns the covariances for areas

Function that returns the covariances for areas based on spatio-temporal
point variograms for use in the spatio-temporal area-to-point kriging

## Usage

``` r
vgmAreaST(x, y = x, model, ndiscrSpace = 16, verbose = FALSE, covariance = TRUE)
```

## Arguments

- x:

  spatio-temporal data frame

- y:

  spatio-temporal data frame

- model:

  spatio-temporal variogram model for point support

- ndiscrSpace:

  number of discretisation in space

- verbose:

  Boolean: default to FALSE, set to TRUE for debugging

- covariance:

  Boolean: whether the covariance shall be evaluated, currently
  disfunction and set to TRUE

## Value

The covariance between 'x' and 'y'.

## Author

Benedikt Graeler

## See also

[`vgmArea`](vgmArea.md)

## Examples

``` r
# see demo('a2pinST')
```
