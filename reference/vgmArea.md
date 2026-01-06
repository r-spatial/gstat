# point-point, point-area or area-area semivariance

Compute point-point, point-area or area-area variogram values from point
model

## Usage

``` r
vgmArea(x, y = x, vgm, ndiscr = 16, verbose = FALSE, covariance = TRUE)
```

## Arguments

- x:

  object of class
  [SpatialPoints](https://edzer.github.io/sp/reference/SpatialPoints.html)
  or
  [SpatialPolygons](https://edzer.github.io/sp/reference/SpatialPolygons.html)

- y:

  object of class
  [SpatialPoints](https://edzer.github.io/sp/reference/SpatialPoints.html)
  or
  [SpatialPolygons](https://edzer.github.io/sp/reference/SpatialPolygons.html)

- vgm:

  variogram model, see [vgm](vgm.md)

- ndiscr:

  number of points to discretize an area, using
  [spsample](https://edzer.github.io/sp/reference/spsample.html)

- verbose:

  give progress bar

- covariance:

  logical; compute covariances, rather than semivariances?

## Value

semivariance or covariance matrix of dimension `length(x)` x `lenght(y)`

## Author

Edzer Pebesma

## Examples

``` r
library(sp)
demo(meuse, ask = FALSE, echo = FALSE)
vgmArea(meuse[1:5,], vgm = vgm(1, "Exp", 1000)) # point-point
#>           [,1]      [,2]      [,3]      [,4]      [,5]
#> [1,] 1.0000000 0.9316129 0.8879422 0.7716384 0.6932850
#> [2,] 0.9316129 1.0000000 0.8679977 0.7536317 0.6958367
#> [3,] 0.8879422 0.8679977 1.0000000 0.8666057 0.7780038
#> [4,] 0.7716384 0.7536317 0.8666057 1.0000000 0.8570468
#> [5,] 0.6932850 0.6958367 0.7780038 0.8570468 1.0000000
vgmArea(meuse[1:5,], meuse.area, vgm = vgm(1, "Exp", 1000)) # point-area
#>           [,1]
#> [1,] 0.1565761
#> [2,] 0.1646099
#> [3,] 0.1648786
#> [4,] 0.1672019
#> [5,] 0.1867704
```
