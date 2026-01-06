# Semivariance Values For a Given Variogram Model

Generates a semivariance values given a variogram model

## Usage

``` r
variogramLine(object, maxdist, n = 200, min = 1.0e-6 * maxdist, 
  dir = c(1,0,0), covariance = FALSE, ..., dist_vector, debug.level = 0)
```

## Arguments

- object:

  variogram model for which we want semivariance function values

- maxdist:

  maximum distance for which we want semivariance values

- n:

  number of points

- min:

  minimum distance; a value slightly larger than zero is usually used to
  avoid the discontinuity at distance zero if a nugget component is
  present

- dir:

  direction vector: unit length vector pointing the direction in x
  (East-West), y (North-South) and z (Up-Down)

- covariance:

  logical; if TRUE return covariance values, otherwise return
  semivariance values

- ...:

  ignored

- dist_vector:

  numeric vector or matrix with distance values

- debug.level:

  gstat internal debug level

## Value

a data frame of dimension (`n` x 2), with columns distance and gamma
(semivariances or covariances), or in case `dist_vector` is a matrix, a
conforming matrix with semivariance/covariance values is returned.

## Note

variogramLine is used to generate data for plotting a variogram model.

## Author

Edzer Pebesma

## See also

[plot.gstatVariogram](plot.gstatVariogram.md)

## Examples

``` r
variogramLine(vgm(5, "Exp", 10, 5), 10, 10)
#>        dist    gamma
#> 1   0.00001 5.000005
#> 2   1.11112 5.525807
#> 3   2.22223 5.996316
#> 4   3.33334 6.417346
#> 5   4.44445 6.794100
#> 6   5.55556 7.131234
#> 7   6.66667 7.432915
#> 8   7.77778 7.702871
#> 9   8.88889 7.944439
#> 10 10.00000 8.160603
# anisotropic variogram, plotted in E-W direction:
variogramLine(vgm(1, "Sph", 10, anis=c(0,0.5)), 10, 10)
#>        dist     gamma
#> 1   0.00001 0.0000030
#> 2   1.11112 0.3278489
#> 3   2.22223 0.6227728
#> 4   3.33334 0.8518530
#> 5   4.44445 0.9821677
#> 6   5.55556 1.0000000
#> 7   6.66667 1.0000000
#> 8   7.77778 1.0000000
#> 9   8.88889 1.0000000
#> 10 10.00000 1.0000000
# anisotropic variogram, plotted in N-S direction:
variogramLine(vgm(1, "Sph", 10, anis=c(0,0.5)), 10, 10, dir=c(0,1,0))
#>        dist     gamma
#> 1   0.00001 0.0000015
#> 2   1.11112 0.1659821
#> 3   2.22223 0.3278475
#> 4   3.33334 0.4814824
#> 5   4.44445 0.6227716
#> 6   5.55556 0.7475999
#> 7   6.66667 0.8518521
#> 8   7.77778 0.9314130
#> 9   8.88889 0.9821674
#> 10 10.00000 1.0000000
variogramLine(vgm(1, "Sph", 10, anis=c(0,0.5)), dir=c(0,1,0), dist_vector = 0.5)
#>   dist     gamma
#> 1  0.5 0.0749375
variogramLine(vgm(1, "Sph", 10, anis=c(0,0.5)), dir=c(0,1,0), dist_vector = c(0, 0.5, 0.75))
#>   dist     gamma
#> 1 0.00 0.0000000
#> 2 0.50 0.0749375
#> 3 0.75 0.1122891
```
