# Meuse river altitude data set

This data set gives a point set with altitudes, digitized from the
1:10,000 topographical map of the Netherlands.

## Usage

``` r
data(meuse.alt)
```

## Format

This data frame contains the following columns:

- x:

  a numeric vector; x-coordinate (m) in RDM (Dutch topographical map
  coordinates)

- y:

  a numeric vector; y-coordinate (m) in RDM (Dutch topographical map
  coordinates)

- alt:

  altitude in m. above NAP (Dutch zero for sea level)

## References

<http://www.gstat.org/>

## See also

[meuse.all](meuse.all.md)

## Examples

``` r
data(meuse.alt)
library(lattice)
xyplot(y~x, meuse.alt, aspect = "iso")
```
