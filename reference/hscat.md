# Produce h-scatterplot

Produces h-scatterplots, where point pairs having specific separation
distances are plotted. This function is a wrapper around xyplot.

## Usage

``` r
hscat(formula, data, breaks, pch = 3, cex = .6, mirror = FALSE, 
  variogram.alpha = 0, as.table = TRUE,...)
```

## Arguments

- formula:

  specifies the dependent variable

- data:

  data where the variable in formula is resolved

- breaks:

  distance class boundaries

- pch:

  plotting symbol

- cex:

  plotting symbol size

- mirror:

  logical; duplicate all points mirrored along x=y? (note that
  correlations are those of the points plotted)

- variogram.alpha:

  parameter to be passed as alpha parameter to
  [variogram](variogram.md); if alpha is specified it will only affect
  xyplot by being passed through ...

- as.table:

  logical; if `TRUE`, panels plot top-to-bottom

- ...:

  parameters, passed to variogram and xyplot

## Value

an object of class trellis; normally the h scatter plot

## Author

Edzer Pebesma

## Note

Data pairs are plotted once, so the h-scatterplot are not symmetric.

## References

<http://www.gstat.org/>

Pebesma, E.J., 2004. Multivariable geostatistics in S: the gstat
package. Computers and Geosciences, 30: 683-691.

## Examples

``` r
library(sp)
data(meuse)
coordinates(meuse) = ~x+y
hscat(log(zinc)~1, meuse, c(0, 80, 120, 250, 500, 1000))
```
