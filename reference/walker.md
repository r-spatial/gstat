# Walker Lake sample and exhaustive data sets

This is the Walker Lake data sets (sample and exhaustive data set), used
in Isaaks and Srivastava's Applied Geostatistics.

## Usage

``` r
data(walker)
```

## Format

This data frame contains the following columns:

- Id:

  Identification Number

- X:

  Xlocation in meter

- Y:

  Ylocation in meter

- V:

  V variable, concentration in ppm

- U:

  U variable, concentration in ppm

- T:

  T variable, indicator variable

## References

Applied Geostatistics by Edward H. Isaaks, R. Mohan Srivastava; Oxford
University Press.

## Note

This data sets was obtained from the data sets on ai-geostats (link no
longer functioning)

## Examples

``` r
library(sp)
data(walker)
summary(walker)
#> Object of class SpatialPointsDataFrame
#> Coordinates:
#>   min max
#> X   8 251
#> Y   8 291
#> Is projected: NA 
#> proj4string : [NA]
#> Number of points: 470
#> Data attributes:
#>        Id              V                U                 T        
#>  Min.   :  1.0   Min.   :   0.0   Min.   :   0.00   Min.   :1.000  
#>  1st Qu.:118.2   1st Qu.: 184.6   1st Qu.:  82.15   1st Qu.:2.000  
#>  Median :235.5   Median : 424.0   Median : 319.30   Median :2.000  
#>  Mean   :235.5   Mean   : 435.3   Mean   : 604.08   Mean   :1.904  
#>  3rd Qu.:352.8   3rd Qu.: 640.9   3rd Qu.: 844.55   3rd Qu.:2.000  
#>  Max.   :470.0   Max.   :1528.1   Max.   :5190.10   Max.   :2.000  
#>                                   NAs    :195                      
summary(walker.exh)
#> Object of class SpatialGridDataFrame
#> Coordinates:
#>   min   max
#> X 0.5 260.5
#> Y 0.5 300.5
#> Is projected: NA 
#> proj4string : [NA]
#> Grid attributes:
#>   cellcentre.offset cellsize cells.dim
#> X                 1        1       260
#> Y                 1        1       300
#> Data attributes:
#>        U                  V         
#>  Min.   :   0.000   Min.   :   0.0  
#>  1st Qu.:   6.674   1st Qu.:  67.8  
#>  Median :  56.902   Median : 221.2  
#>  Mean   : 266.044   Mean   : 278.0  
#>  3rd Qu.: 316.351   3rd Qu.: 429.3  
#>  Max.   :9499.508   Max.   :1631.2  
```
