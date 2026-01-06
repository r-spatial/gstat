# Grid for the NCP, the Dutch part of the North Sea

Gridded data for the NCP (Nederlands Continentaal Plat, the Dutch part
of the North Sea), for a 5 km x 5 km grid; stored as data.frame.

## Usage

``` r
data(ncp.grid)
```

## Format

This data frame contains the following columns:

- x:

  x-coordinate, UTM zone 31

- y:

  y-coordinate, UTM zone 31

- depth:

  sea water depth, m.

- coast:

  distance to the coast of the Netherlands, in km.

- area:

  identifier for administrative sub-areas

## Author

Dutch National Institute for Coastal and Marine Management (RIKZ); data
compiled for R by Edzer Pebesma

## See also

[fulmar](fulmar.md)

## Examples

``` r
data(ncp.grid)
summary(ncp.grid)
#>        x                y               depth           coast       
#>  Min.   :466500   Min.   :5699000   Min.   : 1.00   Min.   :  1.00  
#>  1st Qu.:531500   1st Qu.:5859000   1st Qu.:25.00   1st Qu.: 38.00  
#>  Median :566500   Median :5954000   Median :31.00   Median : 76.00  
#>  Mean   :572625   Mean   :5940147   Mean   :31.79   Mean   : 89.72  
#>  3rd Qu.:606500   3rd Qu.:6024000   3rd Qu.:40.00   3rd Qu.:136.00  
#>  Max.   :736500   Max.   :6129000   Max.   :59.00   Max.   :248.00  
#>       area      
#>  Min.   : 1.00  
#>  1st Qu.: 1.00  
#>  Median : 1.00  
#>  Mean   : 1.87  
#>  3rd Qu.: 2.00  
#>  Max.   :19.00  
```
