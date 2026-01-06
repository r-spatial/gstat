# Fulmaris glacialis data

Airborne counts of Fulmaris glacialis during the Aug/Sept 1998 and 1999
flights on the Dutch (Netherlands) part of the North Sea (NCP,
Nederlands Continentaal Plat).

## Usage

``` r
data(fulmar)
```

## Format

This data frame contains the following columns:

- year:

  year of measurement: 1998 or 1999

- x:

  x-coordinate in UTM zone 31

- y:

  y-coordinate in UTM zone 31

- depth:

  sea water depth, in m

- coast:

  distance to coast of the Netherlands, in km.

- fulmar:

  observed density (number of birds per square km)

## Author

Dutch National Institute for Coastal and Marine Management (RIKZ)

## See also

[ncp.grid](ncp.grid.md)

E.J. Pebesma, R.N.M. Duin, P.A. Burrough, 2005. Mapping Sea Bird
Densities over the North Sea: Spatially Aggregated Estimates and
Temporal Changes. Environmetrics 16, (6), p 573-587.

## Examples

``` r
data(fulmar)
summary(fulmar)
#>       year            x                y               depth     
#>  Min.   :1998   Min.   :476210   Min.   :5694947   Min.   : 1.0  
#>  1st Qu.:1998   1st Qu.:535522   1st Qu.:5806777   1st Qu.:13.0  
#>  Median :1999   Median :568618   Median :5896021   Median :25.0  
#>  Mean   :1999   Mean   :576982   Mean   :5889112   Mean   :23.6  
#>  3rd Qu.:1999   3rd Qu.:604579   3rd Qu.:5946550   3rd Qu.:32.0  
#>  Max.   :1999   Max.   :739042   Max.   :6150942   Max.   :54.0  
#>      coast             fulmar      
#>  Min.   :  1.024   Min.   : 0.000  
#>  1st Qu.:  4.857   1st Qu.: 0.000  
#>  Median : 38.696   Median : 0.000  
#>  Mean   : 56.642   Mean   : 1.005  
#>  3rd Qu.: 89.929   3rd Qu.: 0.000  
#>  Max.   :263.205   Max.   :46.487  
if (FALSE) { # \dontrun{
demo(fulmar)
} # }
```
