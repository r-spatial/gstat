# Coal ash samples from a mine in Pennsylvania

Data obtained from Gomez and Hazen (1970, Tables 19 and 20) on coal ash
for the Robena Mine Property in Greene County Pennsylvania.

## Usage

``` r
data(coalash)
```

## Format

This data frame contains the following columns:

- x:

  a numeric vector; x-coordinate; reference unknown

- y:

  a numeric vector; x-coordinate; reference unknown

- coalash:

  the target variable

## Author

unknown; R version prepared by Edzer Pebesma; data obtained from
<http://homepage.divms.uiowa.edu/~dzimmer/spatialstats/>, Dale
Zimmerman's course page

## References

N.A.C. Cressie, 1993, Statistics for Spatial Data, Wiley.

Gomez, M. and Hazen, K. (1970). Evaluating sulfur and ash distribution
in coal seems by statistical response surface regression analysis. U.S.
Bureau of Mines Report RI 7377.

see also fields manual:
<https://www.image.ucar.edu/GSP/Software/Fields/fields.manual.coalashEX.Krig.shtml>

## Note

data are also present in package fields, as coalash.

## Examples

``` r
data(coalash)
summary(coalash)
#>        x                y            coalash      
#>  Min.   : 1.000   Min.   : 1.00   Min.   : 7.000  
#>  1st Qu.: 5.000   1st Qu.: 8.00   1st Qu.: 8.960  
#>  Median : 7.000   Median :13.00   Median : 9.785  
#>  Mean   : 7.534   Mean   :12.91   Mean   : 9.779  
#>  3rd Qu.:10.000   3rd Qu.:18.00   3rd Qu.:10.568  
#>  Max.   :16.000   Max.   :23.00   Max.   :17.610  
```
