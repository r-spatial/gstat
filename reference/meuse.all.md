# Meuse river data set – original, full data set

This data set gives locations and top soil heavy metal concentrations
(ppm), along with a number of soil and landscape variables, collected in
a flood plain of the river Meuse, near the village Stein. Heavy metal
concentrations are bulk sampled from an area of approximately 15 m x 15
m.

## Usage

``` r
data(meuse.all)
```

## Format

This data frame contains the following columns:

- sample:

  sample number

- x:

  a numeric vector; x-coordinate (m) in RDM (Dutch topographical map
  coordinates)

- y:

  a numeric vector; y-coordinate (m) in RDM (Dutch topographical map
  coordinates)

- cadmium:

  topsoil cadmium concentration, ppm.; note that zero cadmium values in
  the original data set have been shifted to 0.2 (half the lowest
  non-zero value)

- copper:

  topsoil copper concentration, ppm.

- lead:

  topsoil lead concentration, ppm.

- zinc:

  topsoil zinc concentration, ppm.

- elev:

  relative elevation

- om:

  organic matter, as percentage

- ffreq:

  flooding frequency class

- soil:

  soil type

- lime:

  lime class

- landuse:

  landuse class

- dist.m:

  distance to river Meuse (metres), as obtained during the field survey

- in.pit:

  logical; indicates whether this is a sample taken in a pit

- in.meuse155:

  logical; indicates whether the sample is part of the `meuse` (i.e.,
  filtered) data set; in addition to the samples in a pit, an
  sample (139) with outlying zinc content was removed

- in.BMcD:

  logical; indicates whether the sample is used as part of the subset of
  98 points in the various interpolation examples of Burrough and
  McDonnell

## Author

The actual field data were collected by Ruud van Rijn and Mathieu
Rikken; data compiled for R by Edzer Pebesma

## References

P.A. Burrough, R.A. McDonnell, 1998. Principles of Geographical
Information Systems. Oxford University Press.

## Note

`sample` refers to original sample number. Eight samples were left out
because they were not indicative for the metal content of the soil. They
were taken in an old pit. One sample contains an outlying zinc value,
which was also discarded for the meuse (155) data set.

## See also

[meuse.alt](meuse.alt.md)

## Examples

``` r
data(meuse.all)
summary(meuse.all)
#>      sample             x                y             cadmium      
#>  Min.   :  1.00   Min.   :178605   Min.   :329714   Min.   : 0.000  
#>  1st Qu.: 41.75   1st Qu.:179358   1st Qu.:330771   1st Qu.: 0.800  
#>  Median : 82.50   Median :179945   Median :331558   Median : 1.900  
#>  Mean   : 82.50   Mean   :179989   Mean   :331614   Mean   : 3.109  
#>  3rd Qu.:123.25   3rd Qu.:180626   3rd Qu.:332410   3rd Qu.: 3.725  
#>  Max.   :164.00   Max.   :181390   Max.   :333611   Max.   :18.100  
#>                                                                     
#>      copper            lead             zinc             elev       
#>  Min.   : 14.00   Min.   : 27.00   Min.   : 107.0   Min.   : 0.000  
#>  1st Qu.: 23.00   1st Qu.: 68.75   1st Qu.: 191.8   1st Qu.: 7.390  
#>  Median : 29.50   Median :116.00   Median : 307.5   Median : 8.124  
#>  Mean   : 39.42   Mean   :148.55   Mean   : 464.6   Mean   : 7.775  
#>  3rd Qu.: 48.00   3rd Qu.:201.75   3rd Qu.: 662.5   3rd Qu.: 8.915  
#>  Max.   :128.00   Max.   :654.00   Max.   :1839.0   Max.   :10.520  
#>                                                                     
#>      dist.m             om             ffreq            soil      
#>  Min.   :  10.0   Min.   : 1.000   Min.   :1.000   Min.   :1.000  
#>  1st Qu.:  80.0   1st Qu.: 5.000   1st Qu.:1.000   1st Qu.:1.000  
#>  Median : 270.0   Median : 6.550   Median :1.000   Median :1.000  
#>  Mean   : 294.2   Mean   : 7.291   Mean   :1.604   Mean   :1.463  
#>  3rd Qu.: 450.0   3rd Qu.: 8.950   3rd Qu.:2.000   3rd Qu.:2.000  
#>  Max.   :1000.0   Max.   :17.000   Max.   :3.000   Max.   :3.000  
#>                   NA's   :2                                       
#>       lime           landuse     in.pit        in.meuse155      in.BMcD       
#>  Min.   :0.0000   W      :54   Mode :logical   Mode :logical   Mode :logical  
#>  1st Qu.:0.0000   Ah     :42   FALSE:156       FALSE:9         FALSE:66       
#>  Median :0.0000   Am     :22   TRUE :8         TRUE :155       TRUE :98       
#>  Mean   :0.2988   Fw     :10                                                  
#>  3rd Qu.:1.0000   Ab     : 8                                                  
#>  Max.   :1.0000   (Other):27                                                  
#>                   NA's   : 1                                                  
```
