# Oxford soil samples

Data: 126 soil augerings on a 100 x 100m square grid, with 6 columns and
21 rows. Grid is oriented with long axis North-north-west to
South-south-east Origin of grid is South-south-east point, 100m outside
grid.

Original data are part of a soil survey carried out by P.A. Burrough in
1967. The survey area is located on the chalk downlands on the Berkshire
Downs in Oxfordshire, UK. Three soil profile units were recognised on
the shallow Rendzina soils; these are Ia - very shallow, grey calcareous
soils less than 40cm deep over chalk; Ct - shallow to moderately deep,
grey-brown calcareous soils on calcareous colluvium, and Cr: deep,
moderately acid, red-brown clayey soils. These soil profile classes were
registered at every augering.

In addition, an independent landscape soil map was made by interpolating
soil boundaries between these soil types, using information from the
changes in landform. Because the soil varies over short distances, this
field mapping caused some soil borings to receive a different
classification from the classification based on the point data.

Also registered at each auger point were the site elevation (m), the
depth to solid chalk rock (in cm) and the depth to lime in cm. Also, the
percent clay content, the Munsell colour components of VALUE and CHROMA
, and the lime content of the soil (as tested using HCl) were recorded
for the top two soil layers (0-20cm and 20-40cm).

Samples of topsoil taken as a bulk sample within a circle of radius 2.5m
around each sample point were used for the laboratory determination of
Mg (ppm), OM1 %, CEC as mequ/100g air dry soil, pH, P as ppm and K
(ppm).

## Usage

``` r
data(oxford)
```

## Format

This data frame contains the following columns:

- PROFILE:

  profile number

- XCOORD:

  x-coordinate, field, non-projected

- YCOORD:

  y-coordinate, field, non-projected

- ELEV:

  elevation, m.

- PROFCLASS:

  soil class, obtained by classifying the soil profile at the sample
  site

- MAPCLASS:

  soil class, obtained by looking up the site location in the soil map

- VAL1:

  Munsell colour component VALUE, 0-20 cm

- CHR1:

  Munsell colour component CHROMA, 20-40 cm

- LIME1:

  Lime content (tested using HCl), 0-20 cm

- VAL2:

  Munsell colour component VALUE, 0-20 cm

- CHR2:

  Munsell colour component CHROMA, 20-40 cm

- LIME2:

  Lime content (tested using HCl), 20-40 cm

- DEPTHCM:

  soil depth, cm

- DEP2LIME:

  depth to lime, cm

- PCLAY1:

  percentage clay, 0-20 cm

- PCLAY2:

  percentage clay, 20-40 cm

- MG1:

  Magnesium content (ppm), 0-20 cm

- OM1:

  organic matter (%), 0-20 cm

- CEC1:

  CES as mequ/100g air dry soil, 0-20 cm

- PH1:

  pH, 0-20 cm

- PHOS1:

  Phosphorous, 0-20 cm, ppm

- POT1:

  K (potassium), 0-20 cm, ppm

## Author

P.A. Burrough; compiled for R by Edzer Pebesma

## References

P.A. Burrough, R.A. McDonnell, 1998. Principles of Geographical
Information Systems. Oxford University Press.

## Note

`oxford.jpg`, in the gstat package external directory (see example
below), shows an image of the soil map for the region

## Examples

``` r
data(oxford)
summary(oxford)
#>     PROFILE           XCOORD        YCOORD          ELEV       PROFCLASS
#>  Min.   :  1.00   Min.   :100   Min.   : 100   Min.   :540.0   Cr:19    
#>  1st Qu.: 32.25   1st Qu.:200   1st Qu.: 600   1st Qu.:558.0   Ct:36    
#>  Median : 63.50   Median :350   Median :1100   Median :573.0   Ia:71    
#>  Mean   : 63.50   Mean   :350   Mean   :1100   Mean   :573.6            
#>  3rd Qu.: 94.75   3rd Qu.:500   3rd Qu.:1600   3rd Qu.:584.5            
#>  Max.   :126.00   Max.   :600   Max.   :2100   Max.   :632.0            
#>  MAPCLASS      VAL1            CHR1           LIME1            VAL2     
#>  Cr:31    Min.   :2.000   Min.   :1.000   Min.   :0.000   Min.   :4.00  
#>  Ct:36    1st Qu.:3.000   1st Qu.:2.000   1st Qu.:1.000   1st Qu.:4.00  
#>  Ia:59    Median :4.000   Median :2.000   Median :4.000   Median :8.00  
#>           Mean   :3.508   Mean   :2.468   Mean   :2.643   Mean   :6.23  
#>           3rd Qu.:4.000   3rd Qu.:3.000   3rd Qu.:4.000   3rd Qu.:8.00  
#>           Max.   :4.000   Max.   :4.000   Max.   :4.000   Max.   :8.00  
#>       CHR2       LIME2          DEPTHCM         DEP2LIME         PCLAY1     
#>  Min.   :2   Min.   :0.000   Min.   :10.00   Min.   :20.00   Min.   :10.00  
#>  1st Qu.:2   1st Qu.:4.000   1st Qu.:25.00   1st Qu.:20.00   1st Qu.:20.00  
#>  Median :2   Median :5.000   Median :36.00   Median :20.00   Median :24.50  
#>  Mean   :3   Mean   :3.889   Mean   :46.25   Mean   :30.32   Mean   :24.44  
#>  3rd Qu.:4   3rd Qu.:5.000   3rd Qu.:64.75   3rd Qu.:40.00   3rd Qu.:28.00  
#>  Max.   :6   Max.   :5.000   Max.   :91.00   Max.   :90.00   Max.   :37.00  
#>      PCLAY2           MG1              OM1              CEC1      
#>  Min.   :10.00   Min.   : 19.00   Min.   : 2.600   Min.   : 7.00  
#>  1st Qu.:10.00   1st Qu.: 44.00   1st Qu.: 4.100   1st Qu.:12.00  
#>  Median :10.00   Median : 72.00   Median : 5.350   Median :15.00  
#>  Mean   :14.76   Mean   : 93.53   Mean   : 5.995   Mean   :18.88  
#>  3rd Qu.:20.00   3rd Qu.:123.25   3rd Qu.: 7.175   3rd Qu.:25.25  
#>  Max.   :40.00   Max.   :308.00   Max.   :13.100   Max.   :43.00  
#>       PH1            PHOS1             POT1      
#>  Min.   :4.200   Min.   : 1.700   Min.   : 83.0  
#>  1st Qu.:7.200   1st Qu.: 6.200   1st Qu.:127.0  
#>  Median :7.500   Median : 8.500   Median :164.0  
#>  Mean   :7.152   Mean   : 8.752   Mean   :181.7  
#>  3rd Qu.:7.600   3rd Qu.:10.500   3rd Qu.:194.8  
#>  Max.   :7.700   Max.   :25.000   Max.   :847.0  
# open the following file with a jpg viewer:
system.file("external/oxford.jpg", package="gstat")
#> [1] "/home/runner/work/_temp/Library/gstat/external/oxford.jpg"
```
