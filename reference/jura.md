# Jura data set

The jura data set from Pierre Goovaerts' book (see references below). It
contains four `data.frame`s: prediction.dat, validation.dat and
transect.dat and juragrid.dat, and three `data.frame`s with consistently
coded land use and rock type factors, as well as geographic coordinates.
The examples below show how to transform these into spatial (sp) objects
in a local coordinate system and in geographic coordinates, and how to
transform to metric coordinate reference systems.

## Usage

``` r
data(jura)
```

## Format

The `data.frames` prediction.dat and validation.dat contain the
following fields:

- Xloc:

  X coordinate, local grid km

- Yloc:

  Y coordinate, local grid km

- Landuse:

  see book and below

- Rock:

  see book and below

- Cd:

  mg cadmium \\\mbox{kg}^{-1}\\ topsoil

- Co:

  mg cobalt \\\mbox{kg}^{-1}\\ topsoil

- Cr:

  mg chromium \\\mbox{kg}^{-1}\\ topsoil

- Cu:

  mg copper \\\mbox{kg}^{-1}\\ topsoil

- Ni:

  mg nickel \\\mbox{kg}^{-1}\\ topsoil

- Pb:

  mg lead \\\mbox{kg}^{-1}\\ topsoil

- Zn:

  mg zinc \\\mbox{kg}^{-1}\\ topsoil

The `data.frame` juragrid.dat only has the first four fields. In
addition the `data.frame`s jura.pred, jura.val and jura.grid also have
inserted third and fourth fields giving geographic coordinates:

- long:

  Longitude, WGS84 datum

- lat:

  Latitude, WGS84 datum

## Author

Data preparation by David Rossiter (dgr2@cornell.edu) and Edzer Pebesma;
georeferencing by David Rossiter

## References

Goovaerts, P. 1997. Geostatistics for Natural Resources Evaluation.
Oxford Univ. Press, New-York, 483 p. Appendix C describes (and gives)
the Jura data set.

Atteia, O., Dubois, J.-P., Webster, R., 1994, Geostatistical analysis of
soil contamination in the Swiss Jura: Environmental Pollution 86,
315-327

Webster, R., Atteia, O., Dubois, J.-P., 1994, Coregionalization of trace
metals in the soil in the Swiss Jura: European Journal of Soil Science
45, 205-218

## Note

The points data sets were obtained from
http://home.comcast.net/~pgoovaerts/book.html, which seems to be no
longer available; the grid data were kindly provided by Pierre
Goovaerts.

The following codes were used to convert `prediction.dat` and
`validation.dat` to `jura.pred` and `jura.val` (see examples below):

Rock Types: 1: Argovian, 2: Kimmeridgian, 3: Sequanian, 4: Portlandian,
5: Quaternary.

Land uses: 1: Forest, 2: Pasture (Weide(land), Wiese, Grasland), 3:
Meadow (Wiese, Flur, Matte, Anger), 4: Tillage (Ackerland, bestelltes
Land)

Points 22 and 100 in the validation set (`validation.dat[c(22,100),]`)
seem not to lie exactly on the grid originally intended, but are kept as
such to be consistent with the book.

Georeferencing was based on two control points in the Swiss grid system
shown as Figure 1 of Atteia et al. (see above) and further points
digitized on the tentatively georeferenced scanned map. RMSE 2.4 m.
Location of points in the field was less precise.

## Examples

``` r
data(jura)
summary(prediction.dat)
#>       Xloc            Yloc          Landuse           Rock      
#>  Min.   :0.626   Min.   :0.580   Min.   :1.000   Min.   :1.000  
#>  1st Qu.:2.282   1st Qu.:1.487   1st Qu.:2.000   1st Qu.:2.000  
#>  Median :3.043   Median :2.581   Median :3.000   Median :2.000  
#>  Mean   :2.980   Mean   :2.665   Mean   :2.548   Mean   :2.699  
#>  3rd Qu.:3.665   3rd Qu.:3.752   3rd Qu.:3.000   3rd Qu.:3.000  
#>  Max.   :4.920   Max.   :5.690   Max.   :4.000   Max.   :5.000  
#>        Cd               Co               Cr              Cu        
#>  Min.   :0.1350   Min.   : 1.552   Min.   : 8.72   Min.   :  3.96  
#>  1st Qu.:0.6375   1st Qu.: 6.520   1st Qu.:27.44   1st Qu.: 11.02  
#>  Median :1.0700   Median : 9.760   Median :34.84   Median : 17.60  
#>  Mean   :1.3091   Mean   : 9.303   Mean   :35.07   Mean   : 23.73  
#>  3rd Qu.:1.7150   3rd Qu.:11.980   3rd Qu.:42.22   3rd Qu.: 27.82  
#>  Max.   :5.1290   Max.   :17.720   Max.   :67.60   Max.   :166.40  
#>        Ni              Pb               Zn        
#>  Min.   : 4.20   Min.   : 18.96   Min.   : 25.20  
#>  1st Qu.:13.80   1st Qu.: 36.52   1st Qu.: 55.00  
#>  Median :20.56   Median : 46.40   Median : 73.56  
#>  Mean   :19.73   Mean   : 53.92   Mean   : 75.08  
#>  3rd Qu.:25.42   3rd Qu.: 60.40   3rd Qu.: 89.92  
#>  Max.   :53.20   Max.   :229.56   Max.   :219.32  
summary(validation.dat)
#>       Xloc            Yloc          Landuse          Rock            Cd        
#>  Min.   :0.491   Min.   :0.524   Min.   :1.00   Min.   :1.00   Min.   :0.3250  
#>  1st Qu.:2.207   1st Qu.:1.593   1st Qu.:2.00   1st Qu.:2.00   1st Qu.:0.6765  
#>  Median :3.001   Median :2.389   Median :3.00   Median :2.00   Median :1.1865  
#>  Mean   :2.921   Mean   :2.546   Mean   :2.41   Mean   :2.36   Mean   :1.2343  
#>  3rd Qu.:3.716   3rd Qu.:3.339   3rd Qu.:3.00   3rd Qu.:3.00   3rd Qu.:1.6350  
#>  Max.   :4.745   Max.   :5.285   Max.   :4.00   Max.   :5.00   Max.   :3.7800  
#>        Co               Cr              Cu                Ni       
#>  Min.   : 1.652   Min.   : 3.32   Min.   :  3.552   Min.   : 1.98  
#>  1st Qu.: 7.950   1st Qu.:28.44   1st Qu.:  9.150   1st Qu.:15.28  
#>  Median :10.060   Median :34.54   Median : 16.140   Median :21.28  
#>  Mean   : 9.793   Mean   :34.88   Mean   : 23.218   Mean   :20.76  
#>  3rd Qu.:12.490   3rd Qu.:40.59   3rd Qu.: 23.190   3rd Qu.:25.36  
#>  Max.   :20.600   Max.   :70.00   Max.   :154.600   Max.   :43.68  
#>        Pb               Zn        
#>  Min.   : 18.68   Min.   : 25.00  
#>  1st Qu.: 35.31   1st Qu.: 53.19  
#>  Median : 47.00   Median : 73.92  
#>  Mean   : 56.48   Mean   : 77.96  
#>  3rd Qu.: 60.10   3rd Qu.: 90.40  
#>  Max.   :300.00   Max.   :259.84  
summary(transect.dat)
#>        X           Rock.type        Block.Ni            Cd       
#>  Min.   :1.000   Min.   :1.000   Min.   : 6.611   Min.   :0.135  
#>  1st Qu.:2.312   1st Qu.:1.000   1st Qu.:17.938   1st Qu.:0.655  
#>  Median :3.625   Median :2.000   Median :20.242   Median :1.317  
#>  Mean   :3.625   Mean   :2.047   Mean   :19.988   Mean   :1.486  
#>  3rd Qu.:4.938   3rd Qu.:2.000   3rd Qu.:23.565   3rd Qu.:1.961  
#>  Max.   :6.250   Max.   :4.000   Max.   :37.047   Max.   :3.925  
#>                                                   NA's   :96     
#>        Ni       
#>  Min.   : 4.20  
#>  1st Qu.:13.31  
#>  Median :20.52  
#>  Mean   :19.62  
#>  3rd Qu.:24.98  
#>  Max.   :43.68  
#>  NA's   :90     
summary(juragrid.dat)
#>       Xloc            Yloc          Landuse               Rock     
#>  Min.   :0.300   Min.   :0.100   Forest : 986   Argovian    :1185  
#>  1st Qu.:2.050   1st Qu.:1.550   Pasture:1553   Kimmeridgian:2036  
#>  Median :3.000   Median :2.450   Meadow :3247   Sequanian   :1628  
#>  Mean   :2.884   Mean   :2.558   Tillage: 171   Portlandian : 316  
#>  3rd Qu.:3.750   3rd Qu.:3.400                  Quaternary  : 792  
#>  Max.   :5.100   Max.   :5.900                                     

# the following commands were used to create objects with factors instead
# of the integer codes for Landuse and Rock:
if (FALSE) { # \dontrun{
  jura.pred = prediction.dat
  jura.val = validation.dat
  jura.grid = juragrid.dat

  jura.pred$Landuse = factor(prediction.dat$Landuse, 
  labels=levels(juragrid.dat$Landuse))
  jura.pred$Rock = factor(prediction.dat$Rock, 
  labels=levels(juragrid.dat$Rock))
  jura.val$Landuse = factor(validation.dat$Landuse, 
  labels=levels(juragrid.dat$Landuse))
  jura.val$Rock = factor(validation.dat$Rock, 
  labels=levels(juragrid.dat$Rock))
} # }

# the following commands convert data.frame objects into spatial (sp) objects
#   in the local grid:
require(sp)
coordinates(jura.pred) = ~Xloc+Yloc
coordinates(jura.val) = ~Xloc+Yloc
coordinates(jura.grid) = ~Xloc+Yloc
gridded(jura.grid) = TRUE

# the following commands convert the data.frame objects into spatial (sp) objects
#   in WGS84 geographic coordinates
# example is given only for jura.pred, do the same for jura.val and jura.grid
# EPSG codes can be found by searching make_EPSG()
jura.pred <- as.data.frame(jura.pred)
coordinates(jura.pred) = ~ long + lat
proj4string(jura.pred) = CRS("+init=epsg:4326")
#> Warning: GDAL Message 1: +init=epsg:XXXX syntax is deprecated. It might return a CRS with a non-EPSG compliant axis order.
```
