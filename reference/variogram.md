# Calculate Sample or Residual Variogram or Variogram Cloud

Calculates the sample variogram from data, or in case of a linear model
is given, for the residuals, with options for directional, robust, and
pooled variogram, and for irregular distance intervals.

In case spatio-temporal data is provided, the function
[`variogramST`](variogramST.md) is called with a different set of
parameters.

## Usage

``` r
# S3 method for class 'gstat'
variogram(object, ...)
# S3 method for class 'formula'
variogram(object, locations = coordinates(data), data, ...)
# Default S3 method
variogram(object, locations, X, cutoff, width = cutoff/15,
  alpha = 0, beta = 0, tol.hor = 90/length(alpha), tol.ver =
  90/length(beta), cressie = FALSE, dX = numeric(0), boundaries =
  numeric(0), cloud = FALSE, trend.beta = NULL, debug.level = 1,
  cross = TRUE, grid, map = FALSE, g = NULL, ..., projected = TRUE, 
  lambda = 1.0, verbose = FALSE, covariogram = FALSE, PR = FALSE, 
  pseudo = -1)
# S3 method for class 'gstatVariogram'
print(x, ...)
# S3 method for class 'variogramCloud'
print(x, ...)
```

## Arguments

- object:

  object of class `gstat`; in this form, direct and cross (residual)
  variograms are calculated for all variables and variable pairs defined
  in `object`; in case of `variogram.formula`, formula defining the
  response vector and (possible) regressors, in case of absence of
  regressors, use e.g. `z~1`; in case of `variogram.default`: list with
  for each variable the vector with responses (should not be called
  directly)

- data:

  data frame where the names in formula are to be found

- locations:

  spatial data locations. For variogram.formula: a formula with only the
  coordinate variables in the right hand (explanatory variable) side
  e.g. `~x+y`; see examples.

  For variogram.default: list with coordinate matrices, each with the
  number of rows matching that of corresponding vectors in y; the number
  of columns should match the number of spatial dimensions spanned by
  the data (1 (x), 2 (x,y) or 3 (x,y,z)).

- ...:

  any other arguments that will be passed to variogram.default (ignored)

- X:

  (optional) list with for each variable the matrix with
  regressors/covariates; the number of rows should match that of the
  correspoding element in y, the number of columns equals the number of
  regressors (including intercept)

- cutoff:

  spatial separation distance up to which point pairs are included in
  semivariance estimates; as a default, the length of the diagonal of
  the box spanning the data is divided by three.

- width:

  the width of subsequent distance intervals into which data point pairs
  are grouped for semivariance estimates

- alpha:

  direction in plane (x,y), in positive degrees clockwise from positive
  y (North): alpha=0 for direction North (increasing y), alpha=90 for
  direction East (increasing x); optional a vector of directions in
  (x,y)

- beta:

  direction in z, in positive degrees up from the (x,y) plane;

optional a vector of directions

- tol.hor:

  horizontal tolerance angle in degrees

- tol.ver:

  vertical tolerance angle in degrees

- cressie:

  logical; if TRUE, use Cressie”s robust variogram estimate; if FALSE
  use the classical method of moments variogram estimate

- dX:

  include a pair of data points \$y(s_1),y(s_2)\$ taken at locations
  \$s_1\$ and \$s_2\$ for sample variogram calculation only when
  \$\|\|x(s_1)-x(s_2)\|\| \< dX\$ with and \$x(s_i)\$ the vector with
  regressors at location \$s_i\$, and \$\|\|.\|\|\$ the 2-norm. This
  allows pooled estimation of within-strata variograms (use a factor
  variable as regressor, and dX=0.5), or variograms of (near-)replicates
  in a linear model (addressing point pairs having similar values for
  regressors variables)

- boundaries:

  numerical vector with distance interval upper boundaries; values
  should be strictly increasing

- cloud:

  logical; if TRUE, calculate the semivariogram cloud

- trend.beta:

  vector with trend coefficients, in case they are known. By default,
  trend coefficients are estimated from the data.

- debug.level:

  integer; set gstat internal debug level

- cross:

  logical or character; if FALSE, no cross variograms are computed when
  object is of class `gstat` and has more than one variable; if TRUE,
  all direct and cross variograms are computed; if equal to "ST", direct
  and cross variograms are computed for all pairs involving the first
  (non-time lagged) variable; if equal to "ONLY", only cross variograms
  are computed (no direct variograms).

- formula:

  formula, specifying the dependent variable and possible covariates

- x:

  object of class `variogram` or `variogramCloud` to be printed

- grid:

  grid parameters, if data are gridded (not to be called directly; this
  is filled automatically)

- map:

  logical; if TRUE, and `cutoff` and `width` are given, a variogram map
  is returned. This requires package sp. Alternatively, a map can be
  passed, of class SpatialDataFrameGrid (see sp docs)

- g:

  NULL or object of class gstat; may be used to pass settable parameters
  and/or variograms; see example

- projected:

  logical; if FALSE, data are assumed to be unprojected, meaning decimal
  longitude/latitude. For projected data, Euclidian distances are
  computed, for unprojected great circle distances (km). In
  `variogram.formula` or `variogram.gstat`, for data deriving from class
  Spatial, projection is detected automatically using `is.projected`

- lambda:

  test feature; not working (yet)

- verbose:

  logical; print some progress indication

- pseudo:

  integer; use pseudo cross variogram for computing time-lagged spatial
  variograms? -1: find out from coordinates – if they are equal then
  yes, else no; 0: no; 1: yes.

- covariogram:

  logical; compute covariogram instead of variogram?

- PR:

  logical; compute pairwise relative variogram (does NOT check whether
  variable is strictly positive)

## Value

If map is TRUE (or a map is passed), a grid map is returned containing
the (cross) variogram map(s). See package sp.

In other cases, an object of class "gstatVariogram" with the following
fields:

- np:

  the number of point pairs for this estimate; in case of a
  `variogramCloud` see below

- dist:

  the average distance of all point pairs considered for this estimate

- gamma:

  the actual sample variogram estimate

- dir.hor:

  the horizontal direction

- dir.ver:

  the vertical direction

- id:

  the combined id pair

If cloud is TRUE: an object of class `variogramCloud`, with the field
`np` encoding the numbers of the point pair that contributed to a
variogram cloud estimate, as follows. The first point is found by 1 +
the integer division of np by the `.BigInt` attribute of the returned
object, the second point by 1 + the remainder of that division.
as.data.frame.variogramCloud returns no `np` field, but does the
decoding into:

- left:

  for variogramCloud: data id (row number) of one of the data pair

- right:

  for variogramCloud: data id (row number) of the other data in the pair

In case of a spatio-temporal variogram is sought see
[`variogramST`](variogramST.md) for details.

## Note

`variogram.default` should not be called by users directly, as it makes
many assumptions about the organization of the data, that are not fully
documented (but of course, can be understood from reading the source
code of the other `variogram` methods)

Successfully setting `gridded() <- TRUE` may trigger a branch that will
fail unless dx and dy are identical, and not merely similar to within
machine epsilon.

## References

Cressie, N.A.C., 1993, Statistics for Spatial Data, Wiley.

Cressie, N., C. Wikle, 2011, Statistics for Spatio-temporal Data, Wiley.

Pebesma, E.J., 2004. Multivariable geostatistics in S: the gstat
package. Computers and Geosciences, 30: 683-691.

## Author

Edzer Pebesma

## See also

print.gstatVariogram, [plot.gstatVariogram](plot.gstatVariogram.md),
[plot.variogramCloud](plot.variogramCloud.md); for variogram models:
[vgm](vgm.md), to fit a variogram model to a sample variogram:
[fit.variogram](fit.variogram.md) [`variogramST`](variogramST.md) for
details on the spatio-temporal sample variogram.

## Examples

``` r
library(sp)
data(meuse)
# no trend:
coordinates(meuse) = ~x+y
variogram(log(zinc)~1, meuse)
#>     np       dist     gamma dir.hor dir.ver   id
#> 1   57   79.29244 0.1234479       0       0 var1
#> 2  299  163.97367 0.2162185       0       0 var1
#> 3  419  267.36483 0.3027859       0       0 var1
#> 4  457  372.73542 0.4121448       0       0 var1
#> 5  547  478.47670 0.4634128       0       0 var1
#> 6  533  585.34058 0.5646933       0       0 var1
#> 7  574  693.14526 0.5689683       0       0 var1
#> 8  564  796.18365 0.6186769       0       0 var1
#> 9  589  903.14650 0.6471479       0       0 var1
#> 10 543 1011.29177 0.6915705       0       0 var1
#> 11 500 1117.86235 0.7033984       0       0 var1
#> 12 477 1221.32810 0.6038770       0       0 var1
#> 13 452 1329.16407 0.6517158       0       0 var1
#> 14 457 1437.25620 0.5665318       0       0 var1
#> 15 415 1543.20248 0.5748227       0       0 var1
# residual variogram w.r.t. a linear trend:
variogram(log(zinc)~x+y, meuse)
#>     np       dist     gamma dir.hor dir.ver   id
#> 1   57   79.29244 0.1060834       0       0 var1
#> 2  299  163.97367 0.1829983       0       0 var1
#> 3  419  267.36483 0.2264256       0       0 var1
#> 4  457  372.73542 0.2847192       0       0 var1
#> 5  547  478.47670 0.3162418       0       0 var1
#> 6  533  585.34058 0.3571578       0       0 var1
#> 7  574  693.14526 0.3701742       0       0 var1
#> 8  564  796.18365 0.4201289       0       0 var1
#> 9  589  903.14650 0.4216983       0       0 var1
#> 10 543 1011.29177 0.4772549       0       0 var1
#> 11 500 1117.86235 0.5075874       0       0 var1
#> 12 477 1221.32810 0.4617632       0       0 var1
#> 13 452 1329.16407 0.5512305       0       0 var1
#> 14 457 1437.25620 0.4352155       0       0 var1
#> 15 415 1543.20248 0.4556815       0       0 var1
# directional variogram:
variogram(log(zinc)~x+y, meuse, alpha=c(0,45,90,135))
#>     np       dist      gamma dir.hor dir.ver   id
#> 1   12   84.36080 0.04114593       0       0 var1
#> 2   76  165.59800 0.19091543       0       0 var1
#> 3  109  270.29441 0.21867508       0       0 var1
#> 4  134  371.27824 0.23112878       0       0 var1
#> 5  158  478.06480 0.38337565       0       0 var1
#> 6  154  583.35601 0.35513567       0       0 var1
#> 7  159  692.50911 0.35709265       0       0 var1
#> 8  158  797.52941 0.46221222       0       0 var1
#> 9  156  901.86529 0.47081724       0       0 var1
#> 10 156 1011.55318 0.50937290       0       0 var1
#> 11 137 1115.24492 0.57358764       0       0 var1
#> 12 135 1220.31674 0.43193998       0       0 var1
#> 13 109 1328.07859 0.68882673       0       0 var1
#> 14 120 1436.93237 0.53015452       0       0 var1
#> 15  96 1544.68559 0.66909962       0       0 var1
#> 16  11   82.06663 0.07619858      45       0 var1
#> 17  91  165.75829 0.11957011      45       0 var1
#> 18 118  266.93093 0.20557549      45       0 var1
#> 19 136  374.24886 0.27864922      45       0 var1
#> 20 172  479.40618 0.23932562      45       0 var1
#> 21 177  587.53554 0.28038440      45       0 var1
#> 22 209  693.02620 0.34028114      45       0 var1
#> 23 226  796.37554 0.37201935      45       0 var1
#> 24 283  905.25038 0.36146985      45       0 var1
#> 25 264 1012.26326 0.36891951      45       0 var1
#> 26 274 1121.20926 0.36831067      45       0 var1
#> 27 275 1221.63704 0.33875319      45       0 var1
#> 28 282 1330.93431 0.33848846      45       0 var1
#> 29 297 1438.21262 0.31476883      45       0 var1
#> 30 299 1542.75515 0.31707228      45       0 var1
#> 31  16   78.75466 0.07583160      90       0 var1
#> 32  70  160.01667 0.20149652      90       0 var1
#> 33  97  267.68973 0.20686187      90       0 var1
#> 34  98  372.02688 0.28167260      90       0 var1
#> 35 118  479.76226 0.30366429      90       0 var1
#> 36  98  585.85589 0.46344817      90       0 var1
#> 37 115  691.04342 0.36401272      90       0 var1
#> 38 100  796.22142 0.36912878      90       0 var1
#> 39  88  901.26201 0.50261434      90       0 var1
#> 40  72 1004.66642 0.56369456      90       0 var1
#> 41  68 1109.43463 0.77219638      90       0 var1
#> 42  51 1223.73294 0.79679699      90       0 var1
#> 43  44 1322.80887 0.82262644      90       0 var1
#> 44  30 1430.99001 0.80073011      90       0 var1
#> 45  16 1544.27842 1.17421050      90       0 var1
#> 46  18   74.69621 0.19452856     135       0 var1
#> 47  62  163.83075 0.24550456     135       0 var1
#> 48  95  264.21071 0.28119200     135       0 var1
#> 49  89  373.39690 0.37803627     135       0 var1
#> 50  99  475.98691 0.35772223     135       0 var1
#> 51 104  584.05805 0.39065627     135       0 var1
#> 52  91  697.18636 0.46947283     135       0 var1
#> 53  80  792.93648 0.53667425     135       0 var1
#> 54  62  899.44175 0.45817366     135       0 var1
#> 55  51 1014.81674 0.81777411     135       0 var1
#> 56  21 1118.55839 1.03741404     135       0 var1
#> 57  16 1216.88607 1.75971197     135       0 var1
#> 58  17 1323.20745 2.49557308     135       0 var1
#> 59  10 1431.53529 1.77666963     135       0 var1
#> 60   4 1536.74264 2.82057119     135       0 var1
variogram(log(zinc)~1, meuse, width=90, cutoff=1300)
#>     np       dist     gamma dir.hor dir.ver   id
#> 1   41   72.24836 0.1404979       0       0 var1
#> 2  212  142.88031 0.1719093       0       0 var1
#> 3  320  227.32202 0.2554929       0       0 var1
#> 4  371  315.85549 0.3469081       0       0 var1
#> 5  423  406.44801 0.4255276       0       0 var1
#> 6  458  496.09401 0.5042025       0       0 var1
#> 7  455  586.78634 0.5650016       0       0 var1
#> 8  466  677.39566 0.5478706       0       0 var1
#> 9  503  764.55712 0.6076682       0       0 var1
#> 10 480  856.69422 0.6852387       0       0 var1
#> 11 468  944.02864 0.6516089       0       0 var1
#> 12 460 1033.62277 0.6797202       0       0 var1
#> 13 422 1125.63214 0.7001957       0       0 var1
#> 14 408 1212.62350 0.6145586       0       0 var1
#> 15 173 1280.65364 0.6213803       0       0 var1

# GLS residual variogram:
v = variogram(log(zinc)~x+y, meuse)
v.fit = fit.variogram(v, vgm(1, "Sph", 700, 1))
v.fit
#>   model      psill    range
#> 1   Nug 0.08234213    0.000
#> 2   Sph 0.38866509 1098.571
set = list(gls=1)
v
#>     np       dist     gamma dir.hor dir.ver   id
#> 1   57   79.29244 0.1060834       0       0 var1
#> 2  299  163.97367 0.1829983       0       0 var1
#> 3  419  267.36483 0.2264256       0       0 var1
#> 4  457  372.73542 0.2847192       0       0 var1
#> 5  547  478.47670 0.3162418       0       0 var1
#> 6  533  585.34058 0.3571578       0       0 var1
#> 7  574  693.14526 0.3701742       0       0 var1
#> 8  564  796.18365 0.4201289       0       0 var1
#> 9  589  903.14650 0.4216983       0       0 var1
#> 10 543 1011.29177 0.4772549       0       0 var1
#> 11 500 1117.86235 0.5075874       0       0 var1
#> 12 477 1221.32810 0.4617632       0       0 var1
#> 13 452 1329.16407 0.5512305       0       0 var1
#> 14 457 1437.25620 0.4352155       0       0 var1
#> 15 415 1543.20248 0.4556815       0       0 var1
g = gstat(NULL, "log-zinc", log(zinc)~x+y, meuse, model=v.fit, set = set)
variogram(g)
#>     np       dist     gamma dir.hor dir.ver       id
#> 1   57   79.29244 0.1059824       0       0 log-zinc
#> 2  299  163.97367 0.1826061       0       0 log-zinc
#> 3  419  267.36483 0.2256105       0       0 log-zinc
#> 4  457  372.73542 0.2839247       0       0 log-zinc
#> 5  547  478.47670 0.3156087       0       0 log-zinc
#> 6  533  585.34058 0.3566519       0       0 log-zinc
#> 7  574  693.14526 0.3686387       0       0 log-zinc
#> 8  564  796.18365 0.4203337       0       0 log-zinc
#> 9  589  903.14650 0.4212182       0       0 log-zinc
#> 10 543 1011.29177 0.4766290       0       0 log-zinc
#> 11 500 1117.86235 0.5089493       0       0 log-zinc
#> 12 477 1221.32810 0.4637839       0       0 log-zinc
#> 13 452 1329.16407 0.5501712       0       0 log-zinc
#> 14 457 1437.25620 0.4388564       0       0 log-zinc
#> 15 415 1543.20248 0.4580371       0       0 log-zinc

if (require(sf)) {
  proj4string(meuse) = CRS("+init=epsg:28992")
  meuse.ll = sf::st_transform(sf::st_as_sf(meuse), sf::st_crs("+proj=longlat +datum=WGS84"))
# variogram of unprojected data, using great-circle distances, returning km as units
  print(variogram(log(zinc) ~ 1, meuse.ll))
}
#> Loading required package: sf
#> Linking to GEOS 3.12.1, GDAL 3.8.4, PROJ 9.4.0; sf_use_s2() is TRUE
#>     np       dist     gamma dir.hor dir.ver   id
#> 1   57 0.07929104 0.1234479       0       0 var1
#> 2  299 0.16397078 0.2162185       0       0 var1
#> 3  419 0.26736014 0.3027859       0       0 var1
#> 4  457 0.37272890 0.4121448       0       0 var1
#> 5  548 0.47856656 0.4626633       0       0 var1
#> 6  533 0.58553009 0.5646904       0       0 var1
#> 7  573 0.69322813 0.5698718       0       0 var1
#> 8  566 0.79636575 0.6169183       0       0 var1
#> 9  587 0.90330623 0.6489406       0       0 var1
#> 10 544 1.01137212 0.6932668       0       0 var1
#> 11 500 1.11805589 0.7009152       0       0 var1
#> 12 479 1.22176407 0.6044300       0       0 var1
#> 13 450 1.32960776 0.6507175       0       0 var1
#> 14 456 1.43734917 0.5675707       0       0 var1
#> 15 415 1.54317675 0.5748227       0       0 var1
```
