# Generate, or Add to Variogram Model

Generates a variogram model, or adds to an existing model.
`print.variogramModel` prints the essence of a variogram model.

## Usage

``` r
vgm(psill = NA, model, range = NA, nugget, add.to, anis, kappa = 0.5, ..., covtable,
  Err = 0)
# S3 method for class 'variogramModel'
print(x, ...)
# S3 method for class 'variogramModel'
plot(x, cutoff, ..., type = 'l')
as.vgm.variomodel(m)
```

## Arguments

- psill:

  (partial) sill of the variogram model component, or model: see Details

- model:

  model type, e.g. "Exp", "Sph", "Gau", or "Mat". Can be a character
  vector of model types combined with c(), e.g. c("Exp", "Sph"), in
  which case the best fitting is returned. Calling vgm() without a model
  argument returns a data.frame with available models.

- range:

  range parameter of the variogram model component; in case of
  anisotropy: major range

- kappa:

  smoothness parameter for the Matern class of variogram models

- nugget:

  nugget component of the variogram (this basically adds a nugget
  compontent to the model); if missing, nugget component is omitted

- add.to:

  the variogram model to which we want to add a component (structure)

- anis:

  anisotropy parameters: see notes below

- x:

  a variogram model to print or plot

- ...:

  arguments that will be passed to `print`, e.g. `digits` (see
  examples), or to `variogramLine` for the plot method

- covtable:

  if model is `Tab`, instead of model parameters a one-dimensional
  covariance table can be passed here. See covtable.R in tests
  directory, and example below.

- Err:

  numeric; if larger than zero, the measurement error variance component
  that will not be included to the kriging equations, i.e. kriging will
  now smooth the process Y instead of predict the measured Z, where
  Z=Y+e, and Err is the variance of e

- m:

  object of class `variomodel`, see geoR

- cutoff:

  maximum distance up to which variogram values are computed

- type:

  plot type

## Value

If a single model is passed, an object of class `variogramModel`
extending `data.frame`.

In case a vector ofmodels is passed, an object of class
`variogramModelList` which is a list of `variogramModel` objects.

When called without a model argument, a data.frame with available models
is returned, having two columns: short (abbreviated names, to be used as
model argument: "Exp", "Sph" etc) and long (with some description).

as.vgm.variomodel tries to convert an object of class variomodel (geoR)
to vgm.

## Author

Edzer Pebesma

## Details

If only the first argument (`psill`) is given a `character` value/vector
indicating one or more models, as in `vgm("Sph")`, then this taken as a
shorthand form of `vgm(NA,"Sph",NA,NA)`, i.e. a spherical variogram with
nugget and unknown parameter values; see examples below. Read
[fit.variogram](fit.variogram.md) to find out how `NA` variogram
parameters are given initial values for a fitting a model, based on the
sample variogram. Package `automap` gives further options for automated
variogram modelling.

## Note

Geometric anisotropy can be modelled for each individual simple model by
giving two or five anisotropy parameters, two for two-dimensional and
five for three-dimensional data. In any case, the range defined is the
range in the direction of the strongest correlation, or the major range.
Anisotropy parameters define which direction this is (the main axis),
and how much shorter the range is in (the) direction(s) perpendicular to
this main axis.

In two dimensions, two parameters define an anisotropy ellipse, say
`anis = c(30, 0.5)`. The first parameter, `30`, refers to the main axis
direction: it is the angle for the principal direction of continuity
(measured in degrees, clockwise from positive Y, i.e. North). The second
parameter, `0.5`, is the anisotropy ratio, the ratio of the minor range
to the major range (a value between 0 and 1). So, in our example, if the
range in the major direction (North-East) is 100, the range in the minor
direction (South-East) is 0.5 x 100 = 50.

In three dimensions, five values should be given in the form
`anis = c(p,q,r,s,t)`. Now, \$p\$ is the angle for the principal
direction of continuity (measured in degrees, clockwise from Y, in
direction of X), \$q\$ is the dip angle for the principal direction of
continuity (measured in positive degrees up from horizontal), \$r\$ is
the third rotation angle to rotate the two minor directions around the
principal direction defined by \$p\$ and \$q\$. A positive angle acts
counter-clockwise while looking in the principal direction. Anisotropy
ratios \$s\$ and \$t\$ are the ratios between the major range and each
of the two minor ranges. The anisotropy code was taken from GSLIB. Note
that in
[http://www.gslib.com/sec_gb.html](http://www.gslib.com/sec_gb.md) it is
reported that this code has a bug. Quoting from this site: “The third
angle in all GSLIB programs operates in the opposite direction than
specified in the GSLIB book. Explanation - The books says (pp27) the
angle is measured clockwise when looking toward the origin (from the
postive principal direction), but it should be counter-clockwise. This
is a documentation error. Although rarely used, the correct
specification of the third angle is critical if used.”

(Note that `anis = c(p,s)` is equivalent to `anis = c(p,0,0,s,1)`.)

The implementation in gstat for 2D and 3D anisotropy was taken from the
gslib (probably 1992) code. I have seen a paper where it is argued that
the 3D anisotropy code implemented in gslib (and so in gstat) is in
error, but I have not corrected anything afterwards.

## References

Pebesma, E.J., 2004. Multivariable geostatistics in S: the gstat
package. Computers and Geosciences, 30: 683-691.

Deutsch, C.V. and Journel, A.G., 1998. GSLIB: Geostatistical software
library and user's guide, second edition, Oxford University Press.

For the validity of variogram models on the sphere, see Huang, Chunfeng,
Haimeng Zhang, and Scott M. Robeson. On the validity of commonly used
covariance and variogram functions on the sphere. Mathematical
Geosciences 43.6 (2011): 721-733.

## See also

[show.vgms](show.vgms.md) to view the available models,
[fit.variogram](fit.variogram.md), [variogramLine](variogramLine.md),
[variogram](variogram.md) for the sample variogram.

## Examples

``` r
vgm()
#>    short                                      long
#> 1    Nug                              Nug (nugget)
#> 2    Exp                         Exp (exponential)
#> 3    Sph                           Sph (spherical)
#> 4    Gau                            Gau (gaussian)
#> 5    Exc        Exclass (Exponential class/stable)
#> 6    Mat                              Mat (Matern)
#> 7    Ste Mat (Matern, M. Stein's parameterization)
#> 8    Cir                            Cir (circular)
#> 9    Lin                              Lin (linear)
#> 10   Bes                              Bes (bessel)
#> 11   Pen                      Pen (pentaspherical)
#> 12   Per                            Per (periodic)
#> 13   Wav                                Wav (wave)
#> 14   Hol                                Hol (hole)
#> 15   Log                         Log (logarithmic)
#> 16   Pow                               Pow (power)
#> 17   Spl                              Spl (spline)
#> 18   Leg                            Leg (Legendre)
#> 19   Err                   Err (Measurement error)
#> 20   Int                           Int (Intercept)
vgm("Sph")
#>   model psill range
#> 1   Nug    NA     0
#> 2   Sph    NA    NA
vgm(NA, "Sph", NA, NA)
#>   model psill range
#> 1   Nug    NA     0
#> 2   Sph    NA    NA
vgm(, "Sph") # "Sph" is second argument: NO nugget in this case
#>   model psill range
#> 1   Sph    NA    NA
vgm(10, "Exp", 300)
#>   model psill range
#> 1   Exp    10   300
x <- vgm(10, "Exp", 300)
vgm(10, "Nug", 0)
#>   model psill range
#> 1   Nug    10     0
vgm(10, "Exp", 300, 4.5)
#>   model psill range
#> 1   Nug   4.5     0
#> 2   Exp  10.0   300
vgm(10, "Mat", 300, 4.5, kappa = 0.7)
#>   model psill range kappa
#> 1   Nug   4.5     0   0.0
#> 2   Mat  10.0   300   0.7
vgm( 5, "Exp", 300, add.to = vgm(5, "Exp", 60, nugget = 2.5))
#>   model psill range
#> 1   Nug   2.5     0
#> 2   Exp   5.0    60
#> 3   Exp   5.0   300
vgm(10, "Exp", 300, anis = c(30, 0.5))
#>   model psill range ang1 anis1
#> 1   Exp    10   300   30   0.5
vgm(10, "Exp", 300, anis = c(30, 10, 0, 0.5, 0.3))
#>   model psill range ang1 ang2 ang3 anis1 anis2
#> 1   Exp    10   300   30   10    0   0.5   0.3
# Matern variogram model:
vgm(1, "Mat", 1, kappa=.3)
#>   model psill range kappa
#> 1   Mat     1     1   0.3
x <- vgm(0.39527463, "Sph", 953.8942, nugget = 0.06105141)
x
#>   model      psill    range
#> 1   Nug 0.06105141   0.0000
#> 2   Sph 0.39527463 953.8942
print(x, digits = 3);
#>   model  psill range
#> 1   Nug 0.0611     0
#> 2   Sph 0.3953   954
# to see all components, do
print.data.frame(x)
#>   model      psill    range kappa ang1 ang2 ang3 anis1 anis2
#> 1   Nug 0.06105141   0.0000   0.0    0    0    0     1     1
#> 2   Sph 0.39527463 953.8942   0.5    0    0    0     1     1
vv=vgm(model = "Tab",  covtable = 
  variogramLine(vgm(1, "Sph", 1), 1, n=1e4, min = 0, covariance = TRUE))
vgm(c("Mat", "Sph"))
#> [[1]]
#>   model psill range kappa
#> 1   Nug    NA     0   0.0
#> 2   Mat    NA    NA   0.5
#> 
#> [[2]]
#>   model psill range
#> 1   Nug    NA     0
#> 2   Sph    NA    NA
#> 
#> attr(,"class")
#> [1] "variogramModelList" "list"              
vgm(, c("Mat", "Sph")) # no nugget
#> [[1]]
#>   model psill range kappa
#> 1   Mat    NA    NA   0.5
#> 
#> [[2]]
#>   model psill range
#> 1   Sph    NA    NA
#> 
#> attr(,"class")
#> [1] "variogramModelList" "list"              
```
