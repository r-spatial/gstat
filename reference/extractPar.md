# Extracting parameters and their names from a spatio-temporal variogram model

All spatio-temporal variogram models have a different set of parameters.
These functions extract the parameters and their names from the
spatio-temporal variogram model. Note, this function is as well used to
pass the parameters to the optim function. The arguments lower and upper
passed to optim should follow the same structure.

## Usage

``` r
extractPar(model)
extractParNames(model)
```

## Arguments

- model:

  a spatio-temporal variogram model from [`vgmST`](vgmST.md)

## Value

A named numeric vector of parameters or a vector of characters holding
the parameters' names.

## Author

Benedikt Graeler

## See also

[`fit.StVariogram`](fit.StVariogram.md) and [`vgmST`](vgmST.md)

## Examples

``` r
sumMetricModel <- vgmST("sumMetric",
                        space=vgm(30, "Sph", 200,  6),
                        time =vgm(30, "Sph",  15,  7),
                        joint=vgm(60, "Exp",  84, 22),
                        stAni=100)

extractPar(sumMetricModel)
#>    sill.s   range.s  nugget.s    sill.t   range.t  nugget.t   sill.st  range.st 
#>        30       200         6        30        15         7        60        84 
#> nugget.st      anis 
#>        22       100 
extractParNames(sumMetricModel)
#>  [1] "sill.s"    "range.s"   "nugget.s"  "sill.t"    "range.t"   "nugget.t" 
#>  [7] "sill.st"   "range.st"  "nugget.st" "anis"     
```
