# version 2.0-1

* try to fix CRS issue found on OSX

* check for temporal unit in krigeST opposed to covariance function and assign the temporal unit found in krigeST to the spatio-temporal variogram for consistency in case none has been provided by the model; add warning when ST variogram doesn't have a `temporal unit` attribute; #42

# version 2.0-0

* add plot.variogramModel; #40

* support `sf` and `stars` objects, for vector and raster data, respectively; #39

# pre 2.0-0 version

* see [ChangeLog](https://github.com/r-spatial/gstat/blob/master/inst/ChangeLog)
