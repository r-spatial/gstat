# Spatial Interpolation Comparison 1997 data set: Swiss Rainfall

The text below is copied from the data item at ai-geostats, (link no
longer working).

## Usage

``` r
data(sic97) #
```

## Format

The data frames contain the following columns:

- ID:

  this integer value is the number (unique value) of the monitoring
  station

- rainfall:

  rainfall amount, in 10th of mm

## Note

See the pdf that accompanies the original file for a description of the
data. The .dxf file with the Swiss border is not included here.

## Author

Gregoire Dubois and others.

## Examples

``` r
data(sic97) 
image(demstd)
points(sic_full, pch=1)
points(sic_obs, pch=3)
```
