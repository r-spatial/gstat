# rearrange data frame for plotting with levelplot

rearrange data frame for plotting with levelplot

## Usage

``` r
map.to.lev(data, xcol = 1, ycol = 2, zcol = c(3, 4), ns = names(data)[zcol])
```

## Arguments

- data:

  data frame, e.g. output from [krige](krige.md) or
  [predict](predict.gstat.md)

- xcol:

  x-coordinate column number

- ycol:

  y-coordinate column number

- zcol:

  z-coordinate column number range

- ns:

  names of the set of z-columns to be viewed

## Value

data frame with the following elements:

- x:

  x-coordinate for each row

- y:

  y-coordinate for each row

- z:

  column vector with each of the elements in columns `zcol` of `data`
  stacked

- name:

  factor; name of each of the stacked `z` columns

## See also

[image.data.frame](image.md), [krige](krige.md); for examples see
[predict](predict.gstat.md); `levelplot` in package lattice.
