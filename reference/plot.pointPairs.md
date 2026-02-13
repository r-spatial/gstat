# Plot a point pairs, identified from a variogram cloud

Plot a point pairs, identified from a variogram cloud

## Usage

``` r
# S3 method for class 'pointPairs'
plot(x, data, xcol = data$x, ycol = data$y, xlab = "x coordinate", 
ylab = "y coordinate", col.line = 2, line.pch = 0, main = "selected point pairs", ...)
```

## Arguments

- x:

  object of class "pointPairs", obtained from the function
  [plot.variogramCloud](plot.variogramCloud.md), containing point pair
  indices

- data:

  data frame to which the indices refer (from which the variogram cloud
  was calculated)

- xcol:

  numeric vector with x-coordinates of data

- ycol:

  numeric vector with y-coordinates of data

- xlab:

  x-axis label

- ylab:

  y-axis label

- col.line:

  color for lines connecting points

- line.pch:

  if non-zero, symbols are also plotted at the middle of line segments,
  to mark lines too short to be visible on the plot; the color used is
  `col.line`; the value passed to this argument will be used as plotting
  symbol (pch)

- main:

  title of plot

- ...:

  arguments, further passed to `xyplot`

## Value

plots the data locations, with lines connecting the point pairs
identified (and refered to by indices in) x

## Author

Edzer Pebesma

## See also

[plot.variogramCloud](plot.variogramCloud.md)

## Examples

``` r
### The following requires interaction, and is therefore outcommented
#data(meuse)
#coordinates(meuse) = ~x+y
#vgm1 <- variogram(log(zinc)~1, meuse, cloud = TRUE)
#pp <- plot(vgm1, id = TRUE)
### Identify the point pairs
#plot(pp, data = meuse) # meuse has x and y as coordinates
```
