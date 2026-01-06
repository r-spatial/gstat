# Get or set progress indicator

Get or set progress indicator

## Usage

``` r
get_gstat_progress()
set_gstat_progress(value)
```

## Arguments

- value:

  logical

## Value

return the logical value indicating whether progress bars should be
given

## Author

Edzer Pebesma

## Examples

``` r
set_gstat_progress(FALSE)
#> [1] FALSE
get_gstat_progress()
#> [1] FALSE
```
