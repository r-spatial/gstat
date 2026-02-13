# Calculate contrasts from multivariable predictions

Given multivariable predictions and prediction (co)variances, calculate
contrasts and their (co)variance

## Usage

``` r
get.contr(data, gstat.object, X, ids = names(gstat.object$data))
```

## Arguments

- data:

  data frame, output of [predict](predict.gstat.md)

- gstat.object:

  object of class `gstat`, used to extract ids; may be missing if `ids`
  is used

- X:

  contrast vector or matrix; the number of variables in `gstat.object`
  should equal the number of elements in `X` if `X` is a vector, or the
  number of rows in `X` if `X` is a matrix.

- ids:

  character vector with (selection of) id names, present in data

## Details

From data, we can extract the \\n \times 1\\ vector with multivariable
predictions, say \$y\$, and its \\n \times n\\ covariance matrix \$V\$.
Given a contrast matrix in \$X\$, this function computes the contrast
vector \$C=X'y\$ and its variance \$Var(C)=X'V X\$.

## Value

a data frame containing for each row in `data` the generalized least
squares estimates (named beta.1, beta.2, ...), their variances (named
var.beta.1, var.beta.2, ...) and covariances (named cov.beta.1.2,
cov.beta.1.3, ...)

## Author

Edzer Pebesma

## See also

[predict](predict.gstat.md)
