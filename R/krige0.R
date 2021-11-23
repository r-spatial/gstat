extractFormula = function(formula, data, newdata) {
  # extract y and X from data:
  m = model.frame(terms(formula), as(data, "data.frame"), na.action = na.fail)
  y = model.extract(m, "response")
  if (length(y) == 0)
    stop("no response variable present in formula")
  Terms = attr(m, "terms")
  X = model.matrix(Terms, m)
  
  # extract x0 from newdata:
  terms.f = delete.response(terms(formula))
  mf.f = model.frame(terms.f, newdata) #, na.action = na.action)
  x0 = model.matrix(terms.f, mf.f)
  list(y = y, X = X, x0 = x0)
}

idw0 = function(formula, data, newdata, y, idp = 2.0) {
  s = coordinates(data)
  s0 = coordinates(newdata)
  if (missing(y))
    y = extractFormula(formula, data, newdata)$y
  D = 1.0 / (spDists(s0, s) ^ idp)
  sumD = apply(D, 1, sum)
  D %*% y / sumD
}

CHsolve = function(A, b) {
  # solves A x = b for x if A is PD symmetric
  #A = chol(A, LINPACK=TRUE) -> deprecated
  A = chol(A) # but use pivot=TRUE?
  backsolve(A, forwardsolve(A, b, upper.tri = TRUE, transpose = TRUE))
}

krige0 <- function(formula, data, newdata, model, beta, y, ..., 
                   computeVar = FALSE, fullCovariance = FALSE) {
  
  if (inherits(data, "ST"))
  	stopifnot(identical(data@sp@proj4string@projargs, newdata@sp@proj4string@projargs))
  else
  	stopifnot(identical(data@proj4string@projargs, newdata@proj4string@projargs))
  lst = extractFormula(formula, data, newdata)
  X = lst$X
  x0 = lst$x0
  if (missing(y))
    y = lst$y
  ll = (!is.na(is.projected(data)) && !is.projected(data))
  s = coordinates(data)
  s0 = coordinates(newdata)
  if (is(model, "variogramModel")) {
    require(gstat)
    V = variogramLine(model, dist_vector = spDists(s, s, ll),
                      covariance = TRUE)
    v0 = variogramLine(model, dist_vector = spDists(s, s0, ll),
                       covariance = TRUE)
    c0 = variogramLine(model, dist_vector = c(0), covariance = TRUE)$gamma
  } else {
    V = model(data, data, ...)
    v0 = model(data, newdata, ...)
    if (computeVar) {
      if (is(newdata, "SpatialLines") || is(newdata, "SpatialPolygons"))
        stop("varying target support (SpatialLines, SpatialPolygons) is not implemented")
      c0 = as.numeric(model(newdata[1, drop=FALSE],
                            newdata[1, drop=FALSE]))
      # ?check this: provide TWO arguments, so model(x,y) can target
      # eventually Y, instead of measurements Z=Y+e
      # with e measurement error term e
    }
  }
  if (!missing(beta)) { # sk:
    skwts = CHsolve(V, v0)
    if (computeVar)
      var <- c0 - apply(v0*skwts, 2, sum)
  } else { # ok/uk -- need to estimate beta:
    skwts = CHsolve(V, cbind(v0, X))
    ViX = skwts[,-(1:nrow(s0))]
    skwts = skwts[,1:nrow(s0)]
    beta = solve(t(X) %*% ViX, t(ViX) %*% y)
    if (computeVar) { 
      Q = t(x0) - t(ViX) %*% v0
      var <- c0 - apply(v0*skwts, 2, sum) + 
        apply(Q * CHsolve(t(X) %*% ViX, Q), 2, sum)
    }
  }
  pred = x0 %*% beta + t(skwts) %*% (y - X %*% beta)
  if (computeVar) {
    if (fullCovariance) {
      corMat <- cov2cor(variogramLine(model, 
                                      dist_vector = spDists(s0, s0), covariance = TRUE))
      var <- corMat*matrix(sqrt(var) %x% sqrt(var), 
                           nrow(corMat), ncol(corMat))
    }
    list(pred = pred, var = var)
  } else
    pred
}

# define variogram model FUNCTION that can deal with x and y
# being of class SpatialPolygons OR SpatialPoints; SpatialGrid/Pixels are coerced to SpatialPolygons
vgmArea = function(x, y = x, vgm, ndiscr = 16, verbose = FALSE, covariance = TRUE) {
  if (gridded(x))
    x = as(x, "SpatialPolygons")
  if (gridded(y))
    y = as(y, "SpatialPolygons")
  stopifnot(is(x, "SpatialPolygons") || is(x, "SpatialPoints"))
  stopifnot(is(y, "SpatialPolygons") || is(y, "SpatialPoints"))
  stopifnot(is(vgm, "variogramModel"))
  nx = length(x)
  ny = length(y)
  V = matrix(NA, nx, ny)
  if (verbose)
    pb = txtProgressBar(style = 3, max = nx)
  for (i in 1:nx) {
    if (is(x, "SpatialPolygons"))
      px = spsample(x[i,], ndiscr, "regular", offset = c(.5,.5))
    else
      px = x[i,]
    for (j in 1:ny) {
      if (is(y, "SpatialPolygons"))
        py = spsample(y[j,], ndiscr, "regular", offset = c(.5,.5))
      else
        py = y[j,]
      D = spDists(px, py)
      D[D == 0] = 1e-10
      V[i,j] = mean(variogramLine(vgm, dist_vector = D, 
                                  covariance = covariance))
    }
    if (verbose)
      setTxtProgressBar(pb, i)
  }
  if (verbose)
    close(pb)
  V
}
