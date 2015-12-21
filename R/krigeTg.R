# $Id: krigeTg.q,v 1.4 2009-07-07 15:42:39 edzer Exp $

phiInv <- function (x, lambda)
     if (lambda==0) log(x) else (x^lambda-1)/lambda
  
phi <- function(x, lambda)
     if (lambda==0) exp(x) else (x*lambda+1)^(1/lambda)
  
phiPrime <- function (x, lambda)
     if (lambda==0) exp(x) else (x*lambda+1)^(1/lambda-1)
  
phiDouble <- function (x, lambda)
     if (lambda==0) exp(x) else
        lambda * (1/lambda - 1) * (lambda * x + 1)^(1/lambda-2)

krigeTg <- function(formula, locations, newdata, model = NULL, ..., 
	nmax = Inf, nmin = 0, maxdist = Inf, block = numeric(0), 
	nsim = 0, na.action = na.pass, debug.level = 1,
	lambda = 1.0)
{
    m = model.frame(terms(formula), as(locations, "data.frame"))
    Y = model.extract(m, "response")
    if (length(Y) == 0)
        stop("no response variable present in formula")
    Terms = attr(m, "terms")
    X = model.matrix(Terms, m)
    has.intercept = attr(Terms, "intercept")
	if (ncol(X) > 1)
		stop("only formula with intercept allowed, e.g. y ~ 1")
	locations$value = phiInv(Y, lambda)
	locations$value1 = rep(1, length(locations$value))

    OK = krige(value ~ 1, locations, newdata, model, 
		nmax = nmax, nmin = nmin, 
		maxdist = maxdist, block = block, nsim = nsim,
		na.action = na.action, debug.level = debug.level, ...)

   if (nsim > 0) {
     OK@data = as.data.frame(phi(OK@data,lambda))
     return(OK)
   } 

	# else:
	# estimate mu:
    g = gstat(formula = value ~ 1, # locations = locations, 
		data = locations, model = model, nmax = nmax, nmin = nmin, 
		maxdist = maxdist, ...)
    mu = predict(g, newdata = newdata, block = block, nsim = nsim,
		na.action = na.action, debug.level = debug.level, BLUE = TRUE)
	OK$muhat = mu$var1.pred

    SK = krige(value1 ~ 1, locations, newdata,
		model = model, beta = 0.0, nmax = nmax, nmin = nmin, 
		maxdist = maxdist, block = block, nsim = nsim,
		na.action = na.action, debug.level = debug.level, ...)
	
	# find m:
	OK$m = (OK$var1.var - SK$var1.var)/(1 - SK$var1.pred)

	# copy SK output:
	OK$var1SK.pred = SK$var1.pred
	OK$var1SK.var = SK$var1.var

	# compute transGaussian kriging estimate & variance:
	OK$var1TG.pred = phi(OK$var1.pred, lambda) + 
		phiDouble(mu$var1.pred, lambda) * (OK$var1.var/2 - OK$m)
	OK$var1TG.var = phiPrime(mu$var1.pred, lambda)^2 * OK$var1.var
	OK
}
