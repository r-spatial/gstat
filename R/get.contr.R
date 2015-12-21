# $Id: get.contr.q,v 1.8 2006-03-20 15:18:14 edzer Exp $

"get.contr" <-
function (data, gstat.object, X, ids = names(gstat.object$data)) 
{
	contr.fun <- function(x, n, pr.idx, cov.idx, contr) {
		y = matrix(x[pr.idx], n, 1)
		V = matrix(x[cov.idx], n, n)
		beta = t(contr) %*% y
		Vbeta = t(contr) %*% V %*% contr
		ret = c(beta, diag(Vbeta))
		for (j in 1:nrow(Vbeta)) {
	       	if (j > 1)
       			for (k in 1:(j - 1))
        			ret = c(ret, Vbeta[j, k])
		}
		ret
	}
    lti <- function(i, j) { # lower triangular matrix index, when repr as array
        mx = max(i, j) - 1
        mn = min(i, j) - 1
        ((mx) * (mx - 1))/2 + mn + 1
    }
    n = length(ids)
    if (!is.matrix(X)) 
        X = as.matrix(X)
    if (n != nrow(X)) 
        stop("length(ids) should equal nrow(X) or length(X)")
    gstat.names = create.gstat.names(ids)
    names.pr = gstat.names[seq(1, 2 * n, 2)]
    names.cov = matrix("", n, n)
    for (i in 1:n) 
		for (j in 1:n) 
			names.cov[i, j] = ifelse(i == j, gstat.names[2 * i], 
				gstat.names[2 * n + lti(i, j)])
	pr.idx = match(names.pr, names(data))
	cov.idx = match(names.cov, names(data))
	if (any(is.na(pr.idx)) || any(is.na(cov.idx)))
		stop("colunn names in data not matched")

	res = data.frame(t(apply(as.data.frame(data)[names(data)], 1, 
		contr.fun, n = n, pr.idx = pr.idx, cov.idx = cov.idx,
		contr = X)))

	col.names = NULL
    for (j in 1:NCOL(X))
    	col.names = c(col.names, paste("beta", j, sep = "."))
    for (j in 1:NCOL(X))
    	col.names = c(col.names, paste("var.beta", j, sep = "."))
	for (j in 1:NCOL(X)) {
		if (j > 1) {
			for (k in 1:(j - 1)) {
				col.names = c(col.names, paste("cov.beta", 
					k, j, sep = "."))
			}
		}
	}
    names(res) = col.names
    if (is(data, "data.frame"))
		row.names(res) = row.names(data)
	else if (is(data, "SpatialPolygonsDataFrame")) {
		rownames(res) = sapply(data@polygons, function(x) slot(x, "ID"))
		res = SpatialPolygonsDataFrame(as(data, "SpatialPolygons"), res,
				match.ID = TRUE)
	} else if (is(data, "Spatial")) {
		coordinates(res) = coordinates(data)
		gridded(res) = gridded(data)
	}
    res
}
