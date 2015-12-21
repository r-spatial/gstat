# $Id: gstat.q,v 1.28 2009-11-02 21:33:17 edzer Exp $

"cross.name" <- function(id1, id2) {
    paste(id1, id2, sep = ".")
}

"gstat" <-
function (g, id, formula, locations,
	data = NULL, model = NULL, beta,
	nmax = Inf, nmin = 0, omax = 0, maxdist = Inf, force = FALSE,
	dummy = FALSE, set, fill.all = FALSE, fill.cross = TRUE,
	variance = "identity", weights = NULL, merge, degree = 0, vdist = FALSE,
	lambda = 1.0) 
{
	call = match.call()
	if (!missing(locations) && inherits(locations, "formula")) {
		if (!is.null(data))
			coordinates(data) = locations
		# locations = NULL
	} else if (missing(data) && !missing(locations) && is(locations, "Spatial")) {
		data = locations
		locations = NULL
	}
	if (fill.all) {
	# fill all variogram models
		if (missing(g) || is.null(model))
			stop("fill.all assumes object g and model are supplied")
        g.names = names(g$data)
		for (i in 1:length(g.names)) {
           	g$model[[paste(g.names[i])]] = model
			if (fill.cross) {
				for (j in (i+1):length(g.names))
           			g$model[[cross.name(g.names[i], g.names[j])]] = model
			}
		}
        return(g)
	} 
    if (!missing(g) && inherits(g, "gstat") && !missing(id) && 
        !missing(model) && missing(formula) && missing(locations)) {
		# here, only direct or cross variogram model is defined
        g.names = names(g$data)
		if (length(id) == 2) {
           	m1 = match(id[1], g.names)
           	m2 = match(id[2], g.names)
           	if (is.na(m1)) 
               	stop("first id does not match available data")
           	if (is.na(m1)) 
               	stop("second id does not match available data")
           	nm = cross.name(g.names[min(m1, m2)], g.names[max(m1, m2)])
        } else if (length(id) == 1) {
			m1 = match(id, g.names)
        	if (is.na(m1)) 
           		stop("id does not match available data")
			nm = g.names[m1]
		} else
			stop("id should have length 1 or 2")
        g$model[[nm]] = model
        return(g)
    } 
	if (!inherits(formula, "formula"))
        stop("argument formula should be of class formula")
    #if (!inherits(locations, "formula") && !has.coordinates(data))
    #	stop("argument locations should be of class formula or matrix")
    if (missing(beta) || is.null(beta)) 
        beta = numeric(0)
	vfn = pmatch(variance, c("identity", "mu", "mu(1-mu)", "mu^2", "mu^3"))
	if (is.na(vfn))
		stop("unknown value for variance function")
	if (vfn > 1 && length(beta) == 0)
		stop("non-identity variance function only allowed if beta is supplied")
    if (missing(g) || is.null(g)) {
        g = list()
        g[["data"]] = list()
        g[["model"]] = list()
    } else if (!dummy && !identical(proj4string(g$data[[1]]$data), proj4string(data)))
		stop("data items in gstat object have different coordinate reference systems")
    if (missing(id)) 
        id = paste("var", length(g$data) + 1, sep = "")
    g$data[[id]] = list(formula = formula, # locations = locations, 
        data = data, has.intercept = attr(terms(formula), "intercept"),
		beta = beta, nmax = nmax, nmin = nmin, omax = omax,
		maxdist = maxdist, force = force,
		dummy = dummy, vfn = vfn, weights = weights, degree = degree, 
		vdist = vdist, lambda = lambda)
    g$model[[id]] = model
	if (!missing(locations))
		g$locations = locations
    if (!missing(set)) {
        if (!is.list(set)) 
            stop("argument set should be a list")
        g$set = set
    }
	if (!missing(merge))
		g$merge = merge
	g$call = call
    class(g) = c("gstat", "list")
    g
}

"[.gstat" <- function(x, ids) { 
	if (is.numeric(ids)) {
		if (min(ids) < 1 || max(ids) > length(names(x$data)))
			stop("selection index(es) out of bound")
		ids = names(x$data)[ids]
	} else if (any(is.na(match(ids, names(x$data)))))
		stop("selected ids do not match those of gstat object")
	g = list()
	g$data = x$data[ids]
	if (length(ids) > 1) {
		ids.cross = NULL
		for (i in 2:length(ids))
			for (j in 1:(i-1))
				ids.cross = c(ids.cross, cross.name(ids[j], ids[i]))
		g$model = x$model[c(ids, ids.cross)]
	} else
		g$model = x$model[ids]
	if (!is.null(x$set))
		g$set = x$set
	if (!is.null(g$merge))
		g$merge = x$merge
    class(g) = c("gstat", "list")
	g
}
