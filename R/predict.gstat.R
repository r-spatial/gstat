# $Id: predict.gstat.q,v 1.35 2009-11-02 21:33:17 edzer Exp $

predict.gstat <-
function (object, newdata, block = numeric(0), nsim = 0, indicators = FALSE,
	BLUE = FALSE, debug.level = 1, mask, na.action = na.pass, sps.args = 
	list(n = 500, type = "regular", offset = c(.5, .5)), ...) 
{
	if (missing(object) || length(object$data) < 1) 
		stop("no data available")
	if (!inherits(object, "gstat"))
		stop("first argument should be of class gstat")
	if (!is.null(object$locations) && inherits(object$locations, "formula") 
			&& !(is(newdata, "Spatial"))) {
		coordinates(newdata) = object$locations
		return.sp = FALSE
	} else
		return.sp = TRUE
	
	max_dist = getMaxDist(object$data, newdata)
	
	.Call(gstat_init, as.integer(debug.level))
	bl_weights = numeric(0)
	if (!missing(mask)) {
		cat("argument mask is deprecated:")
		stop("use a missing value pattern in newdata instead")
	}
	nvars = length(object$data)
	new.X = NULL
	for (i in 1:length(object$data)) {
		name = names(object$data)[i]
		d = object$data[[i]]
		if (!is.null(d$data)) {
			if (!identical(proj4string(d$data), proj4string(newdata)))
				stop(paste(name, ": data item in gstat object and newdata have different coordinate reference systems"))
		}
		if (d$nmax == Inf) 
			nmax = as.integer(-1)
		else nmax = as.integer(d$nmax)
		nmin = as.integer(max(0, d$nmin))
		if (d$maxdist == Inf)
			maxdist = as.numeric(-1)
		else maxdist = d$maxdist
		if (d$dummy) {
			# tr = terms(d$locations)
			if (is.null(d$beta) || length(d$beta) == 0)
				stop("dummy data should have beta defined")
			if (d$degree != 0)
				stop("dummy data cannot have non-zero degree arg; use formula")
			# loc.dim = length(attr(tr, "term.labels"))
			loc.dim = dim(coordinates(newdata))[[2]]
			.Call(gstat_new_dummy_data, as.integer(loc.dim), 
				as.integer(d$has.intercept), as.double(d$beta), 
				nmax, nmin, maxdist, as.integer(d$vfn), 
				as.integer(is.projected(newdata)), as.integer(d$vdist))
			raw = list(xlevels = NULL)
		} else {
			if (is.null(d$weights))
				w = numeric(0)
			else
				w = d$weights
			raw = gstat.formula(d$formula, d$data)
			.Call(gstat_new_data, as.double(raw$y), 
				as.double(raw$locations), as.double(raw$X),
				as.integer(raw$has.intercept), as.double(d$beta),
				nmax, nmin, maxdist, as.integer(d$force), as.integer(d$vfn),
				as.numeric(w), double(0.0), as.integer(d$degree),
				as.integer(is.projected(d$data)), as.integer(d$vdist),
				as.double(d$lambda), as.integer(d$omax))
		}
		if (!is.null(object$model[[name]])) 
			load.variogram.model(object$model[[name]], c(i - 1, i - 1), max_dist = max_dist)
		raw = gstat.formula.predict(d$formula, newdata, na.action = na.action,
			(length(BLUE) == 2 && BLUE[2]), xlev = raw$xlevels)
		if (is.null(new.X)) 
			new.X = raw$X
		else new.X = cbind(new.X, raw$X)
		if (i > 1) {
			for (j in 1:(i - 1)) {
				cross = cross.name(names(object$data)[j], name)
				if (!is.null(object$model[[cross]])) 
					load.variogram.model(object$model[[cross]], 
						c(i - 1, j - 1), max_dist = max_dist)
			}
		}
	}
	if (!is.null(object$set)) 
		gstat.load.set(object$set)
	if (!is.null(object$merge)) 
		gstat.load.merge(object)
	if (is(newdata, "SpatialPolygons")) {
		pol = newdata@polygons
		if (length(pol) != nrow(raw$locations))
			stop("polygons and center points length mismatch")
		block = matrix(NA, 0, 2)
		nd = as(newdata, "SpatialPolygons")
		block.cols = rep(as.numeric(NA), length(pol))
		for (i in seq(along = pol)) {
			sps.args$x = nd[i]
			cc = coordinates(do.call("spsample", sps.args))
			cc[,1] = cc[,1] - raw$locations[i,1]
			cc[,2] = cc[,2] - raw$locations[i,2]
			block.cols[i] = nrow(block) + 1
			block = rbind(block, cc)
		}
		if (length(pol) == 1)
			block.cols = 2
	} else if (is(newdata, "SpatialLines")) {
		lin = newdata@lines
		if (length(lin) != nrow(raw$locations))
			stop("lines and line midpoints length mismatch")
		block = matrix(NA, 0, 2)
		nd = as(newdata, "SpatialLines")
		block.cols = rep(as.numeric(NA), length(lin))
		for (i in seq(along = lin)) {
			sps.args$x = nd[i]
			cc = coordinates(do.call("spsample", sps.args))
			cc[,1] = cc[,1] - raw$locations[i,1]
			cc[,2] = cc[,2] - raw$locations[i,2]
			block.cols[i] = nrow(block) + 1
			block = rbind(block, cc)
		}
		if (length(lin) == 1)
			block.cols = 2
	} else if (!is.null(dim(block))) { # i.e., block is data.frame or matrix
		if (is.data.frame(block) && !is.null(block$weights)) {
			bl_weights = block$weights
			block$weights = NULL
		} 
		block = data.matrix(block) # converts to numeric
		block.cols = ncol(block)
	} else {
		block = as.numeric(block) # make sure it's not integer
		block.cols = numeric(0)
	} 
	# handle NA's in the parts of newdata used:
	valid.pattern = NULL
	if (any(is.na(raw$locations)) || any(is.na(new.X))) {
		valid.pattern = !(apply(cbind(raw$locations, new.X), 1, 
				function(x) any(is.na(x))))
		raw$locations.all = raw$locations
		raw$locations = as.matrix(raw$locations[valid.pattern, ])
		new.X = as.matrix(new.X[valid.pattern, ])
	} 
	if (nsim) {
		if (indicators == TRUE)
			nsim = -abs(nsim)
	# random path: randomly permute row indices
		perm = sample(seq(along = new.X[, 1]))
		ret = .Call(gstat_predict, as.integer(nrow(as.matrix(new.X))),
			as.double(as.vector(raw$locations[perm, ])),
			as.double(as.vector(new.X[perm,])),
			as.integer(block.cols), as.vector(block), as.vector(bl_weights),
			as.integer(nsim), as.integer(BLUE))[[1]]
		if (nsim == 1)
			colsel = seq(1, by=2, length.out=nvars) # pred1 var1 pred2 var2 ...
		else
			colsel = TRUE
		ret = data.frame(cbind(raw$locations, 
			matrix(ret[order(perm), colsel], nrow(as.matrix(new.X)), abs(nsim) * nvars)))
	}
	else {
		ret = .Call(gstat_predict, as.integer(nrow(as.matrix(new.X))),
			as.double(as.vector(raw$locations)), as.vector(new.X), as.integer(block.cols), 
			as.vector(block), as.vector(bl_weights), as.integer(nsim), as.integer(BLUE))[[1]]
		ret = data.frame(cbind(raw$locations, ret))
	}
	.Call(gstat_exit, NULL)
	if (!is.null(valid.pattern) && any(valid.pattern)) {
		ret.all = data.frame(matrix(NA, length(valid.pattern), ncol(ret)))
		ret.all[, 1:ncol(raw$locations.all)] = raw$locations.all
		ret.all[valid.pattern, ] = ret
		ret = ret.all
	}
	if (abs(nsim) > 0) {
		names.vars = names(object$data)
		if (length(names.vars) > 1) 
			names.vars = paste(rep(names.vars, each = abs(nsim)), 
				paste("sim", 1:abs(nsim), sep = ""), sep = ".")
		else
			names.vars = paste("sim", 1:abs(nsim), sep = "")
	} else
		names.vars = create.gstat.names(names(object$data))
	names(ret) = c(dimnames(raw$locations)[[2]], names.vars)

	if (return.sp) {
		if (is(newdata, "SpatialPolygons")) {
			row.names(ret) = sapply(newdata@polygons, function(x) slot(x, "ID"))
			ret = SpatialPolygonsDataFrame(as(newdata, "SpatialPolygons"), ret,
				match.ID = TRUE)
		} else if (is(newdata, "SpatialLines")) {
			row.names(ret) = sapply(newdata@lines, function(x) slot(x, "ID"))
			ret = SpatialLinesDataFrame(as(newdata, "SpatialLines"), ret,
				match.ID = TRUE)
		} else {
			coordinates(ret) = dimnames(raw$locations)[[2]]
			if (gridded(newdata)) {
				returnFullGrid = fullgrid(newdata)
				fullgrid(newdata) = FALSE
				ret = new("SpatialPixelsDataFrame", 
					new("SpatialPixels", as(ret, "SpatialPoints"),
					grid = newdata@grid, grid.index = newdata@grid.index),
        			data = ret@data, coords.nrs = ret@coords.nrs)
				fullgrid(ret) = returnFullGrid
			}
		}
		proj4string(ret) = CRS(proj4string(newdata))
	}

	return(ret)
}

# call with: create.gstat.names(names(object$data))
# creates the names of the output columns in case of (multivariable) prediction
create.gstat.names <- function(ids, names.sep = ".") {
	nvars = length(ids)
	names.vars = character(nvars * 2 + nvars * (nvars - 1)/2)
	pos = 1
	for (i in 1:length(ids)) {
		name = ids[i]
		names.vars[1 + (i - 1) * 2] = paste(name, "pred", sep = names.sep)
		names.vars[2 + (i - 1) * 2] = paste(name, "var", sep = names.sep)
		if (i > 1) {
			for (j in 1:(i - 1)) {
				cross = paste(ids[j], name, sep = names.sep)
				names.vars[nvars * 2 + pos] = paste("cov", cross, 
					sep = names.sep)
				pos = pos + 1
			}
		}
	}
	return(names.vars)
}

getMaxDist = function(dataLst, newdata) {
	d = apply(bbox(newdata), 1, diff)
	if (!is.null(dataLst[[1]]$data)) {
		d2 = apply(bbox(dataLst[[1]]$data), 1, diff)
		d = apply(rbind(d,d2), 2, max) 
		# there are pathetic cases where this would not be sufficient 
	}
	if (length(d) == 2)
		d = c(d, 0)
	stopifnot(length(d) == 3)
	d
}
