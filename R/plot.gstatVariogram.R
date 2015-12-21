# $Id: plot.gstatVariogram.q,v 1.15 2007-06-08 18:03:25 edzer Exp $

"plot.gstatVariogram" <-
function (x, model = NULL, ylim, xlim, xlab = "distance", 
	ylab = attr(x, "what"), 
	panel = vgm.panel.xyplot, multipanel = TRUE, plot.numbers = FALSE, 
	scales = list(), ids = x$id, group.id = TRUE, skip, layout, ...) 
{
    if (missing(ylim)) {
        ylim = c(min(0, 1.04 * min(x$gamma)), 1.04 * max(x$gamma))
		ylim.set = FALSE
	} else
		ylim.set = TRUE
    if (missing(xlim)) 
        xlim = c(0, 1.04 * max(x$dist))
    labels = NULL
	shift = 0.03
	if (is.numeric(plot.numbers)) {
		shift = plot.numbers
		plot.numbers = TRUE
	} 
    if (plot.numbers == TRUE) 
       	labels = as.character(x$np)
    if (length(unique(x$dir.ver)) > 1 || any(x$dir.ver != 0))
		warning("vertical directions are not dealt with -- yet!")
    if (length(unique(x$dir.hor)) > 1 && group.id == TRUE) { 
	# directional, grouped:
        if (multipanel) {
            if (length(levels(ids)) > 1) { # multivariate directional:
				xyplot(gamma ~ dist | as.factor(dir.hor), data = x, 
					type = c("p", "l"), xlim = xlim, ylim = ylim, xlab = xlab, 
					ylab = ylab, groups = ids, scales = scales, ...)
			} else # univariate directional, multipanel:
				xyplot(gamma ~ dist | as.factor(dir.hor), subscripts = TRUE, 
                	panel = panel, data = x, xlim = xlim, 
                	ylim = ylim, xlab = xlab, ylab = ylab, direction = x$dir.hor, 
                	labels = labels, model = model, shift = shift, 
					mode = "directional", scales = scales, ...)
        } else { # univariate directional, using symbol/color to distinguish
            pch = as.integer(as.factor(x$dir.hor))
            xyplot(gamma ~ dist, data = x, type = c("p", "l"), 
                groups = pch, xlim = xlim, ylim = ylim, xlab = xlab, 
                ylab = ylab, pch = pch, scales = scales, ...)
        }
    } else if (length(unique(ids)) > 1) { # multivariable:
        n = floor(sqrt(2 * length(unique(ids))))
		if (missing(skip)) {
        	skip = NULL
        	for (row in n:1) 
				for (col in 1:n) 
					skip = c(skip, row < col)
		}
		if (missing(layout))
			layout = c(n,n)
        if (missing(scales)) 
            scales = list(y = list(relation = "free"))
		else
			if (!is.null(scales$relation) && scales$relation == "same")
				ylim.set = TRUE
    	if (length(unique(x$dir.hor)) > 1) { # multiv.; directional groups
			if (ylim.set) {
            	xyplot(gamma ~ dist | id, data = x, type = c("p", 
                	"l"), xlim = xlim, ylim = ylim, xlab = xlab, 
                	ylab = ylab, groups = as.factor(x$dir.hor), layout = layout,
                  	skip = skip, scales = scales, ...)
			} else {
            	xyplot(gamma ~ dist | id, data = x, type = c("p", 
                	"l"), xlim = xlim, xlab = xlab, 
                	ylab = ylab, groups = as.factor(x$dir.hor), layout = layout,
                  	skip = skip, scales = scales, ...)
			}
		} else { # non-multi-directional, multivariable
			if (ylim.set) {
        		xyplot(gamma ~ dist | id, data = x, xlim = xlim, 
            		ylim = ylim, xlab = xlab, ylab = ylab, ids = ids, 
            		panel= panel, labels = labels, scales = scales, 
            		layout = layout, skip = skip, prepanel = function(x, y) 
					list(ylim = c(min(0, y), max(0, y))), model = model, 
					direction = c(x$dir.hor[1], x$dir.ver[1]), shift = shift, 
					mode = "cross", ...)
			} else {
        		xyplot(gamma ~ dist | id, data = x, xlim = xlim, 
            		xlab = xlab, ylab = ylab, ids = ids, 
            		panel = panel, labels = labels, scales = scales, 
            		layout = layout, skip = skip, prepanel = function(x, 
                		y) list(ylim = c(min(0, y), max(0, y))), 
					model = model, direction = c(x$dir.hor[1], x$dir.ver[1]), 
					shift = shift, mode = "cross", ...)
			}
		}
    } else  # non multi-directional, univariable -- mostly used of all:
		xyplot(gamma ~ dist, data = x, panel = panel, xlim = xlim, 
			ylim = ylim, xlab = xlab, ylab = ylab, labels = labels, model = model, 
			direction = c(x$dir.hor[1], x$dir.ver[1]), shift = shift, 
			mode = "direct", scales = scales, ...)
}

"plot.variogramMap" <-
function(x, np = FALSE, skip, threshold, ...) {
	x = x$map
	if (!is(x, "SpatialPixelsDataFrame"))
		stop("x should be of class, or extend, SpatialPixelsDataFrame")
	if (np)
		start = 2 
	else
		start = 1 
	idx = seq(start, ncol(x@data), by=2)
	n = floor(sqrt(length(idx) * 2))
    if (missing(skip)) {
        skip = NULL
        for (row in n:1)
            for (col in 1:n)
                skip = c(skip, row < col)
    }
	if (!(missing(threshold)))
		x = x[x[[2]] >= threshold, ]

	levelplot(values ~ dx + dy | ind, as.data.frame(stack(x, select = idx)),
		asp = mapasp(x), layout = c(n, n), skip = skip, ...)
}
