# $Id: vgm.panel.q,v 1.7 2007-06-08 06:45:52 edzer Exp $

"get.direction.unitv" <- function(alpha, beta) {
	cb = cos(beta)
	c(cb * sin(alpha), cb * cos(alpha), sin(beta))
}

"vgm.panel.xyplot" <-
function (x, y, subscripts, type = "p", pch = plot.symbol$pch, 
    col, col.line = plot.line$col, col.symbol = plot.symbol$col, 
    lty = plot.line$lty, cex = plot.symbol$cex, ids, lwd = plot.line$lwd, 
    model = model, direction = direction, labels, shift = shift, mode = mode, ...) 
{
    x <- as.numeric(x)
    y <- as.numeric(y)
    if (length(x) > 0) {
        if (!missing(col)) {
            if (missing(col.line)) 
                col.line <- col
            if (missing(col.symbol)) 
                col.symbol <- col
        }
        plot.symbol <- trellis.par.get("plot.symbol")
        plot.line <- trellis.par.get("plot.line")
        lpoints(x = x, y = y, cex = cex, col = col.symbol, pch = pch, type = type, ...)
        if (!is.null(labels)) 
            ltext(x = x + shift * max(x), y = y, labels = labels[subscripts])

		if (mode == "direct") {
        	if (!missing(model) && !is.null(model)) {
            	ang.hor <- pi * (direction[1]/180)
				ang.ver <- pi * (direction[2]/180)
            	dir <- get.direction.unitv(ang.hor, ang.ver)
            	ret <- variogramLine(model, max(x), dir = dir)
            	llines(x = ret$dist, y = ret$gamma, lty = lty, col = col.line, lwd = lwd)
        	}
		} else if (mode == "cross") {
    		id <- as.character(ids[subscripts][1])
        	if (!missing(model) && !is.null(model)) {
				if (inherits(model, "gstat"))
					m = model$model
				else
					m = model
				if (!is.list(m))
					stop("model argument not of class gstat or list")
				if (is.list(m) && !is.null(m[[id]])) {
                	ang.hor <- pi * (direction[1]/180)
					ang.ver <- pi * (direction[2]/180)
                	dir <- get.direction.unitv(ang.hor, ang.ver)
					ret <- variogramLine(m[[id]], max(x), dir = dir)
					llines(x = ret$dist, y = ret$gamma, lty = lty, col = col.line, lwd = lwd)
				}
			}
        } else if (mode == "directional") {
        	if (!missing(model) && !is.null(model)) {
            	dir <- c(1, 0, 0)
            	if (!missing(direction)) {
                	ang.hor <- pi * (direction[subscripts][1]/180.0)
                	dir <- get.direction.unitv(ang.hor, 0)
            	}
            	ret <- variogramLine(model, max(x), dir = dir)
            	llines(x = ret$dist, y = ret$gamma, lty = lty, col = col.line, lwd = lwd)
        	}
		}

    }
}
