# $Id: plot.variogramCloud.q,v 1.7 2007-10-18 10:13:13 edzer Exp $

"plot.variogramCloud" <-
function (x, identify = FALSE, digitize = FALSE, 
	xlim = c(0, max(x$dist)), ylim, # = c(0, max(x$gamma)), 
	xlab = "distance", ylab = "semivariance", keep = FALSE, ...) 
{
    if (identify || digitize) {
		if (missing(ylim)) ylim = c(0, max(x$gamma))
        plot(x$dist, x$gamma, xlim = xlim, ylim = ylim, xlab = xlab, 
            ylab = ylab, ...)
		.BigInt = attr(x, ".BigInt")
        head = floor(x$np %/% .BigInt) + 1
        tail = floor(x$np %%  .BigInt) + 1
		if (identify) {
			print("mouse-left identifies, mouse-right or Esc stops")
        	labs = paste(head, tail, sep = ",")
        	sel = identify(x$dist, x$gamma, labs, pos = keep)
			ret = data.frame(cbind(head, tail)[sel,, drop = FALSE])
		} else {
			print("mouse-left digitizes, mouse-right closes polygon")
			poly = locator(n = 512, type = "l")
			if (!is.null(poly))
				sel = point.in.polygon(x$dist, x$gamma, poly$x, poly$y)
			else stop("digitized selection is empty")
			ret = data.frame(cbind(head, tail)[sel == 1,,drop = FALSE])
		}
		class(ret) = c("pointPairs", "data.frame")
        if (keep) {
			if (identify) {
				attr(x, "sel") = sel
				attr(x, "text") = labs[sel$ind]
			} else  # digitize
				attr(x, "poly") = poly
			attr(x, "ppairs") = ret
			return(x)
		} else 
        	return(ret)
	} else {
		sel = attr(x, "sel")
		lab = attr(x, "text")
		poly = attr(x, "poly")
		if (!is.null(sel) && !is.null(lab)) {
			if (missing(ylim)) ylim = c(0, max(x$gamma))
        	plot(x$dist, x$gamma, xlim = xlim, ylim = ylim, xlab = xlab, 
            	ylab = ylab, ...)
			text(x$dist[sel$ind], x$gamma[sel$ind], labels=lab, pos= sel$pos)
		} else if (!is.null(poly)) {
			if (missing(ylim)) ylim = c(0, max(x$gamma))
        	plot(x$dist, x$gamma, xlim = xlim, ylim = ylim, xlab = xlab, 
            	ylab = ylab, ...)
			lines(poly$x, poly$y)
		} else {
        	x$np = rep(1, length(x$gamma))
        	plot.gstatVariogram(x, xlim = xlim, ylim = ylim, xlab = xlab, 
        	    ylab = ylab, ...)
		}
    }
}
