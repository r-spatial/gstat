hscat = function(formula, data, breaks, pch = 3, cex = .6, mirror = FALSE, 
		variogram.alpha = 0, as.table = TRUE, ...) {
	stopifnot(!missing(breaks))
	x = variogram(formula, data, cloud = TRUE, cutoff = max(breaks), 
		alpha = variogram.alpha, ...)
	x = as.data.frame(x)
	x$class = cut(x$dist, breaks = breaks)
	y = model.frame(formula, data)[[1]]
	x$xx = y[x$left]
	x$yy = y[x$right]
	if (mirror) 
		x = data.frame(
			xx = c(x$yy, y[x$left]),
			yy = c(x$xx, y[x$right]),
			class = c(x$class, x$class))
	lab = as.character(formula)[2]
	panel = function(x,y,subscripts, ...) {
		xr = c(min(x),max(x))
		llines(xr, xr)
		lpoints(x,y,...)
		ltext(min(x), max(y), paste("r =", signif(cor(x,y),3)), adj=c(0,0.5))
	}
	xyplot(xx~yy|class, x, panel = panel,
		main = "lagged scatterplots", xlab = lab, ylab = lab, 
		as.table = as.table, ...)
}

