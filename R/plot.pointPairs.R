# $Id: plot.pointPairs.q,v 1.4 2006-02-10 19:01:07 edzer Exp $

"plot.pointPairs" <-
function(x, data, xcol = data$x, ycol = data$y, xlab = "x coordinate",
	ylab = "y coordinate", col.line = 2, line.pch = 0, 
	main = "selected point pairs", ...) {

	if (is(data, "SpatialPoints")) {
		cc = coordinates(data) 
		xcol = cc[,1]
		ycol = cc[,2]
		xlab = colnames(cc)[1]
		ylab = colnames(cc)[2]
		asp = mapasp(data)
	} else
		asp = "iso"
	xyplot(ycol ~ xcol, aspect = asp,
		panel = panel.pointPairs, xlab = xlab, ylab = ylab, pairs = x,
		col.line = col.line, line.pch = line.pch, main = main, ...)
}
