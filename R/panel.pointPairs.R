# $Id: panel.pointPairs.q,v 1.3 2008-03-10 10:00:10 edzer Exp $

"panel.pointPairs" <-
function (x, y, type = "p", pch = plot.symbol$pch, col, col.line = 
	plot.line$col, col.symbol = plot.symbol$col, lty = plot.line$lty, 
	cex = plot.symbol$cex, lwd = plot.line$lwd, pairs = pairs, 
	line.pch = line.pch, ...) 
{
    x = as.numeric(x)
    y = as.numeric(y)
    if (length(x) > 0) {
        if (!missing(col)) {
            if (missing(col.line)) 
                col.line = col
            if (missing(col.symbol)) 
                col.symbol = col
        }
        plot.symbol = trellis.par.get("plot.symbol")
        plot.line = trellis.par.get("plot.line")
        lpoints(x = x, y = y, cex = cex, col = col.symbol, pch = pch, ...)
        if (!missing(pairs)) {
			for (i in seq(along = pairs[,1])) {
				xx = c(x[pairs[i,1]], x[pairs[i,2]])
				yy = c(y[pairs[i,1]], y[pairs[i,2]])
            	llines(x = xx, y = yy, lty = lty, col = col.line, lwd = lwd)
				if (line.pch > 0)
					lpoints(mean(xx), mean(yy), pch = line.pch, col = col.line)
			}
        }
    }
}
