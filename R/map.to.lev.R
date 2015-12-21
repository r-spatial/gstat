# $Id: map.to.lev.q,v 1.2 2006-02-10 19:01:07 edzer Exp $

"map.to.lev" <-
function (data, xcol = 1, ycol = 2, zcol = c(3, 4), ns = names(data)[zcol]) 
{
    len = nrow(data)
    d = matrix(nrow = len * length(zcol), ncol = 3)
    xnames = NULL
    if (length(ns) > 1 && length(ns) != length(zcol)) 
        stop("names should have length 1 or equal to length of zcol")
    nr = 1
    for (i in zcol) {
        if (length(ns) == 1) 
            nm = rep(paste(ns, nr), len)
        else nm = rep(ns[nr], len)
        range = (1 + (nr - 1) * len):(nr * len)
        d[range, ] = cbind(data[, xcol], data[, ycol], data[, 
            i])
        xnames = c(xnames, nm)
        nr = nr + 1
    }
	nms <- factor(xnames, levels = unique(xnames))
    d = data.frame(d, nms)
    names(d) = c("x", "y", "z", "name")
    d
}
