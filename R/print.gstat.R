# $Id: print.gstat.q,v 1.7 2006-02-10 19:01:07 edzer Exp $

"print.gstat" <-
function (x, ...) 
{
    if (missing(x) || !inherits(x, "gstat"))
        stop("wrong call")
    data.names <- names(x$data)
    if (length(data.names)) 
        cat("data:\n")
    for (n in data.names) {
        fstr = paste(x$data[[n]]$formula[c(2, 1, 3)], collapse = "")
        #lstr = paste(x$data[[n]]$locations[c(1, 2)], collapse = "")
        cat(n, ": formula =", fstr, ";")
        if (!is.null(x$data[[n]]$data)) {
            data.dim = dim(x$data[[n]]$data)
            cat(" data dim =", data.dim[1], "x", data.dim[2])
        }
        else {
            if (x$data[[n]]$dummy) 
                cat(" dummy data")
            else cat(" NULL data")
        }
		if (x$data[[n]]$nmax != Inf)
			cat(" nmax =", x$data[[n]]$nmax)
		if (x$data[[n]]$nmin > 0)
			cat(" nmin =", x$data[[n]]$nmin)
		if (x$data[[n]]$maxdist < Inf)
			cat(" radius =", x$data[[n]]$maxdist)
		if (x$data[[n]]$vfn > 1)
			cat(" variance function =", 
				c("identity", "mu", "mu(1-mu)", "mu^2", "mu^3")[x$data[[n]]$vfn])
		if (length(x$data[[n]]$beta) > 0)
			cat(" beta =", x$data[[n]]$beta)
		if (x$data[[n]]$degree > 0)
			cat(" degree =", x$data[[n]]$degree)
        cat("\n")
    }
    xx.names = xx = NULL
    for (n in data.names) {
        m = x$model[[n]]
        if (!is.null(m)) {
            xx = rbind(xx, m)
            if (nrow(m) == 1) 
                xx.names = c(xx.names, n)
            else xx.names = c(xx.names, paste(n, "[", 1:nrow(m), 
                "]", sep = ""))
        }
    }
    if (length(data.names) > 1) {
        for (j in 2:length(data.names)) {
            for (i in 1:(j - 1)) {
                n = cross.name(data.names[i], data.names[j])
                m = x$model[[n]]
                if (!is.null(m)) {
                  xx = rbind(xx, m)
                  if (nrow(m) == 1) 
                    xx.names = c(xx.names, n)
                  else xx.names = c(xx.names, paste(n, "[", 1:nrow(m), 
                    "]", sep = ""))
                }
            }
        }
    }
    if (!is.null(xx)) {
        cat("variograms:\n")
        row.names(xx) = xx.names
        print(xx, ...)
    }
    if (!is.null(x$set)) {
        s = gstat.set(x$set)
        for (i in 1:length(s)) cat(s[i], "\n")
    }
    if (!is.null(x$locations))
		print(x$locations)
    invisible(x)
}
