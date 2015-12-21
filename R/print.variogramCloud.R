# $Id: print.variogramCloud.q,v 1.4 2007-10-18 10:13:13 edzer Exp $

as.data.frame.variogramCloud = function(x, row.names, optional, ...) {
	.BigInt = attr(x, ".BigInt")
	x$left = x$np %% .BigInt + 1
	x$right = x$np %/% .BigInt + 1
	x$np = NULL
	class(x) = "data.frame"
	x
}

print.variogramCloud <- function (x, ...) {
    print(as.data.frame(x), ...)
}
