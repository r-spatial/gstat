# $Id: gstat.debug.q,v 1.3 2006-02-10 19:01:07 edzer Exp $

"gstat.debug" <- function(level = 0) {
	invisible(.Call(gstat_debug_level, as.integer(level)))
}
