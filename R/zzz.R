# $Id: zzz.q,v 1.10 2006-02-10 19:01:07 edzer Exp $

### NAMESPACE VERSION:

.onLoad <- function(lib, pkg) {
	# remove the require() call for 2.0.0:
	# require(lattice)
	.Call(gstat_init, as.integer(1))
}

### pre-NAMESPACE VERSION:
## ".First.lib" <-
## function(lib, pkg) {
## 	require(lattice)
## 	library.dynam("gstat", pkg, lib)
## 	.Call(gstat_init, as.integer(1))
## }
 
variogram <- function(object, ...) UseMethod("variogram")
