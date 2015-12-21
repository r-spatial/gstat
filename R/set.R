# $Id: set.q,v 1.7 2008-12-16 14:59:22 edzer Exp $
gstat.set <- function(set) {
	if(!is.list(set))
		stop("set should be a list")
	if (length(set) == 0)
		return(NULL)
	ret = NULL
	n = names(set)
	for (i in (1:length(set))) {
		val = set[[i]]
		if (n[i] == "method")
			str = paste("method: ", val, ";", sep="")
		else {
			if (is.character(val))
				val = paste("'", val, "'", sep = "")
			str = paste("set ", n[i], " = ", val, ";", sep="")
		}
		ret = c(ret, str)
	}
	ret
}

gstat.load.set <- function(set) {
	if(!is.list(set))
		stop("set should be a list")
	if (length(set) == 0)
		return(NULL)
	n = names(set)
	for (i in (1:length(set))) {
		if (n[i] == "method")
			ret = .Call(gstat_set_method, set[[i]])
		else 
			ret = .Call(gstat_set_set, n[i], set[[i]])
	}
	invisible(ret)
}

gstat.load.merge <- function(obj) {

	if (is.character(obj$merge) && length(obj$merge) == 2)
		obj$merge = list(c(obj$merge[1], 1, obj$merge[2], 1))

	if (!is.list(obj$merge))
		stop("merge argument should be list or character vector of lenght 2")

	ret = NULL
	for (i in 1:length(obj$merge)) {
		m = obj$merge[[i]]
		if (is.character(m) && length(m) == 4) {
			id = match(m[c(1,3)], names(obj$data)) - 1 # name ->> id
			if (any(is.na(id)))
				stop(paste("could not match all ids:", m[c(1,3)]))
			col = as.integer(m[c(2,4)]) - 1
			if (any(is.na(col)) || any(col < 0))
				stop("merge: parameters should be positive integers")
			ret = .Call(gstat_set_merge, id[1], col[1], id[2], col[2])
		} else stop(
			"list elements of merge should be lenght 4 character vectors")
	}
	ret
}
