# $Id: gstat.cv.q,v 1.9 2009-10-30 16:11:21 edzer Exp $

"gstat.cv" <-
function (object, nfold = nrow(object$data[[1]]$data), remove.all = FALSE, 
	verbose = interactive(), all.residuals = FALSE, ...) 
{
	if (!inherits(object, "gstat")) 
		stop("first argument should be of class gstat")
	var1 = object$data[[1]]
	data = var1$data
	formula = var1$formula
	if (all.residuals) {
		nc = length(object$data)
		ret = data.frame(matrix(NA, nrow(data), nc))
	} else {
		cc = coordinates(data)
		rownames(cc) = NULL
		df = data.frame(matrix(as.numeric(NA), nrow(data), 2))
		ret = SpatialPointsDataFrame(cc, df)
	}
	if (missing(nfold)) 
		nfold = 1:nrow(data)
	if (length(nfold) == nrow(data))
		fold = nfold
	else if (nfold < nrow(data)) 
		fold = sample(nfold, nrow(data), replace = TRUE)
	else fold = 1:nrow(data)

	if (all.residuals || (remove.all && length(object$data) > 1)) {
		all.data = list()
		for (v in 1:length(object$data))
			all.data[[v]] = object$data[[v]]$data
	}

	if (verbose)
		pb <- txtProgressBar(1,length(unique(fold)),style=3)

	for (i in sort(unique(fold))) {
		if (verbose)
			setTxtProgressBar(pb, i) 
		sel = which(fold == i)
		object$data[[1]]$data = data[-sel, ]
		if (remove.all && length(object$data) > 1) {
			for (v in 2:length(object$data)) {
				varv = object$data[[v]]
				varv$data = all.data[[v]]
				#atv = gstat.formula(varv$formula, varv$data)$locations
				#at1 = gstat.formula(formula, data[sel, ])$locations
				atv = coordinates(varv$data)
				at1 = coordinates(data[sel,])
				cc = rbind(atv, at1)
				rownames(cc) = NULL # as there will be duplicates
				all = SpatialPoints(cc)
				zd = zerodist(all)
				skip = zd[, 1]
				object$data[[v]]$data = varv$data[-skip, ]
			}
		}
		x = predict(object, newdata = data[sel, ], ...)
		if (all.residuals) {
			for (i in 1:length(object$data)) {
				var.i = object$data[[i]]
				data.i = all.data[[i]]
				formula.i = var.i$formula
				observed = gstat.formula(formula.i, data.i)$y[sel]
				pred.name = paste(names(object$data)[i], "pred", sep = ".")
				residual = as.numeric(observed - x[[pred.name]])
				ret[sel, i] = residual
			}
		} else {
			ret[[1]][sel] = x[[1]]
			ret[[2]][sel] = x[[2]]
		}
	}
	if (verbose)
		cat("\n")

	if (! all.residuals) {
		names(ret) = names(x)[1:2]
		ret$observed = gstat.formula(formula, data)$y
		pred.name = paste(names(object$data)[1], "pred", sep = ".")
		ret$residual = ret$observed - ret[[pred.name]]
		var.name = paste(names(object$data)[1], "var", sep = ".")
		ret$zscore = ret$residual/sqrt(ret[[var.name]])
		ret$fold = fold
	} else
		names(ret) = names(object$data)

	if (!is.null(object$locations))
		ret = as.data.frame(ret)

	ret
}
