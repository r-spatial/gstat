# $Id: print.variogramModel.q,v 1.5 2009-02-20 13:53:38 edzer Exp $

"print.variogramModel" =
function (x, ...) 
{
    df = data.frame(x)
	shape.models = c("Mat", "Exc", "Cau", "Ste")
	if (!any(match(df[, "model"], shape.models, nomatch=0)))
		df$kappa = NULL
    if (!any(df[, "anis2"] != 1))  {
		df$anis2 = NULL
		df$ang2 = NULL
		df$ang3 = NULL
		if (!any(df[, "anis1"] != 1))  {
			df$anis1 = NULL
			df$ang1 = NULL
		}
	} 
	if (any(match(df[, "model"], "Tab", nomatch=0))) {
		df$maxdist = df$range
		df$range = NULL
    	print(df, ...)
		cat("covariance table:\n")
		tab = attr(x, "table")
		idx = round(seq(1, length(tab), length=6))
		print(tab[idx])
	} else
    	print(df, ...)
	invisible(x)
}
