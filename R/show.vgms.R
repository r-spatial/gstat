# $Id: show.vgms.q,v 1.6 2008-12-15 14:27:29 edzer Exp $

"show.vgms" <-
function(min = 1e-12 * max, max = 3, n = 50, sill = 1, range = 1,
	models = as.character(vgm()$short[c(1:17)]), nugget = 0, kappa.range = 0.5,
	plot = TRUE, ..., as.groups = FALSE) 
{

	zero.range.models = c("Nug", "Int", "Lin", "Err")
	# print(models)
	L = max(length(sill), length(range), length(nugget), length(models), length(kappa.range))
	sill = rep(sill, length.out = L)
	range = rep(range, length.out = L)
	nugget = rep(nugget, length.out = L)
	i = 0
	if (length(kappa.range) > 1) { # loop over kappa values for Matern model:
		if (missing(models))
			models = "Mat"
		stopifnot(models == "Mat" || models == "Ste" || models == "Exc")
		data = matrix(NA, n * length(kappa.range), 2)
		v.level = rep("", n * length(kappa.range))
		for (kappa in kappa.range) {
			v = vgm(sill[i+1], models, range[i+1], nugget = nugget[i+1], kappa = kappa)
			x = variogramLine(v, 0, 1, 0)
			data[(i*n+1), ] = as.matrix(x)
			x = variogramLine(v, max, n - 1, min)
			data[(i*n+2):((i+1)*n), ] = as.matrix(x)
			m.name = paste("vgm(", sill[i+1], ",\"", models, "\",", range, sep = "")
			if (nugget[i+1] > 0)
				m.name = paste(m.name, ",nugget=", nugget[i+1], sep = "")
			m.name = paste(m.name, ",kappa=", kappa, ")", sep = "")
			v.level[(i*n+1):((i+1)*n)] = rep(m.name, n)
			i =  i + 1
		}
	} else {
		models = rep(models, length.out = L)
		data = matrix(NA, n * length(models), 2)
		v.level = rep("", n * length(models))
		for (m in models) {
			this.range = ifelse(!is.na(pmatch(m, zero.range.models)), 0, range[i+1])
			v = vgm(sill[i+1], m, this.range, nugget = nugget[i+1], kappa = kappa.range)
			x = variogramLine(v, 0, 1, 0)
			data[(i*n+1), ] = as.matrix(x)
			x = variogramLine(v, max, n - 1, min)
			data[(i*n+2):((i+1)*n), ] = as.matrix(x)
			m.name = paste("vgm(", sill[i+1], ",\"", m, "\",", this.range, sep = "")
			if (nugget[i+1] > 0)
				m.name = paste(m.name, ",nugget=", nugget[i+1], sep = "")
			m.name = paste(m.name, ")", sep = "")
			v.level[(i*n+1):((i+1)*n)] = rep(m.name, n)
			i = i + 1
		}
	}
	dframe = data.frame(semivariance = data[,2], distance = data[,1], 
		model = factor(v.level, levels = unique(v.level)))
	vgm.panel = function(x,y, ...) {
		n = length(x)
		lpoints(x[1],y[1])
		llines(x[2:n],y[2:n])
	}
	vgm.panel2 = function(x, y, subscripts, groups, ...) {
		lpoints(0, 0, col = 1)
		panel.superpose(x, y, subscripts, groups, ...)
	}
	if (!plot)
		dframe
	else {
		if (as.groups) {
			model = 0 # avoid NOTE on cran check
			xyplot(semivariance ~ distance, groups = model, dframe[dframe$distance > 0,], 
				panel = vgm.panel2, as.table = TRUE, auto.key = TRUE, type = 'l', ...)
		} else
			xyplot(semivariance ~ distance | model, dframe,
				panel = vgm.panel, as.table = TRUE, ...)
	}
}
