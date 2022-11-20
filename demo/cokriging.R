# $Id: cokriging.R,v 1.4 2006-02-10 19:05:02 edzer Exp $
library(stars)
library(gstat)
data(meuse, package = "sp")
meuse = st_as_sf(meuse, coords = c("x", "y"))
data(meuse.grid, package = "sp")
meuse.grid = st_as_stars(meuse.grid)

# cokriging of the four heavy metal variables
# create gstat object, stepwise:
gstat(id="zn", formula=log(zinc)~1, data=meuse, nmax = 10) |>
  gstat("cu", log(copper)~1, meuse, nmax = 10) |>
  gstat("cd", log(cadmium)~1, meuse, nmax = 10) |>
  gstat("pb", log(lead)~1, meuse, nmax = 10) |>
  gstat(model=vgm(1, "Sph", 900, 1), fill.all=T) -> meuse.g

x <- variogram(meuse.g, cutoff=1000)
meuse.fit = fit.lmc(x, meuse.g)
plot(x, model = meuse.fit)
z <- predict(meuse.fit, newdata = meuse.grid)

z[c(1,3,5,7)] |> merge() |> plot()
# compute & plot standard errors:
z[c(2,4,6,9)] |> setNames(paste0(c("zn", "cu", "pb", "cd"), ": se")) |> 
   merge() |> sqrt() |> plot()

# old-style, with sp:
# indicator cokriging for the 9 percentiles of zinc:
q <- quantile(meuse$zinc, seq(.1,.9,.1))
gstat(id = "zn1", formula = I(zinc < q[1])~1, 
	data = meuse, nmax = 7, beta = .1, set = list(order = 4, zero = 1e-5)) |>
  gstat("zn2", I(zinc < q[2])~1, meuse, nmax = 7, beta=.2) |>
  gstat("zn3", I(zinc < q[3])~1, meuse, nmax = 7, beta=.3) |>
  gstat("zn4", I(zinc < q[4])~1, meuse, nmax = 7, beta=.4) |>
  gstat("zn5", I(zinc < q[5])~1, meuse, nmax = 7, beta=.5) |>
  gstat("zn6", I(zinc < q[6])~1, meuse, nmax = 7, beta=.6) |>
  gstat("zn7", I(zinc < q[7])~1, meuse, nmax = 7, beta=.7) |>
  gstat("zn8", I(zinc < q[8])~1, meuse, nmax = 7, beta=.8) |>
  gstat("zn9", I(zinc < q[9])~1, meuse, nmax = 7, beta=.9) |>
  gstat(model=vgm(1, "Sph", 900, 1), fill.all=T) -> meuse.i

x <- variogram(meuse.i, cutoff=1000)
meuse.fit = fit.lmc(x, meuse.i, correct.diagonal = 1.01)
plot(x, model = meuse.fit)
z <- predict(meuse.fit, newdata = meuse.grid)

z[c(1,3,5,7,9,11,13,15,17)] |> 
	setNames(paste("est.Pr(Zn < ", signif(q,4), ")", sep = "")) |> merge() |> plot()
