# $Id: xyz2img.q,v 1.4 2006-02-10 19:01:07 edzer Exp $

"xyz2img" <-
function (xyz, zcol = 3, xcol = 1, ycol = 2, tolerance = 10 * .Machine$double.eps) 
{
    if (ncol(xyz) < 3) 
        stop("xyz object should have at least three columns")
    z = xyz[, zcol]
    x = xyz[, xcol]
    y = xyz[, ycol]
    xx = sort(unique(x))
    yy = sort(unique(y))
    nx = length(xx)
    ny = length(yy)
    nmax = max(nx, ny)
    difx = diff(xx)
    if (diff(range(unique(difx))) > tolerance) 
        stop("x intervals are not constant")
    dify = diff(yy)
    if (diff(range(unique(dify))) > tolerance) 
        stop("y intervals are not constant")
    dx = mean(difx)
    dy = mean(dify)
    xmin = min(xx)
    xmax = max(xx)
    xrange = xmax - xmin
    ymin = min(yy)
    ymax = max(yy)
    yrange = ymax - ymin
    row = round((x - xmin)/dx) + 1
    col = round((y - ymin)/dy) + 1
	zz = rep(as.numeric(NA), nx * ny)
	zz[row + nx * (col - 1)] = z
	zz = matrix(zz, nrow = nx, ncol = ny)
    list(x = seq(xmin, xmax, dx), y = seq(ymin, ymax, dy), z = zz)
}
