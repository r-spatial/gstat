# $Id: image.data.frame.q,v 1.4 2006-02-10 19:01:07 edzer Exp $

"image.data.frame" <-
function (x, zcol = 3, xcol = 1, ycol = 2, asp = 1, ...)
{
    image.default(xyz2img(xyz = x, zcol = zcol, xcol = xcol, ycol = ycol),
		asp = asp,
		...)
}
