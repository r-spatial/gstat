# $Id: gstat.formula.q,v 1.8 2007-06-08 06:45:52 edzer Exp $

"gstat.formula" <-
function (formula, data)
{
	# check for duplicated pixels; if yes coerce to SpatialPointsDataFrame:
	if (is(data, "SpatialPixels") && anyDuplicated(data@grid.index) != 0)
		gridded(data) = FALSE

    m = model.frame(terms(formula), as(data, "data.frame"), na.action = na.fail)
    Y = model.extract(m, "response")
    if (length(Y) == 0)
        stop("no response variable present in formula")
    Terms = attr(m, "terms")
    X = model.matrix(Terms, m)
    has.intercept = attr(Terms, "intercept")

	if (gridded(data))
		grid = gridparameters(data)
	else
		grid = numeric(0)

	xlevels = .getXlevels(Terms, m)

    list(y = Y, locations = coordinates(data), X = X, call = call,
        has.intercept = has.intercept, grid = as.double(unlist(grid)),
		xlevels = xlevels)
}
