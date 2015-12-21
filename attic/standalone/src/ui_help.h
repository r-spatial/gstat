char *helpmenu[] = {
"   Command              Action",
"",
"    j,<Down>       move down",
"    k,<Up>         move up",
"    l,<Enter>      choose item or get next option",
"    h,<Left>       choose item or get previous option",
"    <Tab>          show plot immediately",
"    c              write gstat commands to file (data, variograms, ...)",
"    e              write sample variogram estimates to file",
"    g              write gnuplot commands to file",
"    p              write plot to PNG file",
"    P              write plot to Postscript (EPS) file",
"    G              write plot to gif file",
"    f              fix, unfix ranges in fit",
"    n              do, do not display of numbers in plot",
"    z              change estimates at h=0 (duplicates)",
"    s              show file (default: sample variogram estimates file)",
"    +,-            next/previous variable",
"    r,^L           redraw screen",
"    ?              print help on current field",
"    q,x            quit, exit (q prompts for saving)",
"    !              shell command",
NULL };

char *ui0_help0[] = {
"Data is read from ASCII files or grid maps. An ASCII tabel has variables",
"on columns and records on rows, and may start with a simplifiedGeoEAS",
"header (line 1: header; line 2: n columns; line 3..n+2 column names)).",
"Rows are separated by newlines, columns by white space (spaces, tabs).",
NULL };

char *ui0_help1[] = {
"Choose variable for which semivariogram or covariogram will be",
"calculated.  Choose pair of variables for which a cross semivariogram or",
"cross covariogram will be calculated",
NULL };

char *ui0_help2[] = {
"If one variable is read: choose from semivariogram or covariogram. If more",
"variables are read: choose from semivariogram, covariogram, (pseudo) cross",
"variogram or cross covariogram. Let z or (z1,z2) be the selected (pair of)",
"variable(s) with mean values m(x) or (m1(x),m2(x)).",
"",
"   Function                         Estimate of",
"",
"semivariogram(h)                  E((Z(x)-m(x)-Z(x+h)-m(x+h))^2)",
"covariogram(h)                    E((Z(x)-m(x))(Z(x+h)-m(x+h)))",
"Pseudo Cross semivariogram(h)     E(((Z1(x)-m1(x))-(Z2(x+h)-m2(x+h)))^2)",
"Cross semivariogram(h)            E((Z1(x)-m1(x)-Z1(x+h)-m1(x+h))",
"                                        (Z2(x)-m2(x)-Z2(x+h)-m2(x+h)))",
"Cross covariogram(h)              E((Z1(x)-m1(x))(Z2(x+h)-m2(x+h)))",
"",
"The choice between cross semivariogram and pseudo cross semivariogram is",
"left to the program. The latter will be chosen when Z1 and Z2 do not share",
"the same locations (and order).",
"",
"m(x) or (m1(x),m2(x)) are ordinary least squares estimates by default",
NULL };

char *ui0_help3[] = {
"Cutoff is the maximum distance: pairs of points with a separation distance",
"less than the cutoff will be included in the variogram estimates.",
"",
"Width denotes the width of the intervals for which pairs of points are merged",
"into subsequent estimate. An interval width of 0 will produce the variogram ",
"cloud",
"",
"Default values are one third of the maximum distance (span), subdivided",
"into 15 equally sized intervals (controllable by `set fraction' and `set",
"intervals')",
NULL };

char *ui0_help4[] = {
"In 2 dimensions, direction is defined (in degrees) as the clockwise angle",
"with the positive y (North):",
"    0     45 ",
"    |   /    ",
"    | /      ",
"----+---- 90 ",
"    | \\     ",
"    |   \\   ",
"   180   135 ",
"",
"In 3 dimensions, the vertical direction is in degrees, positive upwards from",
"the direction in the <x,y> plane.",
"",
"The tolerance (+/-) is the angle deviation from the direction defined, ",
"tolerated for a pair of points to be included in the semivariance ",
"(or covariance) estimate",
NULL };

char *ui0_help5[] = {
"A variogram model consists of (the sum of) one (or more) basic models",
"e.g.,  `1 Nug(0)'  or  `1 Nug(0) + 1 Sph(35)' (without quotes)",
"the first number specifies the (partial) sill, the second (between the",
"parens) the range. Three examples of commonly used basic models:",
"`1.5 Nug(0)' -- nugget model with partial sill (nugget value) of 1.5",
"`1.0 Sph(4)' -- spherical model with partial sill 1 and range 4",
"`2.5 Exp(1)' -- exponential model with partial sill 2.5, range parameter 1",
"                and effective range (where it reaches 95% of its height) 3",
"",
"If this field allready contains a variogram model, then pressing <right>",
"or <enter> allows the modification of individual parameters or model type,",
"whereas pressing <left> allows entering of a completely new variogram model",
NULL };

char *ui0_help6[] = {
"Choose a fit method for fitting the model to the sample variogram,",
"sample covariogram or to the data (reml)",
"",
"wls: iteratively reweighted least squares estimate",
"weights n(h): use the number of point pairs as weights",
"weights n(h)/mod(h)^2: use weights n(h)/(model(h)*model(h))",
"[gnuplot] : let gnuplot do the fit",
"reml: restricted maximum likelihood estimates of sill parameters",
" (equivalent to iterated minque or mivque)",
"",
"Notes",
"1. for non-linear fitting (a fit that includes the ranges)",
"   the Marquardt-Levenberg algorithm is used",
"2. gnuplot writes a report on the fit to the file `fit.log'",
"3. reml takes a long time when applied to more than e.g. 100 data points",
NULL };

char *ui0_help7[] = {
"Calculate sample variogram, fit model if requested and show plot",
"(in order to show something, gnuplot must be present).",
NULL };
