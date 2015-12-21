/* unit basic variogram models */
double fn_nugget(double h, double *r);
double fn_linear(double h, double *r);
double fn_circular(double h, double *r);
double fn_spherical(double h, double *r);
double fn_bessel(double h, double *r);
double fn_gaussian(double h, double *r);
double fn_exclass(double h, double *r);
double fn_matern(double h, double *r);
double fn_matern2(double h, double *r);
double fn_exponential(double h, double *r);
double fn_pentaspherical(double h, double *r);
double fn_periodic(double h, double *r);
double fn_wave(double h, double *r);
double fn_hole(double h, double *r);
double fn_logarithmic(double h, double *r);
double fn_power(double h, double *r);
double fn_spline(double h, double *r);
double fn_legendre(double h, double *r);
double fn_intercept(double h, double *r);

/* the following functions are not all defined */
double da_is_zero(double h, double *r); /* NUG, INT */
double da_fn_linear(double h, double *r);
double da_fn_circular(double h, double *r);
double da_fn_spherical(double h, double *r);
double da_fn_bessel(double h, double *r);
double da_fn_gaussian(double h, double *r);
double da_fn_exponential(double h, double *r);
double da_fn_pentaspherical(double h, double *r);
double da_fn_periodic(double h, double *r);
double da_fn_wave(double h, double *r);
double da_fn_hole(double h, double *r);
double da_fn_logarithmic(double h, double *r);
double da_fn_power(double h, double *r);

/* unit derivative-to-range of basic variogram models */
double da_fn_exponential(double h, double *r);
double da_fn_nugget(double h, double *r);
