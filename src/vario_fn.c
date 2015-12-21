/*
 * vario_fn.c: contains elementary variogram model functions
 */
#include <stdio.h> /* req'd by userio.h */
#include <float.h>
#include <math.h>

#include "defs.h"
/* #define MATHLIB_STANDALONE */
#include <Rmath.h>

#include "userio.h"
#include "utils.h"
#include "vario_fn.h"

/*
 * Copyright (C) 1994, Edzer J. Pebesma
 * 
 * basic variogram functions
 */

#define MIN_BESS 1.0e-3
#ifndef PI
#   define PI 3.14159265359
#endif


double fn_nugget(double h, double *r) {
	return (h == 0.0 ? 0.0 : 1.0);
}

double fn_linear(double h, double *r) {
	if (h == 0)
		return 0.0;
	if (*r == 0)
		return h; /* 1lin() or 1 lin(0): slope 1 (and no range) */
	return (h >= *r ? 1.0 : h/(*r));
}

double da_fn_linear(double h, double *r) {
	if (*r == 0)
		return 0.0; /* 1lin() or 1 lin(0): slope 1 (and no range) */
	if (h > *r)
		return 0.0;
	return -h/((*r) * (*r));
}

double fn_circular(double h, double *r) {
	double hr;
	if (h == 0.0)
		return 0.0;
	if (h >= *r)
		return 1.0;
	hr = h/(*r);
	/*
	 * return 1.0 + (2.0/PI) * (hr * sqrt(1.0 - hr * hr) - acos(hr));
	 * probably equivalent to:
	 */
	return (2.0/PI) * (hr * sqrt(1.0 - hr * hr) + asin(hr));
}

double fn_spherical(double h, double *r) {
	double hr;
	if (h == 0)
		return 0.0;
	if (h >= *r)
		return 1.0;
	hr = h/(*r);
	return hr * (1.5 - 0.5 * hr * hr);
}

double da_fn_spherical(double h, double *r) {
	double hr2;
	if (h > *r)
		return 0.0;
	hr2 = h / ((*r) * (*r));
	return 1.5 * hr2 * (-1.0 + h * hr2);
}

double fn_bessel(double h, double *r) {
	double hr;
	hr = h/(*r);
	if (hr < MIN_BESS) 
		return 0.0;
	/* return 1.0 - hr * bessk1(hr); */
	return 1.0 - hr * bessel_k(hr, 1.0, 1.0);
}

double fn_gaussian(double h, double *r) {
	double hr;
	if (h == 0.0)
		return 0.0;
	hr = h/(*r);
	return 1.0 - exp(-(hr*hr)); 
}

double da_fn_gaussian(double h, double *r) {
	double hr;
	hr = h / (*r);
	return (-hr /(*r)) * exp(-(hr * hr));
}

double fn_exponential(double h, double *r) {
	if (h == 0.0)
		return 0.0;
	return 1.0 - exp(-h/(*r)); 
}

double fn_exclass(double h, double *r) {
	if (h == 0.0)
		return 0.0;
	return 1.0 - exp(-pow(h/r[0],r[1])); 
}

double da_fn_exponential(double h, double *r) {
	double hr;
	hr = -h/(*r);
	return (hr / (*r))  * exp(hr); 
}

double fn_pentaspherical(double h, double *r) {
	double hr, h2r2;
	if (h == 0.0)
		return 0.0;
	if (h >= *r)
		return 1.0;
	hr = h/(*r);
	h2r2 = hr * hr;
	return hr * ((15.0/8.0) + h2r2 * ((-5.0/4.0) + h2r2 * (3.0/8.0)));
}

double da_fn_pentaspherical(double h, double *r) {
	double hr2;
	if (h >= *r)
		return 0.0;
	hr2 = h / ((*r) * (*r));
	return hr2*((-15.0/8.0) + hr2*((15.0/4.0)*h - (15.0/8.0)*h*h*hr2));
}

double fn_periodic(double h, double *r) {
	if (h == 0.0)
		return 0.0;
	return 1.0 - cos(2.0 * PI * h/(*r));
}

double da_fn_periodic(double h, double *r) {
	return (2.0 * PI * h/((*r) * (*r))) * sin(2.0 * PI * h/(*r));
}

double fn_wave(double h, double *r) {
	if (h == 0.0)
		return 0.0;
	return 1.0 - (*r) * sin(PI * h/(*r)) / (PI * h);
}

double da_fn_wave(double h, double *r) {
	return cos(PI * h / (*r)) / (*r) - sin(PI * h / (*r)) / PI * h;
}

double fn_hole(double h, double *r) {
	if (h == 0.0)
		return 0.0;
	return 1.0 - sin(h/(*r))/(h/(*r));
}

double da_fn_hole(double h, double *r) {
	double hr, hr2;
	hr = h/(*r);
	hr2 = h/((*r)*(*r));
	return hr2 * sin(hr) + hr * hr2 * cos(hr);
}

double fn_logarithmic(double h, double *r) {
	if (h == 0.0)
		return 0.0;
	return log(h + *r); 
}

double da_fn_logarithmic(double h, double *r) {
	return 1/(*r);
}

double fn_power(double h, double *r) {
	if (h == 0.0)
		return 0.0;
	return pow(h, *r); 
}

double da_fn_power(double h, double *r) {
	return log(h) * pow(h, *r);
}

double fn_spline(double h, double *r) {
	if (h == 0.0)
		return 0.0;
	return h * h * log(h); 
}

double fn_legendre(double h, double *r) {
	/* *r is range; h is angular lag */
	double r2, ang;
	if (h == 0.0)
		return 0.0;
	r2 = r[0] * r[0];
	ang = h / (PI * 6378.137);
	/* printf("dist: %g, ang: %g\n", h, ang); */
	return 2.0 - (1.0 - r2)/(1 - 2.0 * r[0] * cos(ang) + r2);
}

double fn_intercept(double h, double *r) {
	return 1.0;
}

double da_is_zero(double h, double *r) {
	return 0.0;
}

double fn_matern(double h, double *p) {
	double hr, ans, phi, kappa;

	phi = p[0];
	kappa = p[1];
    if (h == 0.0)
		return(0.0);
	if (h > 600 * phi)
		return(1.0);
	hr = h/phi;
	ans = (pow(2.0, -(kappa - 1.0))/gammafn(kappa)) *
			pow(hr, kappa) * bessel_k(hr, kappa, 1.0);
	/* ans was for correlation; */
    return 1.0 - ans;
}

/* 
 * Ste: M. Stein's representation of the Matern model
 *
 * According to Hannes Kazianka, in R this would be
h=distance matrix
delta=c(RANGE,KAPPA)
maternmodel<-function(h,delta){
	matern<-besselK(2*delta[2]^(1/2)*h/delta[1],delta[2])
	ifelse(matern==0,0,ifelse(!is.finite(matern),1,
	1/(2^(delta[2] -
	1)*gamma(delta[2]))*(2*delta[2]^(1/2)*h/delta[1])^delta[2]*matern))
}
delta<-c(RANGE,KAPPA)
h=Distance Matrix

maternmodel<-function(h,delta){
	matern<-besselK(2*delta[2]^(1/2)*h/delta[1],delta[2])
	multipl<-1/(2^(delta[2] - 1)*gamma(delta[2]))*(2*delta[2]^(1/2)*h/delta[1])^delta[2]
	ifelse(matern==0 | !is.finite(multipl),0,ifelse(!is.finite(matern),1,
	multipl*matern))
}

Now 0*Inf is impossible.

Hannes
*/
double fn_matern2(double h, double *p) {

	double *delta, bes, mult;
	if (h == 0.0)
		return 0.0;
	delta = p;
	h = h / delta[0];
	bes = bessel_k(2.0 * sqrt(delta[1]) * h, delta[1], 1.0);
	if (!isfinite(bes))
		return 0.0;
	if (bes == 0.0)
		return 1.0;
	mult = pow(2.0, 1.0 - delta[1]) / gammafn(delta[1]) * 
		pow((2.0 * sqrt(delta[1]) * h), delta[1]);
	if (!isfinite(mult))
		return 1.0;
	return 1.0 - bes * mult;
}
