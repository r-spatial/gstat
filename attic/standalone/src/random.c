/*
    Gstat, a program for geostatistical modelling, prediction and simulation
    Copyright 1992, 2011 (C) Edzer Pebesma

    Edzer Pebesma, edzer.pebesma@uni-muenster.de
	Institute for Geoinformatics (ifgi), University of Münster 
	Weseler Straße 253, 48151 Münster, Germany. Phone: +49 251 
	8333081, Fax: +49 251 8339763  http://ifgi.uni-muenster.de 

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version. As a special exception, linking 
    this program with the Qt library is permitted.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

    (read also the files COPYING and Copyright)
*/

/*
 * random.c: distribution functions and random number generators
 * uses: error.c userio.h random.h
 * Tue Jun  1 11:25:26 METDST 1993
 * adapted Fri Mar 10 10:47:24 WET 1995
 * (c) E.J. Pebesma
 * time_seed() is called when seed is equal to 0;
 * 0: Marsaglia's congruential random number generator (rn.[ch])
 * 1: drand48() (in stdlib on unix systems)
 */

#include <stdio.h>
#include <stdlib.h> /* rand(), drand48() */
#include <math.h> 	/* sqrt() */

#include "defs.h" /* may define HAVE_DRAND48 */
#include "utils.h"

#ifdef HAVE_LIBGSL
# include <gsl/gsl_rng.h>
# include <gsl/gsl_randist.h>
static gsl_rng *rng = NULL;
#endif /* HAVE_LIBGSL */

#ifdef TIME_WITH_SYS_TIME
# include <sys/time.h> /* gettimeofday */
# include <time.h>     /* time() */
#else
# if HAVE_SYS_TIME
#  include <sys/time.h>
# else
#  include <time.h>
# endif
#endif

#include "debug.h" 

#ifdef PCRCALC /* use functions for a dynamic lybrary in pcrcalc */
# define ER_IMPOSVAL 0
# define ErrMsg(a,b) {printf("Error: %s\n", b); exit(a);}
# define printlog printf
 int debug_level = 0;
 int gl_secure = 0;
 char *gl_mv_string = "NA";
#else
# include "userio.h"
#endif

#include "random.h"

static char start_up[100];

static unsigned long int seed = 0;
static int init = 0;
static unsigned long int time_seed(void);
static void init_random(void);
static void print_start_up(void);

static void start_random_number(int a, int b);
double get_next_random_number(void);
double my_normal(void);

static struct {
	double (*r_unif)(void);
	double (*r_normal)(void);
} my_rng = { NULL, NULL };

double my_gsl_uniform(void);
double my_gsl_normal(void);

#ifdef HAVE_LIBGSL
double my_gsl_uniform(void){
	return(gsl_rng_uniform(rng));
}
double my_gsl_normal(void){
	return gsl_ran_gaussian(rng, 1.0);
}
#endif

void set_rng_functions(
		double (*unif)(void), 
		double (*norm)(void), 
		const char *name) {
	my_rng.r_unif = unif;
	my_rng.r_normal = norm;
	sprintf(start_up, "%s", name);
	init = 1; /* the caller's responsibility, obviously */
	return;
}

unsigned long int get_seed(void) {
	return seed;
}

void set_seed(unsigned long int i) {
	if (i == 0)
		seed = time_seed();
	else
		seed = i;
	init_random();
	return;
}

static unsigned long int time_seed(void) {

#ifdef HAVE_GETTIMEOFDAY
	struct timeval tv;

	if (gettimeofday(&tv, NULL) == 0)
		return (unsigned long int) tv.tv_sec * 1000000 + tv.tv_usec;
	else
		return (unsigned long int) time(NULL);
#else
	return (unsigned long int) time(NULL);
#endif
}

static void print_start_up(void) {
	static int done = 0;

	if (done == 0 && debug_level > DB_NORMAL)
		printlog("%s", start_up);
	done = 1;
	return;
}

static void init_random(void) {
	unsigned int a, b;

	init = 1;
#ifdef HAVE_LIBGSL
	gsl_rng_env_setup();
	rng = gsl_rng_alloc(gsl_rng_default);
	if (getenv("GSL_RNG_SEED") == NULL && seed != 0)
		gsl_rng_set(rng, seed);
	sprintf(start_up, "GSL generator type: %s, seed = %lu\n", 
			gsl_rng_name(rng),
			getenv("GSL_RNG_SEED") == NULL && seed != 0 ?
			seed : gsl_rng_default_seed);
	my_rng.r_unif = my_gsl_uniform;
	my_rng.r_normal = my_gsl_normal;
	printlog("using the GSL random number generator\n");
	return;
#endif

	if (seed > 65536) {
		a = seed << 16;
		a = a >> 16;
		b = seed >> 16;
	} else
		b = a = seed;
	start_random_number( (int) a , (int) b);
	sprintf(start_up, "using Marsaglia's random number generator, seed %lu, a %d, b %d\n", seed, a, b);
	my_rng.r_unif = get_next_random_number;
	my_rng.r_normal = my_normal;
	printlog("using Marsaglia's random number generator\n");
	return; 

#ifdef USE_DRAND48 /* effectively outcommented */
	sprintf(start_up, 
		"using drand48() as random number generator, seed %lu\n", seed);
	srand48(seed);
	printlog("using drand48() as random number generator\n");
	my_rng.r_unif = drand48;
	my_rng.r_normal = my_normal;
	return;
#endif

}

/*
 * [pqr]_uniform: functions for uniform distribution over [0,1]
 */
double p_uniform(double z) {
/*
 * returns the probability of getting a uniform [0,1] distr. value, <= z
 */
	if (z < 0.0)
		return 0.0;
	if (z > 1.0)
		return 1.0;
	return z;
}

double q_uniform(double p) {
/*
 * returns the p-quantile of a uniform distribution
 */
	assert(p >= 0.0 && p <= 1.0);
	return p;
}

double r_uniform(void) {

	assert(my_rng.r_unif != NULL);
	assert(init != 0);
	print_start_up();
	return my_rng.r_unif();
}

#ifndef USING_R
int e_random(int argc, char *argv[]) {

	int un, n, i;

	if (argc != 3) {
		printf("%s [U|N] <n>\n", argv[0]);
		printf("U: uniform, N: standard normal\n");
		printf("<n>: number of values wanted\n");
		exit(0);
	}
	set_seed(0);

	un = *argv[1] == 'U';
	n = atoi(argv[2]);

	set_seed(0);

	if (un)
		for (i = 0; i < n; i++)
			printf("%g\n", r_uniform());
	else
		for (i = 0; i < n; i++)
			printf("%g\n", r_normal());
	return 0;
}
#endif

/*
 * [pqr]_normal: functions for Normal(mean=0,var=1) distribution
 */
double p_normal(double z) {

#define	Z_EPSILON      0.000001       /* accuracy of q_normal approximation */
#define	Z_MAX          6.0            /* maximum meaningful z value */

/*
 * NOTE: both p_normal() and q_normal() are from Gary Perlman's z.c
 * (not copyrighted)
 *
 * FUNCTION p_normal: probability of normal z value
 * ALGORITHM:
 *  Adapted from a polynomial approximation in:
 *   Ibbetson D, Algorithm 209
 *   Collected Algorithms of the CACM 1963 p. 616
 *  Note:
 *   This routine has six digit accuracy, so it is only useful for absolute
 *   z values < 6.  For z values >= to 6.0, p_normal() returns 0.0.
 */

	double	y, x, w;
	
	if (z == 0.0)
		x = 0.0;
	else {
		y = 0.5 * fabs (z);
		if (y >= (Z_MAX * 0.5))
			x = 1.0;
		else if (y < 1.0) {
			w = y*y;
			x = ((((((((0.000124818987 * w
				-0.001075204047) * w +0.005198775019) * w
				-0.019198292004) * w +0.059054035642) * w
				-0.151968751364) * w +0.319152932694) * w
				-0.531923007300) * w +0.797884560593) * y * 2.0;
		} else {
			y -= 2.0;
			x = (((((((((((((-0.000045255659 * y
				+0.000152529290) * y -0.000019538132) * y
				-0.000676904986) * y +0.001390604284) * y
				-0.000794620820) * y -0.002034254874) * y
				+0.006549791214) * y -0.010557625006) * y
				+0.011630447319) * y -0.009279453341) * y
				+0.005353579108) * y -0.002141268741) * y
				+0.000535310849) * y +0.999936657524;
		}
	}
	return (z > 0.0 ? ((x + 1.0) * 0.5) : ((1.0 - x) * 0.5));
}

double q_normal(double p) {
/*
 * FUNCTION q_normal: compute critical z value to produce given probability
 * (the p-quantile of the normal distribution)
 * ALGORITHM
 *  Begin with upper and lower limits for z values (maxz and minz)
 *  set to extremes.  Choose a z value (zval) between the extremes.
 *  Compute the probability of the z value.  Set minz or maxz, based
 *  on whether the probability is less than or greater than the
 *  desired p.  Continue adjusting the extremes until they are
 *  within Z_EPSILON of each other.
 */
	double	minz = -Z_MAX;    /* minimum of range of z */
	double	maxz = Z_MAX;     /* maximum of range of z */
	double	zval = 0.0;       /* computed/returned z value */
	double	pval;     /* prob (z) function, pval := p_normal(zval) */
	
	assert(p > 0.0 && p < 1.0);
	while (maxz - minz > Z_EPSILON) {
		pval = p_normal(zval);
		if (pval > p)
			maxz = zval;
		else
			minz = zval;
		zval = (maxz + minz) * 0.5;
	}
	return zval;
}

double r_normal(void) {
	assert(my_rng.r_normal != NULL);
	print_start_up();
	return(my_rng.r_normal());
}

double my_normal(void) {

	static int iset = 0;
	static double gset;
	double fac, r, v1, v2;
	if  (! iset) {
		do {
			v1 = 2.0 * r_uniform() - 1.0;
			v2 = 2.0 * r_uniform() - 1.0;
			r = v1 * v1 + v2 * v2;
		} while (r >= 1.0 || r == 0.0);
		fac = sqrt(-2.0 * log(r) / r);
		gset = v1 * fac;
		iset = 1;
		return(v2 * fac);
	} else {
		iset = 0;
		return(gset);
	}
}

/* 
 * [pqr]_triangular: functions for symmetric triangular distribution
 * (min=0,max=1)
 */
double p_triangular(double z) {
/*
 * returns the probability of getting a triangular [0,1] distr. value, <= z
 */
 	if (z < 0.5)
 		return 2.0 * z * z;
 	z = 1.0 - z;
 	return 1.0 - (2.0 * z * z);
}

double q_triangular(double p) {
/*
 * returns the p-quantile of a triangular distribution
 */
 	assert(p >= 0.0 && p <= 1.0);
	if (p < 0.5)
		return sqrt(0.5 * p);
	else
		return 1.0 - sqrt(0.5 * (1.0 - p));
}

double r_triangular(void) {
/*
 * return a random number, following a triangular distribution in [0,1].
 * To modify it to a triangularly distributed random number in [min,max]
 * call: t = min + (max - min) * r_triangular();
 * This triangular distribution has max. prob (2) at 0.5, prob 0 at 0 and 1
 */
 	double t;
 	int side;

 	do {
 		t = r_uniform();
 	} while (t <= 0.0 || t > 1.0);
 	/*
 	 * create a half-triangular distributed number on [0,0.5],
 	 * with max. probability (2) on the 0.5 boundary
 	 */
 	t = 0.5 * sqrt(t);
 	/* 
 	 * now decide on which side of the top we're gona put this one
 	 */
 	side = r_uniform() < 0.5;
 	return (side == 1 ? t : 1.0 - t);
}

/*
 * EJP: the following was found on
 * bugs.nosc.mil:/pub/ada/random/random_number.[ch]
 */

/*
 *
 *  Title: 	random_number
 *  Last Mod: 	Fri Mar 18 07:28:57 1988
 *  Author: 	Vincent Broman
 *		<broman@schroeder.nosc.mil>
 */  

/*  
 *  This package makes available Marsaglia's highly portable generator 
 *  of uniformly distributed pseudo-random numbers.
 *  
 *  The sequence of 24 bit pseudo-random numbers produced has a period 
 *  of about 2**144, and has passed stringent statistical tests 
 *  for randomness and independence.
 *  
 *  Supplying two seeds to start_random_number is required once
 *  at program startup before requesting any random numbers, like this:
 *      start_random_number(101, 202);
 *      r := get_next_random_number();
 *  The correspondence between pairs of seeds and generated sequences 
 *  of pseudo-random numbers is many-to-one.
 *  
 *  This package should compile and run identically on any 
 *  machine/compiler which supports >=16 bit integer arithmetic
 *  and >=24 bit floating point arithmetic.
 *  
 *  References:
 *      M G Harmon & T P Baker, ``An Ada Implementation of Marsaglia's
 *      "Universal" random_number Number Generator'', Ada Letters, late 1987.
 *      
 *      G Marsaglia, ``Toward a universal random number generator'',
 *      to appear in the Journal of the American Statistical Association.
 *  
 *  George Marsaglia is at the Supercomputer Computations Research Institute
 *  at Florida State University.
 */  

/*
 *  Title: 	random_number
 *  Last Mod: 	Fri Mar 18 08:52:13 1988
 *  Author: 	Vincent Broman
 *		<broman@schroeder.nosc.mil>
 */

#define P 179
#define PM1 (P - 1)
#define Q (P - 10)
#define STATE_SIZE 97
#define MANTISSA_SIZE 24
#define RANDOM_REALS 16777216.0
#define INIT_C    362436.0
#define INIT_CD  7654321.0
#define INIT_CM 16777213.0

static unsigned int ni;
static unsigned int nj;
static double u[STATE_SIZE];
static double c, cd, cm;


static unsigned int collapse (int anyint, unsigned int size)
/*
 * return a value between 0 and size-1 inclusive.
 * this value will be anyint itself if possible, 
 * otherwise another value in the required interval.
 */
{
	if (anyint < 0)
		anyint = - (anyint / 2);
	while (anyint >= size)
		anyint /= 2;
	return (anyint);
}

static void start_random_number (int seed_a, int seed_b)
/*
 * This procedure initialises the state table u for a lagged 
 * Fibonacci sequence generator, filling it with random bits 
 * from a small multiplicative congruential sequence.
 * The auxilliaries c, ni, and nj are also initialized.
 * The seeds are transformed into an initial state in such a way that
 * identical results are guaranteed across a wide variety of machines.
 */
{
double s, bit;
unsigned int ii, jj, kk, mm;
unsigned int ll;
unsigned int sd;
unsigned int elt, bit_number;

	sd = collapse (seed_a, PM1 * PM1);
	ii = 1 + sd / PM1;
	jj = 1 + sd % PM1;
	sd = collapse (seed_b, PM1 * Q);
	kk = 1 + sd / PM1;
	ll = sd % Q;
	if (ii == 1 && jj == 1 && kk == 1)
		ii = 2;

	ni = STATE_SIZE - 1;
	nj = STATE_SIZE / 3;
	c  = INIT_C;
	c /= RANDOM_REALS;		/* compiler might mung the division itself */
	cd = INIT_CD;
	cd /= RANDOM_REALS;
	cm = INIT_CM;
	cm /= RANDOM_REALS;

	for (elt = 0; elt < STATE_SIZE; elt += 1) {
		s = 0.0;
		bit = 1.0 / RANDOM_REALS;
		for (bit_number = 0; bit_number < MANTISSA_SIZE; bit_number += 1) {
	    	mm = (((ii * jj) % P) * kk) % P;
	    	ii = jj;
	    	jj = kk;
	    	kk = mm;
	    	ll = (53 * ll + 1) % Q;
	    	if (((ll * mm) % 64) >= 32)
				s += bit;
	    	bit += bit;
		}
		u[elt] = s;
	}
}
	
double get_next_random_number(void)
/*
 * Return a uniformly distributed pseudo random number
 * in the range 0.0 .. 1.0-2**(-24) inclusive.
 * There are 2**24 possible return values.
 * Side-effects the non-local variables: u, c, ni, nj.
 */
{
double uni;
	    
	if (u[ni] < u[nj])
		uni = u[ni] + (1.0 - u[nj]);
	else
		uni = u[ni] - u[nj];
	u[ni] = uni;

	if (ni > 0)
		ni -= 1;
	else
		ni = STATE_SIZE - 1;

	if (nj > 0)
		nj -= 1;
	else
		nj = STATE_SIZE - 1;

	if (c < cd)
		c = c + (cm - cd);
	else
		c = c - cd;

	if (uni < c)
		return (uni + (1.0 - c));
	else
		return (uni - c);
}
