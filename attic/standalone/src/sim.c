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
 * sim.c: functions, special to (un)conditional simulation
 */
#include <stdio.h>
#include <math.h>

#include "matrix2.h"
#include "defs.h"
#include "random.h"
#include "debug.h"
#include "glvars.h" /* gl_nsim */
#include "userio.h"
#include "data.h"
#include "gls.h"
#include "utils.h"
#include "lm.h"
#include "sim.h"

static void simulate_mvn(const double *est, VEC *result, const int *is_datum);
static void simulate_uniform(double *est, VEC *result, int orc);
static unsigned int n_orvc = 0, n_total = 0;

const double *cond_sim(double *est, int dim, METHOD m, int *is_datum, int orc) {
/*
 * distributions come as { est1, var1, est2, var2,..., cov1-1,.. }
 * and return as { sim-1, sim-2, sim-3, ..., sim-dim }
 */
	int i;
	static VEC *result = VNULL;
	static double *sim = NULL;
	static int max_dim = -1;

	assert(dim > 0);
	assert(est != NULL);
	assert(dim >= 1);
	assert(is_simulation(m));
	assert(est != NULL);
	assert(is_datum != NULL);

	if (dim > max_dim) {
		sim = (double *) erealloc(sim, dim * sizeof(double));
		max_dim = dim;
	}
	result = v_resize(result, dim);

	for (i = 0; i < dim; i++)
		is_datum[i] = (fabs(est[2 * i + 1]) < gl_zero); /* `zero' variance */

	if (m == GSI) /* generate zero mean correlated noise: */
		simulate_mvn(est, result, is_datum); 
	else { /* m == IS */
		correct_orv(est, dim, orc); /* affects est */
		simulate_uniform(est, result, orc);
	}
	for (i = 0; i < dim; i++)
		sim[i] = result->ve[i]; 
	return sim;
}

static void simulate_mvn(const double *est, VEC *result, const int *is_datum) {
	static MAT *M = MNULL;
	static VEC *ind = VNULL, *sim = VNULL;
	static PERM *p = PNULL;
	int i, j;
	volatile int dim;

	p = px_resize(p, result->dim);
	for (i = dim = 0; i < result->dim; i++) {
		if (is_datum[i] == 0) /* non-``zero'' variance -- do simulate */
			p->pe[dim++] = i; /* index of simulation point */
	}
	p->size = dim;
	/*
	 * now dim is the number of pos. variances,
	 * p points their position
	 */
	M = m_resize(M, dim, dim);
	for (i = 0; i < dim; i++) {
		M->me[i][i] = est[2 * p->pe[i] + 1]; /* variances on diagonal */
		for (j = 0; j < i; j++) /* off-diagonal: covariances */
			M->me[j][i] = M->me[i][j] =
				est[2 * result->dim + LTI2(p->pe[j],p->pe[i])];
	}
	if (DEBUG_COV) {
		printlog("# simulation covariance matrix:\n");
		m_logoutput(M);
	}
	/* decompose M: */
	M = CHfactor(M);
	if (DEBUG_COV) {
		printlog("# decomposed error covariance matrix:\n");
		m_logoutput(M);
	}
	/* zero upper triangle: */
	for (i = 0; i < M->m; i++) 
		for (j = i + 1; j < M->m; j++)
			M->me[i][j] = 0.0;
	/* make ind a iid N(0,1) vector */
	ind = v_resize(ind, dim);
	for (i = 0; i < dim; i++)
		ind->ve[i] = r_normal(); /* generate N(0,1) independent samples */
	/* make MVN */
	sim = v_resize(sim, dim);
	sim = mv_mlt(M, ind, sim); /* create zero mean correlated noise */
	if (DEBUG_COV) {
		printlog("# correlated noise vector:\n");
		v_logoutput(sim);
	}
	/* fill result vector: */
	for (i = j = 0; i < result->dim; i++) {
		if (j < dim && i == p->pe[j]) { /* simulated */
			result->ve[i] = est[2 * i] + sim->ve[j];
			j++;
		} else 
			result->ve[i] = est[2 * i];
	}
	if (DEBUG_COV) {
		printlog("\n# simulated values:\n");
		if (is_datum != NULL) {
			for (i = 0; i < result->dim; i++) {
				printlog("%g (%s)\n", result->ve[i], 
					is_datum[i] ? "datum point" : "simulated");
			}
		} else {
			for (i = 0; i < result->dim; i++)
				printlog(" %g", result->ve[i]);
			printlog("\n");
		}
	}
}

static void simulate_uniform(double *est, VEC *result, int orc) {
/*
 * depending on gl_is_pdf, assume indicative (1) or cumulative (0)
 * densitiy function in est;
 * sum of all estimates should equal about 1, which is not
 * the case when cumulative, when last estimate should be
 * close to 1.
 */
	int i, hit;
	double cdf, rn;
	static double *pdf = NULL;

	if (result->dim == 1) {
		result->ve[0] = 1.0 * (r_uniform() < est[0]);
		if (DEBUG_ORDER && (est[0] < 0.0 || est[0] > 1.0))
			pr_warning("order relation violation: P %g not in [0,1]\n", est[0]);
		return;
	}
	if (pdf == NULL)
		pdf = (double *) emalloc(result->dim * sizeof(double));
	for (i = 0; i < result->dim; i++)
		pdf[i] = est[2 * i]; /* local copy */
	if (orc == 4) /* cumulative indicators: make raw pdf */
		for (i = 1; i < result->dim; i++)
			pdf[i] -= pdf[i-1]; /* make pdf from cdf */
	rn = r_uniform();
	hit = 0;
	cdf = pdf[0];
	while (rn > cdf) {
		hit++;
		if (hit < result->dim)
			cdf += pdf[hit];
		else
			break; /* fix hit at this interval */
	}
	/* hit now denotes the class in [ 0 ... dim ] in which rn falls */
	for (i = 0; i < result->dim; i++)
		if (orc <= 3)
			result->ve[i] = (hit == i ? 1.0 : 0.0);
		else
			result->ve[i] = (hit <= i ? 1.0 : 0.0); 
}

void correct_orv(double *est, int n_vars, int orc) {
/*
 * does ``order relation violations corrections'' (acc. to GSLIB pages 77-81);
 * if IS_PDF: correct < 0, > 1, sum <= 1;
 * if not IS_PDF: average upward/downward correction as in GSLIB's ik3d.
 */
 	int i;
 	static int violation = -1, size = -1;
 	static double *down = NULL, *up = NULL, *ori = NULL;
 	double sum = 0.0;

 	if (down == NULL || size < n_vars) {
 		down = (double *) erealloc(down, n_vars * sizeof(double));
 		up = (double *) erealloc(up, n_vars * sizeof(double));
 		ori = (double *) erealloc(ori, n_vars * sizeof(double));
		size = n_vars;
	}

	/* save original */
	for (i = 0; i < n_vars; i++)
		ori[i] = est[2 * i];

 	if (orc <= 3) {
 		for (i = 0; i < n_vars; i++) {
 			est[2 * i] = MAX(est[2 * i], 0.0);
 			est[2 * i] = MIN(est[2 * i], 1.0);
 			sum += est[2 * i];
 		}
 		if (orc == 3 && sum != 1.0) {
			if (DEBUG_ORDER)
 				printlog("sum!=1: ");
 			for (i = 0; i < n_vars; i++)
 				est[2 * i] /= sum;
 		} else if (orc == 2 && sum > 1.0) {
			if (DEBUG_ORDER)
 				printlog("sum>1: ");
 			for (i = 0; i < n_vars; i++)
 				est[2 * i] /= sum;
		}
 	} else { /* orc == 4: cdf; do the GSLIB upward/downward averaging: */
 		/* upward correction: */
 		up[0] = MAX(0.0, MIN(1.0, est[0]));
 		for (i = 1; i < n_vars; i++) /* don't go down && stay < 1 */
 			up[i] = MAX(up[i-1], MIN(1.0, est[2 * i])); 
 		/* downward correction: */
 		down[n_vars-1] = MAX(0.0, MIN(1.0, est[2 * (n_vars-1)]));
 		for (i = n_vars - 2; i >= 0; i--) /* don't go up && stay > 0 */
 			down[i] = MIN(down[i+1], MAX(0.0, est[2 * i])); 
 		/* average upward/downward: */
 		for (i = 0; i < n_vars; i++)
 			est[2 * i] = 0.5 * (down[i] + up[i]);
 	}

	if (n_total == 0) { /* first time */
		if (DEBUG_ORDER) {
			printlog(
	"order relation violation:\n(before correction) --> (after correction)\n");
#ifndef USING_R
 			atexit(print_orvc);
#endif
		}
	}
 	n_total++;

	for (i = violation = 0; !violation && i < n_vars; i++)
 		if (ori[i] != est[2 * i])
 			violation = 1;

 	n_orvc += violation;

	if (DEBUG_ORDER) {
		if (violation) { /* print the order correction */
			for (i = 0; i < n_vars; i++)
 				printlog("%g ", ori[i]);
			printlog("--> ");
			for (i = 0; i < n_vars; i++)
 				printlog("%g ", est[2 * i]);
			printlog("\n");
		}
	}
}

void print_orvc(void) {
	if (n_total > 0) {
		if (n_orvc > 0)
			printlog(
			"number of corrected order relation violations: %u of %u (%.1f%%)\n",
			n_orvc, n_total, (100.0 * n_orvc)/n_total);
		else
			printlog("no order relation violations\n");
		n_orvc = 0; 
		n_total = 0;
	}
}
