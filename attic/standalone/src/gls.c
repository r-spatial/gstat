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
gls.c: module for multivariable generalized least squares prediction
given a linear model Y = X Beta + e, gls() calculates from
the selection d[0]->sel,...,d[n_vars-1]->sel, the multivariable
 BLUE estimation of x0 beta (and Cov(x0[i] Beta, x0[j] Beta)), or
 BLUP of Y(x0) (multivariable universal kriging) (and MSPE Y(x0))
 BLP of Y(x0) (multivariable simple kriging) (and MSPE Y(x0))
 (BLUE: Generalized Least Squares (GLS) best linear unbiased estimate;
 BLUP: GLS best linear unbiased prediction; BLP: GLS best linear prediction)

 UPDATE: only update y, beta (BLUP);

References:
a. Christensen, Ronald, 1991. Linear Models for Multivariate, Time Series,
  and Spatial Data. Springer Verlag, New York. (Ch. VI.2, VI.3)

b. Ver Hoef, Jay M., Noel Cressie, 1993. Multivariable Spatial Prediction.
  Mathematical Geology, 25 (2), pp. 219--240. [[[Note: eq. (18) should have
  a `+' between SIGMA_0 and C'SIGMA-1 C ]]]

Notation (variables and comment):
  v'  the transpose of (matrix or vector) v
  C-1 the inverse of matrix v (`Cinv' in a variable name)
  n   total number of observations in data selections (rows_C)
  m   total number of variables [n_vars]
  p   sum of number of X columns (for non-empty data selections) [cols of X]
  y   data vector (y[0]',...,y[n_vars-1]')', y[i] is the i-th variable [n x 1]
  C   (generalized) covariance matrix, C[i][j] = Cov(element(y,i),element(y,j))
      (element(y,i) is the i-th element in y, not the i-th variable) [n x n]
  X   design matrix [n x p]; 
beta  parameter vector [p x 1]: E(y) = X beta
  C0  (generalized) covariance matrix with y(x0):
      C0[i][j] = Cov(element(y,i),y[j](x0)) [n x m]
  X0  the value that X would have at location "where" [p x m]
      (Note: Christensen's x0 is Ver Hoefs X0'; I use x0: e.g. x0'beta).
*/

#include <stdio.h>
#include <math.h>
#include "config.h"
#include "matrix2.h"
#ifdef HAVE_SPARSE
# include "sparse2.h"
#endif
#include "defs.h"
#include "data.h"
#include "utils.h"
#include "select.h"
#include "lm.h"
#include "vario.h"
#include "vario_io.h"
#include "glvars.h"
#include "userio.h"
#include "plot.h"
#include "debug.h"
#include "gls.h"

static void fill_est(DATA **d, VEC *blup, MAT *MSPE, int n_vars, double *est);
static void debug_result(VEC *blup, MAT *MSPE, enum GLS_WHAT pred);
static VEC *get_mu(VEC *mu, const VEC *y, DATA **d, int nvars);
static MAT *get_corr_mat(MAT *C, MAT *R);

#ifndef USING_R
static void plot_weights(DATA **d, int nvars, DPOINT *where, MAT *weights);
#endif

#define M_DEBUG(a,b) { if (DEBUG_COV){printlog("\n# %s:\n", b); \
	m_logoutput(a);}}
#define V_DEBUG(a,b) {if (DEBUG_COV){printlog("\n# %s:\n", b); \
	v_logoutput(a);}}
#define UPDATE_BLP (pred == UPDATE && last_pred == GLS_BLP)
#define UPDATE_BLUP (pred == UPDATE && last_pred == GLS_BLUP)

#ifdef HAVE_SPARSE
void    sp_logoutput(SPMAT *A);
#define SM_DEBUG(a,b) { if (DEBUG_COV){printlog("\n# %s:\n", b); \
	sp_logoutput(a);}}
#endif

/*
static MAT *convert_vmuC(MAT *C, DATA *d);
static MAT *convert_vmuC0(MAT *C0, DATA *d, double v);
static MAT *convert_vmuC00(MAT *MSPE, double v);
*/

static void convert_C(MAT *C, VEC *mu, double (*fn)(double));
static void convert_C0(MAT *C0, VEC *mu, VEC *mu0, double (*fn)(double));
static void convert_MSPE(MAT *MSPE, VEC *mu0, double (*fn)(double));

typedef struct {
	MAT *C,        /* (Generalized) Covariance matrix */
		*X,        /* design matrix, y = X beta + e */
		*CinvX,    /* C-1 X */
		*XCinvX;   /* X' C-1 X */
#ifdef HAVE_SPARSE
	SPMAT *spC;    /* sparse version of C */
#else
	void *spC;
#endif
	VEC *y,        /* measurement vector */
		*mu,       /* mu vector, E(y) */
		*mu0,      /* mu at loc x0 */
		*beta;     /* parameter vector */
} GLM; /* structure is locally defined, will be held in void *glm in DATA */

static GLM *new_glm(void);

/*
 * n_vars is the number of variables to be considered,
 * d is the data array of variables d[0],...,d[n_vars-1],
 * pred determines which estimate is required: BLUE, BLUP, or BLP
 */
void gls(DATA **d /* pointer to DATA array */,
		int n_vars, /* length of DATA array (to consider) */
		enum GLS_WHAT pred, /* what type of prediction is requested */
		DPOINT *where, /* prediction location */
		double *est /* output: array that holds the predicted values and variances */)
{
	GLM *glm = NULL; /* to be copied to/from d */
	static MAT *X0 = MNULL, *C0 = MNULL, *MSPE = MNULL, *CinvC0 = MNULL,
		*Tmp1 = MNULL, *Tmp2 = MNULL, *Tmp3, *R = MNULL;
	static VEC *blup = VNULL, *tmpa = VNULL, *tmpb = VNULL;
	volatile unsigned int i, rows_C;
	unsigned int j, k, l = 0, row, col, start_i, start_j, start_X, global,
		one_nbh_empty;
	VARIOGRAM *v = NULL;
	static enum GLS_WHAT last_pred = GLS_INIT; /* the initial value */
	double c_value, *X_ori;

	if (d == NULL) { /* clean up */
		if (X0 != MNULL) M_FREE(X0); 
		if (C0 != MNULL) M_FREE(C0);
		if (MSPE != MNULL) M_FREE(MSPE);
		if (CinvC0 != MNULL) M_FREE(CinvC0);
		if (Tmp1 != MNULL) M_FREE(Tmp1);
		if (Tmp2 != MNULL) M_FREE(Tmp2);
		if (Tmp3 != MNULL) M_FREE(Tmp3);
		if (R != MNULL) M_FREE(R);
		if (blup != VNULL) V_FREE(blup);
		if (tmpa != VNULL) V_FREE(tmpa);
		if (tmpb != VNULL) V_FREE(tmpb);
		last_pred = GLS_INIT;
		return;
	}
#ifndef HAVE_SPARSE
	if (gl_sparse) {
		pr_warning("sparse matrices not supported: compile with --with-sparse");
		gl_sparse = 0;
	}
#endif

	if (DEBUG_COV) {
		printlog("we're at %s X: %g Y: %g Z: %g\n",
			IS_BLOCK(where) ? "block" : "point",
			where->x, where->y, where->z);
	}

	if (pred != UPDATE) /* it right away: */
		last_pred = pred;

	assert(last_pred != GLS_INIT);

	if (d[0]->glm == NULL) { /* allocate and initialize: */
		glm = new_glm();
		d[0]->glm = (void *) glm;
	} else
		glm = (GLM *) d[0]->glm;

	glm->mu0 = v_resize(glm->mu0, n_vars);
	MSPE = m_resize(MSPE, n_vars, n_vars);
	if (pred == GLS_BLP || UPDATE_BLP) {
		X_ori = where->X;
		for (i = 0; i < n_vars; i++) { /* mu(0) */
			glm->mu0->ve[i] = calc_mu(d[i], where);
			blup = v_copy(glm->mu0, v_resize(blup, glm->mu0->dim));
			where->X += d[i]->n_X; /* shift to next x0 entry */
		}
		where->X = X_ori; /* ... and set back */
		for (i = 0; i < n_vars; i++) { /* Cij(0,0): */
			for (j = 0; j <= i; j++) {
				v = get_vgm(LTI(d[i]->id,d[j]->id));
				MSPE->me[i][j] = MSPE->me[j][i] = COVARIANCE0(v, where, where, d[j]->pp_norm2);
			}
		}
		fill_est(NULL, blup, MSPE, n_vars, est); /* in case of empty neighbourhood */
	}
	/* xxx */
	/*
	logprint_variogram(v, 1);
	*/

/* 
 * selection dependent problem dimensions: 
 */
	for (i = rows_C = 0, one_nbh_empty = 0; i < n_vars; i++) {
		rows_C += d[i]->n_sel;
		if (d[i]->n_sel == 0)
			one_nbh_empty = 1;
	}

	if (rows_C == 0 /* all selection lists empty */
			|| one_nbh_empty == 1) { /* one selection list empty */
		if (pred == GLS_BLP || UPDATE_BLP)
			debug_result(blup, MSPE, pred);
		return;
	}

	for (i = 0, global = 1; i < n_vars && global; i++)
		global = (d[i]->sel == d[i]->list 
				&& d[i]->n_list == d[i]->n_original
				&& d[i]->n_list == d[i]->n_sel);

/*
 * global things: enter whenever (a) first time, (b) local selections or
 * (c) the size of the problem grew since the last call (e.g. simulation)
 */
	if ((glm->C == NULL && glm->spC == NULL) || !global || rows_C > glm->C->m) {
/* 
 * fill y: 
 */
		glm->y = get_y(d, glm->y, n_vars);

		if (pred != UPDATE) {
			if (! gl_sparse) {
				glm->C = m_resize(glm->C, rows_C, rows_C);
				m_zero(glm->C);
			} 
#ifdef HAVE_SPARSE
			else {
				if (glm->C == NULL) {
					glm->spC = sp_get(rows_C, rows_C, gl_sparse);
					/* d->spLLT = spLLT = sp_get(rows_C, rows_C, gl_sparse); */
				} else {
					glm->spC = sp_resize(glm->spC, rows_C, rows_C);
					/* d->spLLT = spLLT = sp_resize(spLLT, rows_C, rows_C); */
				}
				sp_zero(glm->spC);
			} 
#endif
			glm->X = get_X(d, glm->X, n_vars);
			M_DEBUG(glm->X, "X");
			glm->CinvX = m_resize(glm->CinvX, rows_C, glm->X->n);
			glm->XCinvX = m_resize(glm->XCinvX, glm->X->n, glm->X->n);
			glm->beta = v_resize(glm->beta, glm->X->n);
			for (i = start_X = start_i = 0; i < n_vars; i++) { /* row var */
				/* fill C, mu: */
				for (j = start_j = 0; j <= i; j++) { /* col var */
					v = get_vgm(LTI(d[i]->id,d[j]->id));
					for (k = 0; k < d[i]->n_sel; k++) { /* rows */
						row = start_i + k;
						for (l = 0, col = start_j; col <= row && l < d[j]->n_sel; l++, col++) {
							if (pred == GLS_BLUP)
								c_value = GCV(v, d[i]->sel[k], d[j]->sel[l]);
							else
								c_value = COVARIANCE(v, d[i]->sel[k], d[j]->sel[l]);
							/* on the diagonal, if necessary, add measurement error variance */
							if (d[i]->colnvariance && i == j && k == l)
								c_value += d[i]->sel[k]->variance;
							if (! gl_sparse)
								glm->C->me[row][col] = c_value;
#ifdef HAVE_SPARSE
							else {
								if (c_value != 0.0)
									sp_set_val(glm->spC, row, col, c_value);
							} 
#endif
						} /* for l */
					} /* for k */
					start_j += d[j]->n_sel;
				} /* for j */
				start_i += d[i]->n_sel;
				if (d[i]->n_sel > 0)
					start_X += d[i]->n_X - d[i]->n_merge;
			} /* for i */

			/*
			if (d[0]->colnvmu)
				glm->C = convert_vmuC(glm->C, d[0]);
			*/
			if (d[0]->variance_fn) {
				glm->mu = get_mu(glm->mu, glm->y, d, n_vars);
				convert_C(glm->C, glm->mu, d[0]->variance_fn);
			}

			if (DEBUG_COV && pred == GLS_BLUP)
				printlog("[using generalized covariances: max_val - semivariance()]");
			if (! gl_sparse) {
				M_DEBUG(glm->C, "Covariances (x_i, x_j) matrix C (lower triangle only)");
			}
#ifdef HAVE_SPARSE
			else {
				SM_DEBUG(glm->spC, "Covariances (x_i, x_j) sparse matrix C (lower triangle only)")
			}
#endif
/* check for singular C: */
			if (! gl_sparse && gl_cn_max > 0.0) {
				for (i = 0; i < rows_C; i++) /* row */ 
					for (j = i+1; j < rows_C; j++) /* col > row */
						glm->C->me[i][j] = glm->C->me[j][i]; /* fill symmetric */
				if (is_singular(glm->C, gl_cn_max)) {
					pr_warning("Covariance matrix (nearly) singular at location [%g,%g,%g]: skipping...",
						where->x, where->y, where->z);
					m_free(glm->C); glm->C = MNULL; /* assure re-entrance if global */
					return;
				}
			}
/* 
 * factorize C: 
 */
			if (! gl_sparse)
				LDLfactor(glm->C);
#ifdef HAVE_SPARSE
			else {
				sp_compact(glm->spC, 0.0);
				spCHfactor(glm->spC);
			}
#endif
		} /* if (pred != UPDATE) */
		if (pred != GLS_BLP && !UPDATE_BLP) { /* C-1 X and X'C-1 X, beta */
/* 
 * calculate CinvX: 
 */
    		tmpa = v_resize(tmpa, rows_C);
    		for (i = 0; i < glm->X->n; i++) {
				tmpa = get_col(glm->X, i, tmpa);
				if (! gl_sparse)
					tmpb = LDLsolve(glm->C, tmpa, tmpb);
#ifdef HAVE_SPARSE
				else
					tmpb = spCHsolve(glm->spC, tmpa, tmpb);
#endif
				set_col(glm->CinvX, i, tmpb);
			}
/* 
 * calculate X'C-1 X: 
 */
			glm->XCinvX = mtrm_mlt(glm->X, glm->CinvX, glm->XCinvX); /* X'C-1 X */
			M_DEBUG(glm->XCinvX, "X'C-1 X");
			if (gl_cn_max > 0.0 && is_singular(glm->XCinvX, gl_cn_max)) {
				pr_warning("X'C-1 X matrix (nearly) singular at location [%g,%g,%g]: skipping...",
					where->x, where->y, where->z);
				m_free(glm->C); glm->C = MNULL; /* assure re-entrance if global */
				return;
			}
			m_inverse(glm->XCinvX, glm->XCinvX);
/* 
 * calculate beta: 
 */
			tmpa = vm_mlt(glm->CinvX, glm->y, tmpa); /* X'C-1 y */
			glm->beta = vm_mlt(glm->XCinvX, tmpa, glm->beta); /* (X'C-1 X)-1 X'C-1 y */
			V_DEBUG(glm->beta, "beta");
			M_DEBUG(glm->XCinvX, "Cov(beta), (X'C-1 X)-1");
			M_DEBUG(R = get_corr_mat(glm->XCinvX, R), "Corr(beta)");
		} /* if pred != GLS_BLP */
	} /* if redo the heavy part */

	if (pred != GLS_BLP && !UPDATE_BLP) { /* but BLUE or BLUP */
		X0 = get_X0(d, X0, where, n_vars);
		M_DEBUG(X0, "X0 (X values at prediction location x0)");
		blup = vm_mlt(X0, glm->beta, blup); /* X0' beta = beta'X0 -> vm_ */
		if (pred == GLS_BLUP)
			V_DEBUG(blup, "BLUE(mu), E(y(x0)) = X0'beta");
	}

	if (pred == GLS_BLUE) { /* we did the blue, it's in blup */
/* 
 * BLUE = X0 beta; Cov(X0 beta)= X0'(X'C-1X)-1 X0 
 */
		Tmp1 = mtrm_mlt(X0, glm->XCinvX, Tmp1);
		m_mlt(Tmp1, X0, MSPE); /* X0'(X'C-1X)-1 X0 */
		fill_est(d, blup, MSPE, n_vars, est);
		debug_result(blup, MSPE, pred);
		return; /* Quit function */
	}

/* 
 * now the part that's got to be done every time, for x0 and where change: 
 * resize matrices (BLP, BLUP): 
 */
	C0 = m_resize(C0, rows_C, n_vars);
	CinvC0 = m_resize(CinvC0, rows_C, n_vars);
/* 
 * fill C0: 
 */
	for (i = 0; i < n_vars; i++) { /* cols */
		for (j = 0, start_j = 0; j < n_vars; j++) { /* rows */
			v = get_vgm(LTI(d[i]->id, d[j]->id));
			for (k = 0; k < d[j]->n_sel; k++) {
				if (pred == GLS_BLUP)
					C0->me[start_j+k][i] = GCV0(v, d[j]->sel[k], where, d[j]->pp_norm2);
				else
					C0->me[start_j+k][i] = COVARIANCE0(v, d[j]->sel[k], where, d[j]->pp_norm2);
			}
			start_j += d[j]->n_sel;
		}
	}

	/*
	if (d[0]->colnvmu) {
		X0 = get_X0(d, X0, where, n_vars);
		C0 = convert_vmuC0(C0, d[0], X0->me[0][0]);
	}
	*/
	if (d[0]->variance_fn)
		convert_C0(C0, glm->mu, glm->mu0, d[0]->variance_fn);

	M_DEBUG(C0, "Covariances (x_i, x_0), C0");

/* 
 * calculate CinvC0: 
 */
   	tmpa = v_resize(tmpa, rows_C);
   	for (i = 0; i < n_vars; i++) {
		tmpa = get_col(C0, i, tmpa);
		if (! gl_sparse)
			tmpb = LDLsolve(glm->C, tmpa, tmpb);
#ifdef HAVE_SPARSE
		else
			tmpb = spCHsolve(glm->spC, tmpa, tmpb);
#endif
		set_col(CinvC0, i, tmpb);
	}
	M_DEBUG(CinvC0, "C-1 C0");

	if (pred == GLS_BLP || UPDATE_BLP) {
/* 
 * BLP = mu_0 + C0'C-1 (y-mu_i) 
 */
		V_DEBUG(glm->y, "data values y");

		if (DEBUG_COV) {
			printlog("beta is:\n");
			for (i = 0; i < n_vars; i++)
				for (j = 0; j < d[i]->beta->size; j++)
					printlog("%g%s", d[i]->beta->val[j], 
						j == d[i]->beta->size - 1 ? "\n" : " ");
		}

		glm->mu = get_mu(glm->mu, glm->y, d, n_vars);
		V_DEBUG(glm->mu, "mean values (mu_i)");

		/* y - mu_i: */
		v_sub(glm->y, glm->mu, tmpa);
		V_DEBUG(tmpa, "Residual vector (y - mu_i)");

		tmpb = vm_mlt(CinvC0, tmpa, tmpb); /* C0'C-1 (y-mu_i) */
		V_DEBUG(tmpb, "Weighted res. vector, C0'C-1 (y-mu_i)");
		v_add(blup, tmpb, blup); /* mu_0 + C0'C-1 (y-mu_i) */
/* 
 * MSPE = C(0,0) - C0'C-1 C0 
 */
		Tmp2 = mtrm_mlt(CinvC0, C0, Tmp2);  /* C0'C-1 C0 */

		/*
		if (d[0]->colnvmu)
			MSPE = convert_vmuC00(MSPE, X0->me[0][0]);
		*/
		if (d[0]->variance_fn)
			convert_MSPE(MSPE, glm->mu0, d[0]->variance_fn);

		m_sub(MSPE, Tmp2, MSPE); /* C(0,0) - C0'C-1 C0 */
		fill_est(d, blup, MSPE, n_vars, est);
		debug_result(blup, MSPE, pred);
		return; /* Quit function */
	}

/* 
 * GLS_BLUP, universal kriging estimate remains: 
 */
	tmpa = mv_mlt(glm->X, glm->beta, tmpa); /* X beta */
	tmpb = v_sub(glm->y, tmpa, tmpb); /* y - X beta */
	tmpa = vm_mlt(CinvC0, tmpb, tmpa); /* c0'C-1 (Y - X beta) */
	blup = v_add(blup, tmpa, blup); /* x0 beta + c0'C-1 (Y - X beta) */
/* 
 * universal kriging MSPE: 
 * (a) Cov_ij(0,0): 
 */
	for (i = 0; i < n_vars; i++) {
		for (j = 0; j <= i; j++) {
			v = get_vgm(LTI(d[i]->id, d[j]->id));
			MSPE->me[i][j] = MSPE->me[j][i] = GCV0(v, where, where, d[j]->pp_norm2);
		}
	}
	M_DEBUG(MSPE, "[a] Cov_ij(B,B) or Cov_ij(0,0)");
/* 
 * (c) (x0-X'C-1 c0)'(X'C-1X)-1 (x0-X'C-1 c0): 
 */
	Tmp1 = mtrm_mlt(glm->CinvX, C0, Tmp1); /* X'C-1 c0 */
	Tmp1 = m_sub(X0, Tmp1, Tmp1); /* (x0 - X'C-1 c0) */
	Tmp2 = m_copy(Tmp1, m_resize(Tmp2, Tmp1->m, Tmp1->n));
	Tmp3 = m_mlt(glm->XCinvX, Tmp1, Tmp3); /* (X'C-1 X)-1 (x0 - X'C-1 c0) */
	Tmp1 = mtrm_mlt(Tmp2, Tmp3, Tmp1);
	M_DEBUG(Tmp1, "[c] (x0-X'C-1 c0)'(X'C-1 X)-1(x0-X'C-1 c0)");
/* 
 * (b) c0'C-1 c0: 
 */
	Tmp2 = mtrm_mlt(C0, CinvC0, Tmp2);
	M_DEBUG(Tmp2, "[b] c0'C-1 c0");
/* 
 * (a - b + c) =
 * Cov_ij(0,0) - c0'C-1 c0 + (x0-X'C-1 c0)'(X'C-1 X)-1 (x0-X'C-1 c0): 
 */
	m_sub(MSPE, Tmp2, MSPE); /* a-b */
	m_add(MSPE, Tmp1, MSPE); /* +c */
/* 
 * done: 
 */
	fill_est(d, blup, MSPE, n_vars, est);
	debug_result(blup, MSPE, pred);
	if (DEBUG_COV || plotfile) { /* calculate kriging weights explicitly: */
		/* Tmp3' * glm->CinvX' + CinvC0 */
		Tmp1 = m_mlt(glm->CinvX, Tmp3, Tmp1);
		Tmp2 = m_add(Tmp1, CinvC0, Tmp2);
		M_DEBUG(Tmp2, "kriging weights");
#ifndef USING_R
		if (plotfile)
			plot_weights(d, n_vars, where, Tmp2);
#endif
		if (DEBUG_COV)
			printlog("\n\n");
	}
	return;
}

static void fill_est(DATA **d, VEC *blup, MAT *MSPE, int n_vars, double *est)
{
	int i, j, n_filled;
	static IVEC *v = IVNULL;

	if (n_vars == 1) {
		est[0] = blup->ve[0];
		est[1] = MSPE->me[0][0];
		return;
	}
	v = iv_resize(v, n_vars);
	if (d == NULL) { /* GLS_BLP, initializing */
		for (i = 0; i < n_vars; i++)
			v->ive[i] = i;
		n_filled = n_vars;
	} else { /* n_vars > 1: avoid possibly empty variables -> NA */
		for (i = j = 0; i < n_vars; i++) {
			if (d[i]->n_sel > 0) {
				v->ive[j] = i;
				j++;
			}
		}
		n_filled = j;
	}
	for (i = 0; i < n_filled; i++) { /* only adress non-empty variables */
		est[2 * v->ive[i]] = blup->ve[v->ive[i]];
		est[2 * v->ive[i] + 1] = MSPE->me[v->ive[i]][v->ive[i]];
		for (j = 0; j < i; j++)
			est[2 * n_vars + LTI2(v->ive[i], v->ive[j])] =
				MSPE->me[v->ive[i]][v->ive[j]];
	}
	return;
}

static void debug_result(VEC *blup, MAT *MSPE, enum GLS_WHAT pred) {

	if (! DEBUG_COV)
		return;
	switch (pred) {
	case GLS_BLP:
		V_DEBUG(blup, "Best Linear Predictor");
		M_DEBUG(MSPE, "Prediction Covariances");
		break;
	case GLS_BLUE:
		V_DEBUG(blup, "Best Linear Unbiased Estimate (X0'beta)");
		M_DEBUG(MSPE, "Estimation Covariances, Cov(X0'beta)");
		break;
	case GLS_BLUP:
		V_DEBUG(blup, "Best Linear Unbiased Predictor");
		M_DEBUG(MSPE, "MSPE ([a]-[b]+[c])");
		break;
	case UPDATE:
		V_DEBUG(blup, "Updated predictor");
		M_DEBUG(MSPE, "MSPE (updated)");
		break;
	case GLS_INIT:
		ErrMsg(ER_IMPOSVAL, "invalid value for pred");
		break;
	}
}

double *make_gls(DATA *d, int calc_residuals) {
/* 
 * if calc_residuals == 0, return value is allocated, but not freed 
 */
	int i, j, size;
	double *est = NULL;
	DATA **data;
	GLM *glm;

	glm = (GLM *) d->glm;
	if (glm == NULL) {
		data = get_gstat_data();
		glm = (GLM *) data[0]->glm;
	}
	if (glm && (glm->C || glm->spC)) { /* renew: variogram may have changed */
		if (! gl_sparse)
			m_free((MAT *) glm->C);
#ifdef HAVE_SPARSE
		else
			sp_free((SPMAT *) glm->spC);
#endif
		glm->C = MNULL;
		glm->spC = NULL;
	} 
	select_at(d, NULL); /* where == NULL --> global selection */
	if (calc_residuals) {
		est = (double *) emalloc(get_n_outfile() * sizeof(double));
		for (i = 0; i < d->n_list; i++) {
			gls(&d, 1, GLS_BLUE, d->list[i], est);
			glm = (GLM *) d->glm;
			d->list[i]->attr = glm->y->ve[i] - est[0];
		}
		efree(est);
		est = NULL;
	} else { /* no residuals -- return beta & Cov(beta) */
		size = d->n_X * (1 + d->n_X);
		est = (double *) emalloc(size * sizeof(double));
		/* fill the GLM stuff: */
		gls(&d, 1, GLS_BLUE, d->list[0], est); 
		glm = (GLM *) d->glm;
		for (i = 0; i < glm->beta->dim; i++) {
			est[2 * i] = glm->beta->ve[i];
			est[2 * i + 1] = glm->XCinvX->me[i][i];
			for (j = 0; j < i; j++)
				est[2 * glm->beta->dim + LTI2(i,j)] = glm->XCinvX->me[i][j];
		}
	}
	gls(NULL, 0, GLS_INIT, NULL, NULL);
	return est; /* possibly NULL */
}

double *make_gls_mv(DATA **d, int n_vars) {
/* 
 * allocates memory for est (return value) but does not free it 
 */
	int i, j, sum_X, index, size = 0;
	double *est = NULL;
	GLM *glm;
	DPOINT where;

	for (i = sum_X = 0; i < n_vars; i++) {
		select_at(d[i], NULL); /* where == NULL --> global selection */
		sum_X += d[i]->n_X;
	}
	where = *d[0]->list[0];
	where.X = (double *) emalloc(sum_X * sizeof(double)); /* replace */
	for (i = 0; i < sum_X; i++) /* fill with nonsense values: */
		where.X[i] = 0.0; 
	size = sum_X + (sum_X * (sum_X + 1))/2;
	est = (double *) emalloc(size * sizeof(double));
	/* fill the GLM stuff: */
	gls(d, n_vars, GLS_BLUE, &where, est); 
	glm = (GLM *) d[0]->glm;
	assert(glm != NULL);
	for (i = 0; i < glm->beta->dim; i++) {
		assert((2 * i + 1) < size);
		est[2 * i] = glm->beta->ve[i];
		est[2 * i + 1] = glm->XCinvX->me[i][i];
		for (j = 0; j < i; j++) {
			index = 2 * glm->beta->dim + LTI2(i,j);
			assert(index < size);
			est[index] = glm->XCinvX->me[i][j];
		}
	}
	gls(NULL, 0, GLS_INIT, NULL, NULL);
	efree(where.X);
	return est;
}

static GLM *new_glm(void) {
	GLM *glm;

	glm = (GLM *) emalloc(sizeof(GLM));
	glm->X = glm->C = glm->CinvX = glm->XCinvX = MNULL;
	glm->spC = NULL;
	glm->y = glm->mu = glm->mu0 = glm->beta = VNULL;
	return glm;
}

void free_glm(void *v_glm) {
	GLM *glm;

	if (v_glm == NULL)
		return;

	glm = (GLM *) v_glm;
	if (glm->X)
		m_free(glm->X);
	if (glm->C)
		m_free(glm->C);
	if (glm->CinvX)
		m_free(glm->CinvX);
	if (glm->XCinvX)
		m_free(glm->XCinvX);
	if (glm->y)
		v_free(glm->y);
	if (glm->beta)
		v_free(glm->beta);
	if (glm->mu0)
		v_free(glm->mu0);
	if (glm->mu)
		v_free(glm->mu);
	/* EJPXX
	if (glm->mu)
		v_free(glm->mu);
	if (glm->mu0)
		v_free(glm->mu0);
	*/
	free(glm);
}

static void convert_C(MAT *C, VEC *mu, double (*fn)(double)) {
	int i, j;
	double sqrtfni;

	assert(C && mu);
	assert(C->m == mu->dim);

	for (i = 0; i < mu->dim; i++) {
		/* assert(mu->ve[i] >= 0.0); */
		/* be more friendly: */
		if (mu->ve[i] < 0.0)
			ErrMsg(ER_IMPOSVAL, "can not take square root of negative mean values!");
		C->me[i][i] *= fn(mu->ve[i]);
		sqrtfni = sqrt(fn(mu->ve[i]));
		for (j = 0; j < i; j++)
			C->me[i][j] *= sqrtfni * sqrt(fn(mu->ve[j]));
	}
}

static void convert_C0(MAT *C0, VEC *mu, VEC *mu0, double (*fn)(double)) {
	int i, j;
	double sqrtfni;

	assert(C0 && mu && mu0);
	assert(C0->m == mu->dim);
	assert(C0->n == mu0->dim);

	for (i = 0; i < mu->dim; i++) {
		assert(mu->ve[i] >= 0.0);
		sqrtfni = sqrt(fn(mu->ve[i]));
		for (j = 0; j < mu0->dim; j++)
			C0->me[i][j] *= sqrtfni * sqrt(fn(mu0->ve[j]));
	}
}

static void convert_MSPE(MAT *MSPE, VEC *mu0, double (*fn)(double)) {
	int i, j;
	double sqrtfni;

	assert(MSPE && mu0);
	assert(MSPE->m == mu0->dim);

	for (i = 0; i < mu0->dim; i++) {
		assert(mu0->ve[i] >= 0.0);
		MSPE->me[i][i] *= fn(mu0->ve[i]);
		sqrtfni = sqrt(fn(mu0->ve[i]));
		for (j = 0; j < i; j++) {
			MSPE->me[i][j] *= sqrtfni * sqrt(fn(mu0->ve[j]));
			MSPE->me[j][i] = MSPE->me[i][j];
		}
	}
}

static VEC *get_mu(VEC *mu, const VEC *y, DATA **d, int n_vars) {
 	int i, start_j, j;

 	mu = v_resize(mu, y->dim);
	for (i = start_j = 0; i < n_vars; i++) { 
 		for (j = 0; j < d[i]->n_sel; j++)
 			mu->ve[start_j + j] = calc_mu(d[i], d[i]->sel[j]);
		start_j += d[i]->n_sel;
	}
	return mu;
}

static MAT *get_corr_mat(MAT *C, MAT *R) {
	int i, j;

	assert(C);
	assert(C->m == C->n);

	R = m_copy(C, m_resize(R, C->m, C->n));
	for (i = R->m - 1; i >= 0; i--) {
		assert(R->me[i][i] > 0.0);
		for (j = 0; j < i; j++)
			R->me[i][j] /= sqrt(R->me[i][i] * R->me[j][j]);
		for (j = i + 1; j < R->m; j++)
			R->me[i][j] = R->me[j][i];
	}

	for (i = 0; i < R->m; i++)
		R->me[i][i] = 1.0;
	return(R);
}

#ifndef USING_R
static void plot_weights(DATA **d, int nvars, DPOINT *where, MAT *weights) {
	int i, j, k, l, ps;
	static int ps_min = 0, ps_max = 0, *index = NULL /*, gif = 0 */ ;

	if (gl_plotweights == 0) /* only plotfile was specified: */
		gl_plotweights = 1;
	/* plot labels: weights */
	fprintf(plotfile, "# at location (%g, %g):\n", where->x, where->y);
	fprintf(plotfile, "reset; set size ratio -1; set key outside\n");
	switch (gl_plotweights) {
		case 1:
			for (k = 0; k < nvars; k++) {
				for (i = l = 0; i < nvars; i++) {
					fprintf(plotfile, "# variable %d:\n", i);
					for (j = 0; j < d[i]->n_sel; j++)
						fprintf(plotfile, "set label \" %.3f\" at %g,%g left\n",
							weights->me[l++][k], d[i]->sel[j]->x, 
							d[i]->sel[j]->y);
				}
				/* plot command: */
				fprintf(plotfile, "plot ");
				for (i = 0; i < nvars; i++)
					fprintf(plotfile, "'-' title '%s',", name_identifier(i));
				/* plot where: */
				fprintf(plotfile, "'-' title 's_0'\n");
				/* inline data: */
				for (i = 0; i < nvars; i++) {
					for (j = 0; j < d[i]->n_sel; j++)
						fprintf(plotfile, "%g %g\n", d[i]->sel[j]->x, 
							d[i]->sel[j]->y);
					fprintf(plotfile, "e\n");
				}
				/* end & data where: */
				fprintf(plotfile, "%g %g\ne\n", where->x, where->y);
				fprintf(plotfile, "pause -1\nreset\n");
			} /* for k */
			break;
		default: /* for var 0 only: */
			/* fprintf(plotfile, "set term gif; set out 'foo%05d.gif'\n", 
				gif++); */
			assert(weights->m >= d[0]->n_sel);
			if (index == NULL) /* point sizes & sign -> pos/neg weights */
				index = (int *) emalloc(sizeof(int) * d[0]->n_list);
			for (i = 0; i < d[0]->n_list; i++)
				index[i] = 0;
			/* classify points in neighbourhood selection: */
			for (i = 0; i < d[0]->n_sel; i++) {
				ps = (int) floor(fabs(weights->me[i][0]) * gl_plotweights);
				ps++;
				if (weights->me[i][0] < 0)
					ps = -ps;
				if (ps < ps_min)
					ps_min = ps;
				if (ps > ps_max)
					ps_max = ps;
				index[GET_INDEX(d[0]->sel[i])] = ps;
			}
			fprintf(plotfile, "plot ");
			for (i = ps_min; i < 0; i++) 
				fprintf(plotfile, "'-' title '[%g - %g]' lt 2 pt 1 ps %d,", 
					(i+1.0)/gl_plotweights, (1.0 * i)/gl_plotweights, -i);
			for (i = 1; i <= ps_max; i++) 
				fprintf(plotfile, "'-' title '[%g - %g]' lt 1 pt 1 ps %d,", 
					(i-1.0)/gl_plotweights, (1.0 * i)/gl_plotweights, i);
			fprintf(plotfile, "'-' title 'not considered' lt 0 pt 2 ps 1,");
			fprintf(plotfile, "'-' title 's_0' lt 3 pt 3 ps 1\n");
			for (i = ps_min; i < 0; i++) {
				for (j = 0; j < d[0]->n_list; j++)
					if (index[j] == i)
						fprintf(plotfile, "%g %g\n", d[0]->list[j]->x, 
							d[0]->list[j]->y);
				fprintf(plotfile, "e\n");
			}
			for (i = 1; i <= ps_max; i++) {
				for (j = 0; j < d[0]->n_list; j++)
					if (index[j] == i)
						fprintf(plotfile, "%g %g\n", d[0]->list[j]->x, 
							d[0]->list[j]->y);
				fprintf(plotfile, "e\n");
			}
			for (j = 0; j < d[0]->n_list; j++)
				if (index[j] == 0)
					fprintf(plotfile, "%g %g\n", d[0]->list[j]->x, 
						d[0]->list[j]->y);
			fprintf(plotfile, "e\n");
			fprintf(plotfile, "%g %g\ne\n", where->x, where->y);
			break;
	}
	fprintf(plotfile, "\n"); /* newline before next point location */
	return;
}
#endif /* USING_R */

#ifdef HAVE_SPARSE
/* sp_foutput -- output sparse matrix A to file/stream fp */
void    sp_logoutput(SPMAT *A)
{
	int     i, j_idx, m /* , n */;
	SPROW  *rows;
	row_elt *elts;

	printlog("SparseMatrix: ");
	if ( A == SMNULL )
	{
		printlog("*** NULL ***\n");
		error(E_NULL,"sp_foutput");    return;
	}
	printlog("%d by %d\n",A->m,A->n);
	m = A->m;       /* n = A->n; */
	if ( ! (rows=A->row) )
	{
		printlog("*** NULL rows ***\n");
		error(E_NULL,"sp_foutput");    return;
	}

	for ( i = 0; i < m; i++ )
	{
		printlog("row %d: ",i);
		if ( ! (elts=rows[i].elt) )
		{
			printlog("*** NULL element list ***\n");
			continue;
		}
		for ( j_idx = 0; j_idx < rows[i].len; j_idx++ )
		{
			printlog("%d:%-20.15g ",elts[j_idx].col,
							elts[j_idx].val);
			if ( j_idx % 3 == 2 && j_idx != rows[i].len-1 )
				printlog("\n     ");
		}
		printlog("\n");
	}
	printlog("#\n");	/* to stop looking beyond for next entry */
}
#endif
