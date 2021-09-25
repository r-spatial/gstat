/*
 * lm.c: contains routines for linear model y=Xb+e, e independent.
 */
#include <stdio.h>
#include <stddef.h> /* size_t for DJGPP */
#include <float.h> /* DBL_MAX */
#include <math.h> /* sqrt() */

#include "mtrx.h"
#include "defs.h"
#include "userio.h"
#include "data.h"
#include "utils.h"
#include "debug.h"
#include "select.h"
#include "glvars.h"
#include "lm.h"

static void predict_lm(LM *lms, MAT *X0, double *est);
static void create_lm(DATA **d, int nvars);
static VEC *get_weights(DATA **d, VEC *weights, int nvars);
static int get_colX_nr(DATA **d, int var, int x);
static int zero_weights_count(LM *lm);

void pred_lm(DATA **data, int n_vars, DPOINT *where, double *est) {
	int i = 0, global = 1;
	LM *lm;
	static MAT *X0 = MNULL;

	global = 1;
	while (global && i < n_vars) {
		if (data[i]->sel != data[i]->list) /* local selection */
			global = 0; /* and jump out of this loop */
		i++;
	}
	if (data[0]->lm == NULL || !global) {
		create_lm(data, n_vars);
		if (DEBUG_FIT) {
			printlog("at location [%g,%g,%g]:\n",
				where->x, where->y, where->z);
			logprint_lm(data[0], data[0]->lm);
		} 
	} 
	lm = (LM *) data[0]->lm;
	if (lm == NULL || lm->y->dim == 0 || lm->is_singular) {
		for (i = 0; i < n_vars; i++) {
			set_mv_double(&(est[2 * i]));
			set_mv_double(&(est[2 * i + 1]));
		}
		if (lm && lm->is_singular)
			pr_warning("singular X matrix at x[%g], y[%g], z[%g]:",
			where->x, where->y, where->z);
	} else {
		X0 = get_X0(data, X0, where, n_vars);
		if (DEBUG_COV) {
			printlog("#X0 is ");
			m_logoutput(X0);
		}
		predict_lm(lm, X0, est);
	}
	return;
}

double *make_ols(DATA *d) {
/* return value is allocated, but not freed */
	int i, j, size;
	double *est = NULL;
	DATA **data;
	LM *lm;

	lm = (LM *) d->lm;
	if (lm == NULL) {
		data = get_gstat_data();
		lm = (LM *) data[0]->lm;
	}
	select_at(d, NULL); /* where == NULL --> global selection */
	/* return beta & Cov(beta) */
	size = d->n_X * (1 + d->n_X);
	est = (double *) emalloc(size * sizeof(double));
	for (i = 0; i < size; i++)
		set_mv_double(&(est[i]));
	/* fill the LM stuff: */
	create_lm(&d, 1);
	if (DEBUG_FIT)
		logprint_lm(d, d->lm);
	lm = (LM *) d->lm;
	if (lm->is_singular)
		return est; /* all NA's */
	for (i = 0; i < lm->beta->dim; i++) {
		est[2 * i] = lm->beta->ve[i];
		est[2 * i + 1] = ME(lm->Cov, i, i);
		for (j = 0; j < i; j++)
			est[2 * lm->beta->dim + LTI2(i,j)] = ME(lm->Cov, i, j);
	}
	free_lm(d->lm);
	d->lm = NULL;
	return est;
}

LM *init_lm(LM *lm) {
	if (lm == NULL)
		lm = (LM *) emalloc(sizeof(LM));
   	lm->y = VNULL;
   	lm->weights = VNULL;
   	lm->beta = VNULL;
   	lm->Xty = VNULL;
   	lm->X = MNULL;
   	lm->Cov = MNULL;
   	lm->Chol = MNULL;
   	lm->SSReg = lm->MSReg = DBL_MAX;
   	lm->SSErr = lm->MSErr = DBL_MAX;
	lm->is_singular = 0;
   	return lm;
}

void free_lm(LM *lm) {
	if (lm->y)
		v_free(lm->y);
	if (lm->weights)
		v_free(lm->weights);
	if (lm->beta)
		v_free(lm->beta);
	if (lm->Xty)
		v_free(lm->Xty);
	if (lm->X)
		m_free(lm->X);
	if (lm->Chol)
		m_free(lm->Chol);
	if (lm->Cov)
		m_free(lm->Cov);
	efree(lm);
}

static void create_lm(DATA **data, int nvars) {
/*
 * obtain base necessities (y, X, weights), and calculate the rest
 */
	LM *lm;
	int i;

	lm = (LM *) data[0]->lm;
	if (lm == NULL) /* create */
		data[0]->lm = lm = init_lm(NULL);
	lm->X = get_X(data, lm->X, nvars);
	lm->y = get_y(data, lm->y, nvars);
	lm->weights = get_weights(data, lm->weights, nvars);
/* 
 * check for intercept:
 */
 	if (nvars == 1) {
		lm->has_intercept = (data[0]->colX[0] == 0);
		for (i = 1; i < data[0]->n_X && ! lm->has_intercept; i++)
			lm->has_intercept = (data[0]->colX[i] == 0);
	}
/*
 * calculate:
 */
	data[0]->lm = (void *) calc_lm(lm);
	return;
}

void make_residuals_lm(DATA *d) {
/*
 * make (local or global) residuals
 */
	int i;
	double est[2];
	static MAT *X0 = MNULL;

	if (d->is_residual)
		return;
	if (d->beta) {
		/*
		pr_warning("calculating residuals with respect to pre-defined mean"); 
		*/
		for (i = 0; i < d->n_list; i++)
			d->list[i]->attr -= calc_mu(d, d->list[i]);
	} else {
		select_at(d, NULL);
		create_lm(&d, 1);
		if (DEBUG_FIT)
			logprint_lm(d, d->lm);
		for (i = 0; i < d->n_list; i++) {
			X0 = get_X0(&d, X0, d->list[i], 1);
			predict_lm((LM *) d->lm, X0, est);
			d->list[i]->attr -= est[0]; /* y = Xb + e, so e^ = y - Xb^ */
		}
	}
	d->is_residual = 1;
	return;
}

static void predict_lm(LM *lm, MAT *X0, double *est) {
/*
 * make a prediction for x0, store pred. and variance in est[0] and est[1]
 */
	VEC *blup = VNULL;
	MAT *tmp = MNULL, *ans = MNULL;

	if (lm->beta == NULL)
		ErrMsg(ER_IMPOSVAL, "lm->beta NULL: sample too small?");
	blup = vm_mlt(X0, lm->beta, blup); /* X0' beta = beta'X0 -> vm_ */
	/*
	 * Cov(y0^ - y0) = x0'X'X-1 x0 MSErr, or x0'(X'WX)-1 x0 MSErr
	 * for this, solve X'X b = x0 for b, then b = X'X-1 x0; ans = x0'b
	 * MSErr added when prediction is pointwise 
	 */
	tmp = CHsolve(lm->Chol, X0, tmp, PNULL);
	ans = mtrm_mlt(X0, tmp, ans);
	ans = sm_mlt(lm->MSErr, ans, ans);
	for (int i = 0; i < ans->m; i++) {
		est[2 * i] = blup->ve[i];
		est[2 * i + 1] = ME(ans, i, i);
		if (max_block_dimension(1) == 0.0) /* pointwise prediction */
			est[2 * i + 1] += lm->MSErr;
		for (int j = 0; j < i; j++)
			est[2 * ans->m + LTI2(i,j)] = ME(ans, i, j);
	}
	v_free(blup); 
	m_free(tmp); 
	m_free(ans);
	return;
}

MAT *get_X(DATA **d, MAT *X, int nvars) {
	int i, j, k, rows, cols, colX;

	for (i = rows = cols = 0; i < nvars; i++) {
		rows += d[i]->n_sel;
		if (d[i]->n_sel > 0)
			cols += d[i]->n_X - d[i]->n_merge;
	}
	/*
	if (rows <= 0 || cols <= 0)
		ErrMsg(ER_IMPOSVAL, "get_X: size <= 0");
	*/
	X = m_resize(X, rows, cols);
	m_zero(X);
	for (i = rows = 0; i < nvars; i++) { /* i: var block index */
		for (j = 0; d[i]->n_sel && j < d[i]->n_X; j++)  { /* j: column index */
			colX = get_colX_nr(d, i, j);
			for (k = 0; k < d[i]->n_sel; k++) /* k: row index */
				ME(X, rows + k, colX) = d[i]->sel[k]->X[j];
		}
		rows += d[i]->n_sel;
	}
	return X;
}

MAT *get_X0(DATA **d, MAT *X0, DPOINT *where, int nvars) {
	int i, j, start_i, cols, colX;

	for (i = cols = 0; i < nvars; i++) {
		if (d[i]->n_sel > 0)
			cols += d[i]->n_X - d[i]->n_merge;
	}
	X0 = m_resize(X0, cols, nvars);
	m_zero(X0); /* initialize; now fill non-zero entries: */
	for (i = start_i = 0; i < nvars; i++) {
		for (j = 0; d[i]->n_sel && j < d[i]->n_X; j++) { /* for X-column k+.. */
			colX = get_colX_nr(d, i, j);
			ME(X0, colX, i) = where->X[start_i+j]; /* colX is row index: X0' */
		}
		start_i += d[i]->n_X;
	}
	return X0;
}

VEC *get_y(DATA **d, VEC *y, int nvars) {
	int i, j, offset, size;

	for (i = size = 0; i < nvars; i++)
		size += d[i]->n_sel;
	/*
	if (size <= 0)
		ErrMsg(ER_IMPOSVAL, "get_X: size <= 0");
	*/
	y = v_resize(y, size);
	for (j = offset = 0; j < nvars; j++) {
		for (i = 0; i < d[j]->n_sel; i++) 
			y->ve[offset + i] = d[j]->sel[i]->attr;
		offset += d[j]->n_sel;
	}
	return y;
}

static int get_colX_nr(DATA **d, int var, int this_x) {
	int offset_x = 0, colX, i, j;

	for (i = 0; i < var; i++)
		if (d[i]->n_sel)
			offset_x += d[i]->n_X - d[i]->n_merge;
	if (d[var]->n_merge == 0)
		return offset_x + this_x;
	for (i = 0; i < d[var]->n_merge; i++) { 
		if (d[var]->mtbl[i].col_this_X == this_x) { /* hit: merge this one! */
			colX = d[var]->mtbl[i].col_other_X;
			if (d[var]->mtbl[i].to_var > 0)
				for (j = 0; j < d[var]->mtbl[i].to_var - 1; j++)
					colX += d[j]->n_X - d[j]->n_merge;
			return colX;
		} 
	}
	/* so, we're at least at offset_x, but how much should we add? */
	colX = offset_x + this_x; /* i.e. the maximum... */
	for (i = 0; i < d[var]->n_merge; i++)
		for (j = 0; j < this_x; j++)
			if (d[var]->mtbl[i].col_this_X == j)
				colX--; /* ...minus all previously merged entries */
	return colX;
}

static VEC *get_weights(DATA **d, VEC *weights, int nvars) {
	int i, j, size;

	for (i = size = 0; i < nvars; i++) {
		if (d[i]->colnvariance <= 0) /* no weights */
			return VNULL;
		if (d[i]->n_sel > 0)
			size += d[i]->n_sel;
	}
	if (size <= 0)
		return VNULL;
	weights = v_resize(weights, size);
	for (i = size = 0; i < nvars; i++) {
		for (j = 0; j < d[i]->n_sel; j++)
			weights->ve[size + j] = 1.0 / d[i]->sel[j]->variance;
			/* ->variance > 0 is assured in data.c */
		size += d[i]->n_sel;
	}
	return weights;
}

void logprint_lm(DATA *d, LM *lm) {
	double SSTot;
	char line[] = "-----------------------------------------------------------";
	int i, coords = 0;

	if (lm->dfReg <= 0)
		return;
	SSTot = lm->SSReg + lm->SSErr;

	if (d != NULL) {
		printlog("\nmodel: %s = ", d->variable);
		for (i = 0; i < d->n_X; i++) {
			if (i > 0) {
				printlog(" + ");
				if ((i + 2) % 5 == 0)
					printlog("\n");
			}
			printlog("%g", lm->beta->ve[i]);
			if (d->colX[i] > 0)
				printlog(" [col %d]", d->colX[i]);
			if (d->colX[i] < 0) {
				printlog(" %s", POLY_NAME(d->colX[i]));
				coords = 0;
			}
		}
		printlog(" + e\n");
		if (coords)
			printlog(
"(Note: coordinate coefficients apply to normalized coordinates)\n\n");
	}
	printlog("Summary statistics (model %s intercept):\n",
		lm->has_intercept ? "with" : "without");
	printlog("Source            df         SS           MS           F\n");
	printlog("%s\n", line);
	printlog("Regression       %3d %12.6g %12.6g",
		lm->dfReg, lm->SSReg, lm->MSReg);
	if (lm->MSErr <= 0.0)
		printlog("      Inf\n");
	else
		printlog(" %12.6g\n", lm->MSReg/lm->MSErr);
	printlog("Error            %3d %12.6g %12.6g\n",
		lm->dfE, lm->SSErr, lm->MSErr);
	printlog("%s\nTotal, %s %3d %12.6g\n%s\n\n", line,
		lm->has_intercept ? "corrected" : "uncorr.  ",
		lm->dfReg + lm->dfE, SSTot, line);
	/* if (SSTot > 0.0)
		printlog("R2: %g\n", lm->SSReg / SSTot); */

	return;
}

LM *calc_lm(LM *lm) {
/*
 * calculate Chol,Xty,beta,SSErr,SSReg,MSErr^...
 * ASSUMES lm->X, lm->y and optionally lm->weights to be filled!
 */
	double y_mean = 0, w;
	int i, j;
	/* static MAT *QR = MNULL; */
	static VEC *tmp = VNULL;

	if (lm->X == MNULL || lm->y == VNULL)
		ErrMsg(ER_NULL, "calc_lm(): y or X");
	if (lm->X->m != lm->y->dim) {
		message("size: %d %d %d\n", lm->X->m, lm->y->dim, lm->X->n);
		ErrMsg(ER_IMPOSVAL, "calc_lm: matrices wrong size");
	}
	if (lm->X->m < lm->X->n) {
		lm->is_singular = 1;
		return lm;
	}
/* 
 * allocate structures: 
 */
	lm->is_singular = 0;
	lm->beta = v_resize(lm->beta, lm->X->n);
	lm->Xty = v_resize(lm->Xty, lm->X->n);
	tmp = v_resize(tmp, lm->X->n);
	if (lm->X->n == 0 || lm->y->dim == 0)
		return lm;

	if (DEBUG_COV) {
		printlog("#y is "); v_logoutput(lm->y);
		printlog("#X is "); m_logoutput(lm->X);
		if (lm->weights) {
			printlog("#w is "); v_logoutput(lm->weights);
		}
	}
/*
 * create, in case of weights, V^{-1/2}X and V^{-1/2}y:
 */
	if (lm->weights != VNULL) {
		for (i = 0; i < lm->X->m; i++) { /* rows */
			w = sqrt(lm->weights->ve[i]);
			for (j = 0; j < lm->X->n; j++) /* cols */
				ME(lm->X, i, j) *= w;
			lm->y->ve[i] *= w;
		}
	}
/* 
 * create Chol = X'X (or X'WX) and XtY = (y'X)' = X'y  (X'Wy)
 */
	lm->Xty = vm_mlt(lm->X, lm->y, lm->Xty);
	if (DEBUG_COV) {
		printlog("#X'y is "); v_logoutput(lm->Xty);
	}

	lm->Chol = mtrm_mlt(lm->X, lm->X, lm->Chol);
	if (DEBUG_COV) {
		printlog("#X'X is ");
		m_logoutput(lm->Chol);
	}
	lm->Cov = m_copy(lm->Chol, lm->Cov); /* save copy of X'X */

	int info;
	lm->Chol = CHfactor(lm->Chol, PNULL, &info);
	if (info != 0) {
		lm->is_singular = 1;
		return lm;
	}
	lm->beta = CHsolve1(lm->Chol, lm->Xty, lm->beta, PNULL);
	if (DEBUG_COV) {
		printlog("#beta is "); 
		v_logoutput(lm->beta);
	}
/*
 * estimate error variance:
 */
	tmp = mv_mlt(lm->X, lm->beta, tmp);  /* X b */
	tmp = v_sub(lm->y, tmp, tmp); /* e = y - X b */
	if (lm->weights) {
		for (i = 0, lm->SSErr = 0.0; i < lm->X->m; i++)
			lm->SSErr += lm->weights->ve[i] * SQR(tmp->ve[i]);
	} else
		lm->SSErr = in_prod(tmp, tmp); /* e'e */
	if (DEBUG_COV)
		printlog("#SSErr is %g\n", lm->SSErr);
/*
 * estimate SSReg (Draper & Smith, p. 81, Section 2.2)
 * (unweighted, corrected): beta' X'X beta - n * y_mean^2
 * (weighted, corrected): sum weight_i * (x_i beta - beta_0)^2
 * (unweighted): beta' X'X beta
 * (weighted): see below
 */
	tmp = v_resize(tmp, lm->X->n);
	tmp = vm_mlt(lm->Cov, lm->beta, tmp);
	lm->SSReg = in_prod(lm->beta, tmp);
	if (lm->has_intercept) {
		for (i = 0, y_mean = 0.0; i < lm->y->dim; i++)
			y_mean += lm->y->ve[i];
		y_mean /= lm->y->dim;
		lm->SSReg -= SQR(y_mean) * lm->y->dim;
	}
	lm->dfReg = lm->X->n;
	if (lm->has_intercept)
		lm->dfReg -= 1;
	if (lm->dfReg > 0)
		lm->MSReg = lm->SSReg/lm->dfReg;
	else
		lm->MSReg = DBL_MAX;
 	lm->dfE = lm->X->m - zero_weights_count(lm) - lm->X->n;
	if (lm->dfE == 0)
		lm->MSErr = DBL_MAX;
	else
		lm->MSErr = lm->SSErr/lm->dfE;
	lm->Cov = m_inverse(lm->Cov, &info); /* (X'X)-1 */
	if (info != 0)
		pr_warning("linear model has singular covariance matrix");
	/* next, multiply with sigma^2 */
	sm_mlt(lm->MSErr, lm->Cov, lm->Cov); /* in situ mlt */
	return lm;
}

static int zero_weights_count(LM *lm) {
	int i, n_zero = 0;

	if (lm->weights == VNULL)
		return 0;
	for (i = 0; i < lm->weights->dim; i++)
		if (lm->weights->ve[i] < gl_zero)
			n_zero++;
	return n_zero;
}

double calc_mu(const DATA *d, const DPOINT *where) {
	int i;
	double mu, *from;

	assert(d->beta);

	mu = 0.0;
	from = where->X;
	for (i = 0; i < d->beta->size; i++)
		mu += from[i] * d->beta->val[i];
	return mu;
}

