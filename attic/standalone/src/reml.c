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
 * iterated mivque (reml) estimate of covariance components, following
 * Christensen's 1994 derivations Math.Geol. 25 (5), 541-558, and
 * Kitanidis' 1985 derivation (for estimate covariances) Math.Geol. 17 (2).
 *
 * Pre: data pointer (later: extend to 2 or more, crossvariograms?)
 *      initial variogram for *data
 * Post: returns pointer to the updated variogram
 * if output != NULL, all iteration steps are printed to log device,
 * with final variogram model parameter covariances
 * Aug 1994; 
 * Mod. (I-Aw) Nov 1998
 */
#include <stdio.h>
#include <math.h> /* fabs() */

#include "matrix2.h"

#include "defs.h"
#include "userio.h"
#include "debug.h"
#include "data.h"
#include "utils.h"
#include "vario.h"
#include "glvars.h"
#include "select.h"
#include "lm.h" /* get_X, get_y */
#include "reml.h"

static int reml(VEC *Y, MAT *X, MAT **Vk, int n_k, int max_iter,
	double fit_limit, VEC *teta);
static double trace_matrix(MAT *m);
static MAT *calc_VinvIminAw(MAT *Vw, MAT *X, MAT *VinvIminAw, int calc_Aw);
static void calc_rhs_Tr_m(int n_models, MAT **Vk,MAT *VinvIminAw, 
	VEC *Y, VEC *rhs, MAT *Tr_m);
static double calc_ll(MAT *Vw, MAT *X, VEC *y, int n);

static MAT *IminAw = MNULL;

VARIOGRAM *reml_sills(DATA *data, VARIOGRAM *vp) {
	int i, j, k;
	MAT **Vk = NULL, *X = MNULL;
	VEC *Y = VNULL, *init = VNULL;
	DPOINT *dpa, *dpb;
	double dx, dy = 0.0, dz = 0.0, dzero2;

	if (data == NULL  || vp == NULL)
		ErrMsg(ER_NULL, "reml()");
	select_at(data, (DPOINT *) NULL);
	if (vp->n_models <= 0)
		ErrMsg(ER_VARNOTSET, "reml: please define initial variogram model");
/*
 * create Y, X, Vk's only once:
 */
	Y = get_y(&data, Y, 1);
	X = get_X(&data, X, 1);
	Vk = (MAT **) emalloc(vp->n_models * sizeof(MAT *));
	init = v_resize(init, vp->n_models);

	for (i = 0; i < vp->n_models; i++) {
		init->ve[i] = vp->part[i].sill; /* remember init. values for updating */
		vp->part[i].sill = 1;
		Vk[i] = m_resize(MNULL, X->m, X->m);
	}
	dzero2 = gl_zero * gl_zero;
	for (i = 0; i < data->n_list; i++) {
		for (j = 0; j < vp->n_models; j++) /* fill diagonals */
			Vk[j]->me[i][i] = Covariance(vp->part[j], 0.0, 0.0, 0.0);
		for (j = 0; j < i; j++) { /* off-diagonal elements: */
			dpa = data->list[i];
			dpb = data->list[j];
			/* 
		 	 * if different points coincide on a locations, shift them,
		 	 * or the covariance matrix will become singular
		 	 */
			dx = dpa->x - dpb->x;
			dy = dpa->y - dpb->y;
			dz = dpa->z - dpb->z;
			if (data->pp_norm2(dpa, dpb) < dzero2) {
				if (data->mode & X_BIT_SET)
					dx = (dx >= 0 ? gl_zero : -gl_zero);
				if (data->mode & Y_BIT_SET)
					dy = (dy >= 0 ? gl_zero : -gl_zero);
				if (data->mode & Z_BIT_SET)
					dz = (dz >= 0 ? gl_zero : -gl_zero);
			}
			for (k = 0; k < vp->n_models; k++)
				Vk[k]->me[i][j] = Vk[k]->me[j][i] = 
						Covariance(vp->part[k], dx, dy, dz);
		}
	}
	if (reml(Y, X, Vk, vp->n_models, gl_iter, gl_fit_limit, init))
		vp->ev->refit = 0;
	else /* on convergence */
		pr_warning("no convergence while fitting variogram");
	for (i = 0; i < vp->n_models; i++)
		vp->part[i].sill = init->ve[i];
	update_variogram(vp);
	if (DEBUG_VGMFIT)
		logprint_variogram(vp, 1);
	for (i = 0; i < vp->n_models; i++)
		m_free(Vk[i]); 
	efree(Vk);
	m_free(X);
	v_free(Y);
	v_free(init);
	return vp;
}

static int reml(VEC *Y, MAT *X, MAT **Vk, int n_k, int max_iter,
	double fit_limit, VEC *teta) {
 	volatile int n_iter = 0;
 	int i;
	volatile double rel_step = DBL_MAX;
	VEC *rhs = VNULL;
	VEC *dteta = VNULL;
	MAT *Vw = MNULL, *Tr_m = MNULL, *VinvIminAw = MNULL;

	Vw = m_resize(Vw, X->m, X->m);
	VinvIminAw = m_resize(VinvIminAw, X->m, X->m);
	rhs = v_resize(rhs, n_k);
	Tr_m = m_resize(Tr_m, n_k, n_k);
	dteta = v_resize(dteta, n_k);
	while (n_iter < max_iter && rel_step > fit_limit) {
		print_progress(n_iter, max_iter);
		n_iter++;
		dteta = v_copy(teta, dteta);
		/* fill Vw, calc VinvIminAw, rhs; */
		for (i = 0, m_zero(Vw); i < n_k; i++)
			ms_mltadd(Vw, Vk[i], teta->ve[i], Vw); /* Vw = Sum_i teta[i]*V[i] */
		VinvIminAw = calc_VinvIminAw(Vw, X, VinvIminAw, n_iter == 1);
		calc_rhs_Tr_m(n_k, Vk, VinvIminAw, Y, rhs, Tr_m);
		/* Tr_m * teta = Rhs; symmetric, solve for teta: */
		LDLfactor(Tr_m);
		LDLsolve(Tr_m, rhs, teta);
		if (DEBUG_VGMFIT) {
			printlog("teta_%d [", n_iter);
			for (i = 0; i < teta->dim; i++)
				printlog(" %g", teta->ve[i]);
			printlog("] -(log.likelyhood): %g\n",
				calc_ll(Vw, X, Y, n_k));
		}
		v_sub(teta, dteta, dteta); /* dteta = teta_prev - teta_curr */
		if (v_norm2(teta) == 0.0)
			rel_step = 0.0;
		else
			rel_step = v_norm2(dteta) / v_norm2(teta);
	} /* while (n_iter < gl_iter && rel_step > fit_limit) */

	print_progress(max_iter, max_iter);
	if (n_iter == gl_iter)
		pr_warning("No convergence after %d iterations", n_iter);

	if (DEBUG_VGMFIT) { /* calculate and report covariance matrix */
		/* first, update to current est */
		for (i = 0, m_zero(Vw); i < n_k; i++)
			ms_mltadd(Vw, Vk[i], teta->ve[i], Vw); /* Vw = Sum_i teta[i]*V[i] */
		VinvIminAw = calc_VinvIminAw(Vw, X, VinvIminAw, 0);
		calc_rhs_Tr_m(n_k, Vk, VinvIminAw, Y, rhs, Tr_m);
		m_inverse(Tr_m, Tr_m);
		sm_mlt(2.0, Tr_m, Tr_m); /* Var(YAY)=2tr(AVAV) */
		printlog("Lower bound of parameter covariance matrix:\n");
		m_logoutput(Tr_m);
		printlog("# Negative log-likelyhood: %g\n", calc_ll(Vw, X, Y, n_k));
	}
	m_free(Vw);
	m_free(VinvIminAw);
	m_free(Tr_m);
	v_free(rhs);
	v_free(dteta);
	return (n_iter < max_iter && rel_step < fit_limit); /* converged? */
}

static MAT *calc_VinvIminAw(MAT *Vw, MAT *X, MAT *VinvIminAw, int calc_Aw) {
/*
 * calculate V_w^-1(I-A_w) (==VinvIminAw),
 * A = X(X'X)^-1 X' (AY = XBeta; Beta = (X'X)^-1 X'Y)
 *
 * on second thought (Nov 1998 -- more than 4 years later :-))
 * calc (I-Aw) only once and keep this constant during iteration.
 */
 	MAT *tmp = MNULL, *V = MNULL;
 	VEC *b = VNULL, *rhs = VNULL;
 	int i, j;

	if (X->m != Vw->n || VinvIminAw->m != X->m)
		ErrMsg(ER_IMPOSVAL, "calc_VinvIminAw: sizes don't match");
	
	if (calc_Aw) {
		IminAw = m_resize(IminAw, X->m, X->m);
		tmp = m_resize(tmp, X->n, X->n);
		tmp = mtrm_mlt(X, X, tmp); /* X'X */
		m_inverse(tmp, tmp); /* (X'X)-1 */
		/* X(X'X)-1 -> X(X'X)-1 X') */
		IminAw = XVXt_mlt(X, tmp, IminAw);
		for (i = 0; i < IminAw->m; i++) /* I - Aw */
			for (j = 0; j <= i; j++)
				if (i == j)
					IminAw->me[i][j] = 1.0 - IminAw->me[i][j];
				else
					IminAw->me[i][j] = IminAw->me[j][i] = -IminAw->me[i][j];
	}

	V = m_copy(Vw, V);
	LDLfactor(V);

	rhs = v_resize(rhs, X->m);
	b = v_resize(b, X->m);

	for (i = 0; i < X->m; i++) { /* solve Vw X = (I-A) for X -> V-1(I-A) */
		rhs = get_col(IminAw, i, rhs);
		LDLsolve(V, rhs, b);
		set_col(VinvIminAw, i, b);
	}
	v_free(rhs);
	v_free(b);
	m_free(V);

	if (tmp) 
		m_free(tmp);

	return VinvIminAw;
}

static void calc_rhs_Tr_m(int n_models, MAT **Vk,MAT *VinvIminAw, 
		VEC *y, VEC *rhs, MAT *Tr_m) {
	int j, k;
	MAT **Pr = NULL, *Tmp = MNULL;
	VEC *v_tmp = VNULL, *v_tmp2;

	Pr = (MAT **) emalloc(n_models * sizeof(MAT *));
	v_tmp2 = vm_mlt(VinvIminAw, y, VNULL); /* Vw-(I-Aw)Y == Y'(I-Aw)'Vw- */
	for (j = 0; j < n_models; j++) {
		Pr[j] = m_mlt(Vk[j], VinvIminAw, MNULL);
		Tmp = m_mlt(Pr[j], Pr[j], Tmp);
		Tr_m->me[j][j] = trace_matrix(Tmp); /* diagonal */
		/* using Tr(A B) == Tr(B A) */
		for (k = 0; k < j; k++) { /* we did Pr[k] and Pr[j], so */
			Tmp = m_mlt(Pr[j], Pr[k], Tmp); /* off-diagonal */
			Tr_m->me[j][k] = Tr_m->me[k][j] = trace_matrix(Tmp);
		}
		v_tmp = vm_mlt(Vk[j], v_tmp2, v_tmp); /* Vw-1(I-Aw)Y */
		rhs->ve[j] = in_prod(v_tmp2, v_tmp);
	}
	for (j = 0; j < n_models; j++)
		m_free(Pr[j]);
	efree(Pr);
	m_free(Tmp);
	v_free(v_tmp);
	v_free(v_tmp2);
	return;
}

static double calc_ll(MAT *Vw, MAT *X, VEC *y, int n) {
/*
 * calculate negative log-likelyhood
 */
	static MAT *M1 = MNULL;
	static VEC *res = VNULL, *tmp = VNULL;
	double zQz;
	volatile double ldet;
	int i;

	IminAw->m -= n;

	/* |B'(I-A)Vw(I-A)'B|, pretty inefficiently, can 4 x as fast: */
	/* M1 = m_mlt(IminAw, Vw, M1); M2 = mmtr_mlt(M1, IminAw, M2); */
	M1 = XVXt_mlt(IminAw, Vw, M1);
	LDLfactor(M1);
	for (i = 0, ldet = 0.0; i < M1->m; i++) {
		assert(M1->me[i][i] > 0.0);
		ldet += log(M1->me[i][i]);
	}
	/* y'B'A'(B'A'Vw A B)-1 A B y */
	res = mv_mlt(IminAw, y, res); /* the m-n residuals B(I-A)'Y */
	tmp = LDLsolve(M1, res, tmp);
		/* M1 tmp = res -> tmp = M1-1 res */
	zQz = in_prod(res, tmp);  /* res' M1inv res */

	IminAw->m += n;

	return 0.5 * ((Vw->m - n)*log(2*PI) + ldet + zQz);
}

static double trace_matrix(MAT *m) {
/* returns trace of a square matrix */
	int i;
	double trace;

	if (m == NULL)
		ErrMsg(ER_NULL, "trace_matrix: NULL argument");
	if (m->m != m->n)
		ErrMsg(ER_IMPOSVAL, "trace_matrix: non-square matrix");
	for (i = 0, trace = 0.0; i < m->m; i++)
		trace += m->me[i][i];
	return trace;
}

MAT *XtVX_mlt(MAT *X, MAT *V, MAT *out) {
/* for a symmetric matrix V, return X' V X */
	static MAT *VX = MNULL;
	int i, j, k;

	if (X==(MAT *)NULL || V==(MAT *)NULL )
		error(E_NULL, "XtVX_mlt");
	if (X->m != V->m)
		error(E_SIZES, "XtVX_mlt");
	if (V->m != V->n)
		error(E_SQUARE, "XtVX_mlt");

	out = m_resize(out, X->n, X->n);
	VX = m_resize(VX, V->m, X->n);
	m_zero(out);

	VX = m_mlt(V, X, VX);
	for (i = 0; i < X->n; i++) {
		for (j = i; j < X->n; j++)
			for (k = 0; k < X->m; k++)
				out->me[i][j] += X->me[k][i] * VX->me[k][j];
		for (j = 0; j <= i; j++) /* symmetry */
			out->me[i][j] = out->me[j][i];
	}
	return out;
}

MAT *XVXt_mlt(MAT *X, MAT *V, MAT *out) {
/* for a symmetric matrix V, return X V X' */
	static MAT *VXt = MNULL;
	int i, j, k;

	if (X==(MAT *)NULL || V==(MAT *)NULL )
		error(E_NULL, "XtVX_mlt");
	if (X->n != V->m)
		error(E_SIZES, "XtVX_mlt");
	if (V->m != V->n)
		error(E_SQUARE, "XtVX_mlt");

	out = m_resize(out, X->m, X->m);
	VXt = m_resize(VXt, V->m, X->n);
	m_zero(out);

	VXt = mmtr_mlt(V, X, VXt);
	for (i = 0; i < X->m; i++) {
		for (j = i; j < X->m; j++)
			for (k = 0; k < X->n; k++)
				out->me[i][j] += X->me[i][k] * VXt->me[k][j];
		for (j = 0; j <= i; j++) /* symmetry */
			out->me[i][j] = out->me[j][i];
	}
	return out;
}

MAT *XdXt_mlt(MAT *X, VEC *d, MAT *out) {
/* for a diagonal matrix in d, return X d X' */
	int i, j, k;

	if (X==(MAT *)NULL || d==(VEC *)NULL )
		error(E_NULL, "XVXt_mlt");
	if (X->n != d->dim)
		error(E_SIZES, "XVXt_mlt");

	out = m_resize(out, X->n, X->n);
	m_zero(out);

	for (i = 0; i < X->m; i++) {
		for (j = i; j < X->m; j++)
			for (k = 0; k < X->n; k++)
				out->me[i][j] += X->me[i][k] * X->me[j][k] * d->ve[k];
		for (j = 0; j <= i; j++) /* symmetry */
			out->me[i][j] = out->me[j][i];
	}
	return out;
}

MAT *XtdX_mlt(MAT *X, VEC *d, MAT *out) {
/* for a diagonal matrix in d, return X' d X */
	int i, j, k;

	if (X==(MAT *)NULL || d==(VEC *)NULL )
		error(E_NULL, "XtVX_mlt");
	if (X->m != d->dim)
		error(E_SIZES, "XtVX_mlt");

	out = m_resize(out, X->n, X->n);
	m_zero(out);

	for (i = 0; i < X->n; i++) {
		for (j = i; j < X->n; j++)
			for (k = 0; k < X->m; k++)
				out->me[i][j] += X->me[k][i] * X->me[k][j] * d->ve[k];
		for (j = 0; j <= i; j++) /* symmetry */
			out->me[i][j] = out->me[j][i];
	}
	return out;
}
