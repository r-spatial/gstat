/* interface roughly follows meschach; implementation rewritten from scratch */
#include <string.h> /* memcpy, memset */
#include <math.h> /* fabs */

#define USE_FC_LEN_T
#include <R_ext/Lapack.h>
#ifndef FCONE
# define FCONE
#endif

#include <R.h>

#include "defs.h" /* CDECL */
#include "utils.h" /* efree, emalloc */
#include "userio.h" /* ErrMsg */
#include "glvars.h" /* gl_blas */
#include "debug.h"
#include "mtrx.h"

/* get rid of -0.000000 output: */
#define _zero_(x) (fabs(x) < 1.e-7 ? 0.0 : x)

/* 0. book keeping: initialisation, memory allocation, zero, copy, print */

void m_free(MAT *m) {
	efree(m->v);
	efree(m);
}

void v_free(VEC *v) {
	efree(v->ve);
	efree(v);
}

void iv_free(IVEC *iv) {
	efree(iv->ive);
	efree(iv);
}

void px_free(PERM *p) {
	efree(p->pe);
	efree(p);
}

MAT *m_init(void) {
	MAT *mat = emalloc(sizeof(MAT));
	mat->n = mat->m = mat->max = 0;
	mat->v = (double *) NULL;
	return(mat);
}

MAT *m_resize(MAT *mat, size_t m, size_t n) {
	if (mat == MNULL)
		mat = m_init();
	if (m * n > mat->max) {
		mat->max = m * n;
		mat->v = (double *) erealloc(mat->v, mat->max * sizeof(double)); /* takes care of NULL m */
	}
	mat->m = m;
	mat->n = n;
	return(mat);
}

VEC *v_init(void) {
	VEC *v = emalloc(sizeof(VEC));
	v->dim = v->max = 0;
	v->ve = NULL;
	return(v);
}

VEC *v_resize(VEC *v, size_t n) {
	if (v == NULL)
		v = v_init();
	if (n > v->max) {
		v->ve = erealloc(v->ve, n * sizeof(double));
		v->max = n;
	}
	v->dim = n;
	return(v);
}

PERM *p_init(void) {
	PERM *p = emalloc(sizeof(PERM));
	p->size = p->max = 0;
	p->pe = (int *) NULL;
	return(p);
}

PERM *px_resize(PERM *p, size_t n) {
	if (p == PNULL)
		p = p_init();
	if (n > p->max) {
		p->pe = erealloc(p->pe, n * sizeof(size_t));
		p->max = n;
	}
	p->size = n;
	return(p);
}

IVEC *iv_init(void) {
	IVEC *iv = emalloc(sizeof(IVEC));
	iv->size = iv->max = 0;
	iv->ive = (int *) NULL;
	return(iv);
}

IVEC *iv_resize(IVEC *iv, size_t n) {
	if (iv == IVNULL)
		iv = iv_init();
	if (n > iv->max) {
		iv->ive = erealloc(iv->ive, n * sizeof(int));
		iv->max = n;
	}
	iv->size = n;
	return(iv);
}

MAT *m_zero(MAT *m) {
	if (m != MNULL)
		memset(m->v, 0x00, m->m * m->n * sizeof(double));
	return(m);
}

VEC *v_zero(VEC *v) {
	if (v != VNULL)
		memset(v->ve, 0x00, v->dim * sizeof(double));
	return(v);
}

MAT *m_copy(MAT *in, MAT *out) {
	if (in == out)
		return(out);
	out = m_resize(out, in->m, in->n);
	memcpy(out->v, in->v, in->m * in->n * sizeof(double));
	return(out);
}

VEC *v_copy(VEC *in, VEC *out) {
	if (in == out)
		return(out);
	out = v_resize(out, in->dim);
	memcpy(out->ve, in->ve, in->dim * sizeof(double));
	return(out);
}

void m_logoutput(MAT * a) {
	unsigned int i, j, tmp;

	if (a == (MAT *) NULL) {
		printlog("#Matrix: NULL\n");
		return;
	}
	printlog("#Matrix: %d by %d\n", a->m, a->n);
	if (a->v == NULL) {
		printlog("NULL\n");
		return;
	}
	printlog("rbind(\n");
	for (i = 0; i < a->m; i++) {	/* for each row... */
		printlog("c(");
		for (j = 0, tmp = 2; j < a->n; j++, tmp++) {
			/* for each col in row: */
			printlog("%9f", _zero_(ME(a, i, j)));
			if (j + 1 < a->n)
				printlog(", ");
			else 
				printlog(")");
		}
		if (i + 1 < a->m)
			printlog(", ");
		else
			printlog("  ");
		printlog("# row %u\n", i + 1);
	}
	printlog(")\n");
}

void v_logoutput(VEC * x) {
	unsigned int i, tmp;

	if (x == (VEC *) NULL) {
		printlog("#Vector: NULL\n");
		return;
	}
	printlog("#Vector: dim: %d\n", x->dim);
	if (x->ve == NULL) {
		printlog("NULL\n");
		return;
	}
	printlog("c(");
	for (i = 0, tmp = 0; i < x->dim; i++, tmp++) {
		printlog("%9f", _zero_(x->ve[i]));
		if (i + 1 < x->dim)
			printlog(", ");
	}
	printlog(")");
}

/* 1: vector-scalar, vector-vector (BLAS-1) */
VEC *sv_mlt(double s, VEC *v, VEC *out) { /* out <- s * v */
	out = v_resize(out, v->dim);
	for (int i = 0; i < v->dim; i++)
		out->ve[i] = s * v->ve[i];
	return(out);
}

double v_norm2(VEC *v) { /* 2-norm  */
	return(in_prod(v, v));
}

VEC *v_add(VEC *v1, VEC *v2, VEC *out) { /* out = v1 + v2 */
	if (v1->dim != v2->dim)
		ErrMsg(ER_IMPOSVAL, "v_sub size mismatch");
	out = v_resize(out, v1->dim); 
	for (int i = 0; i < out->dim; i++)
		out->ve[i] = v1->ve[i] + v2->ve[i];
	return(out);
}

VEC *v_sub(VEC *v1, VEC *v2, VEC *out) { /* out = v1 - v2 = -1 * v2 + v1 */
	if (v1->dim != v2->dim)
		ErrMsg(ER_IMPOSVAL, "v_sub size mismatch");
	out = v_resize(out, v1->dim);
	for (int i = 0; i < out->dim; i++)
		out->ve[i] = v1->ve[i] - v2->ve[i];
	return(out);
}

double in_prod(VEC *a, VEC *b) { /* a'b */
	if (a->dim != b->dim)
		ErrMsg(ER_IMPOSVAL, "in_prod: dimensions don't match");

	if (! gl_blas) {
		double d = 0.0;
		for (int i = 0; i < a->dim; i++)
			d += a->ve[i] * b->ve[i];
		return(d);
	} else {
		int one = 1;
		return(F77_CALL(ddot)((int *) &(a->dim), a->ve, &one, b->ve, &one));
	}
}

/* 2: vector-matrix (BLAS-2) */
VEC *vm_mlt(MAT *m, VEC *v, VEC *out) { /* out <- v m */
	if (m->m != v->dim)
		ErrMsg(ER_IMPOSVAL, "vm_mlt: dimensions");
	out = v_zero(v_resize(out, m->n));
	if (! gl_blas) {
    	for (size_t i = 0; i < m->n; i++)
        	for (size_t j = 0; j < v->dim; j++)
				out->ve[i] += v->ve[j] * ME(m, j, i);
	} else {
		double alpha = 1.0, beta = 0.0;
		int one = 1;
		F77_CALL(dgemv)("T", (int *) &(m->m), (int *) &(m->n), &alpha, m->v, 
			(int *) &(m->m), v->ve, &one, &beta, out->ve, &one FCONE);
	}
	return(out);
}

VEC *mv_mlt(MAT *m, VEC *v, VEC *out) { /* out <- m v */
	if (v == out)
		ErrMsg(ER_IMPOSVAL, "mv_mlt in situ");
	if (m->n != v->dim)
		ErrMsg(ER_IMPOSVAL, "mv_mlt non-matching sizes");
	out = v_zero(v_resize(out, m->m));
	if (! gl_blas) {
		for (int j = 0; j < m->m; j++)
			for (int i = 0; i < m->n; i++)
				out->ve[j] += ME(m, j, i) * v->ve[i];
	} else {
		double alpha = 1.0, beta = 0.0;
		int one = 1;
		F77_CALL(dgemv)("N", (int *) &(m->m), (int *) &(m->n), &alpha, m->v, 
			(int *) &(m->m), v->ve, &one, &beta, out->ve, &one FCONE);
	}
	return(out);
}

/* 3: matrix-matrix (BLAS-3) */
MAT *m_mlt(MAT *m1, MAT *m2, MAT *out) { /* out <- m1 %*% m2 */
	if (m1->n != m2->m)
		ErrMsg(ER_IMPOSVAL, "mv_mlt non-matching sizes");

	if (! gl_blas) {
		out = m_zero(m_resize(out, m1->m, m2->n));
		for (int i = 0; i < m1->m; i++)
			for (int j = 0; j < m2->n; j++)
				for (int k = 0; k < m1->n; k++)
					ME(out, i, j) += ME(m1, i, k) * ME(m2, k, j);
	} else {
		double alpha = 1.0, beta = 0.0;
		out = m_resize(out, m1->m, m2->n);
		F77_CALL(dgemm)("N", "N", (int *) &(m1->m), (int *) &(m2->n), (int *) &(m1->n), &alpha, 
			m1->v, (int *)&(m1->m), 
			m2->v, (int *)&(m2->m), 
			&beta, out->v, (int *) &(m1->m) FCONE FCONE);
	}
	return(out);
}

MAT *mtrm_mlt(MAT *m1, MAT *m2, MAT *out) { /* out <- t(m1) %*% m2 */
	if (m1->m != m2->m)
		ErrMsg(ER_IMPOSVAL, "mtrm_mlt non-matching m arrays");

	out = m_zero(m_resize(out, m1->n, m2->n));
	if (! gl_blas) {
		for (int i = 0; i < m1->n; i++)
			for (int j = 0; j < m2->n; j++)
				for (int k = 0; k < m1->m; k++)
					ME(out, i, j) += ME(m1, k, i) * ME(m2, k, j);
	} else {
		double alpha = 1.0, beta = 0.0;
		F77_CALL(dgemm)("T", "N", (int *) &(m1->n), (int *) &(m2->n), (int *) &(m1->m), &alpha, 
			m1->v, (int *)&(m1->m), 
			m2->v, (int *)&(m2->m), 
			&beta, out->v, (int *) &(m1->n) FCONE FCONE);
	}
	return(out);
}

MAT *mmtr_mlt(MAT *m1, MAT *m2, MAT *out) { /* out <- m1 m2' */
	if (m1->n != m2->n)
		ErrMsg(ER_IMPOSVAL, "mmtr_mlt non-matching m arrays");
	out = m_zero(m_resize(out, m1->m, m2->m));
	if (! gl_blas) {
		for (int i = 0; i < m1->m; i++)
			for (int j = 0; j < m2->m; j++)
				for (int k = 0; k < m1->n; k++)
					ME(out, i, j) += ME(m1, i, k) * ME(m2, j, k);
	} else {
		double alpha = 1.0, beta = 0.0;
		F77_CALL(dgemm)("N", "T", (int *) &(m1->m), (int *) &(m2->m), (int *) &(m1->n), &alpha, 
			m1->v, (int *)&(m1->m), 
			m2->v, (int *)&(m2->m), 
			&beta, out->v, (int *) &(m1->m) FCONE FCONE);
	}
	return(out);
}

MAT *ms_mltadd(MAT *m1, MAT *m2, double s, MAT *out) { /* out <- m1 + s * m2 */
	/* return m1 + s * m2 */
	if (m1->m != m2->m || m1->n != m2->n)
		ErrMsg(ER_IMPOSVAL, "ms_mltadd: dimension mismatch");
	out = m_resize(out, m1->m, m1->n);
	for (int j = 0; j < m1->n; j++)
		for (int i = 0; i < m1->m; i++)
			ME(out, i, j) = ME(m1, i, j) + s * ME(m2, i, j);
	return(out);
}

MAT *sm_mlt(double s, MAT *m1, MAT *out) { /* out <- s * m1 */
	out = m_resize(out, m1->m, m1->n);
	for (int j = 0; j < m1->n; j++)
		for (int i = 0; i < m1->m; i++)
			ME(out, i, j) = s * ME(m1, i, j);
	return(out);
}

MAT *m_add(MAT *m1, MAT *m2, MAT *out) { /* out <- m1 + m2 */
	if (m1->m != m2->m || m1->n != m2->n)
		ErrMsg(ER_IMPOSVAL, "m_add size mismatch");
	out = m_resize(out, m1->m, m1->n);
	for (int j = 0; j < m1->n; j++)
		for (int i = 0; i < m1->m; i++)
			ME(out, i, j) = ME(m1, i, j) + ME(m2, i, j);
	return(out);
}

MAT *m_sub(MAT *m1, MAT *m2, MAT *out) { /* out <- m1 - m2 */
	if (m1->m != m2->m || m1->n != m2->n)
		ErrMsg(ER_IMPOSVAL, "m_sub size mismatch");
	out = m_resize(out, m1->m, m1->n);
	for (int j = 0; j < m1->n; j++)
		for (int i = 0; i < m1->m; i++)
			ME(out, i, j) = ME(m1, i, j) - ME(m2, i, j);
	return(out);
}

/* 4: matrix factorisation, solving systems of equations */
MAT *CHfactor(MAT *m, PERM *piv, int *info) {
    if (m->m != m->n) 
		Rf_error("CHfactor: 'm' must be a square matrix");

	for (int i = 1; i < m->m; i++)
		for (int j = 0; j < i; j++)
			ME(m, i, j) = 0.0; /* zero lower triangle of Fortran order */

	if (piv == PNULL) { /* Choleski: */
		F77_CALL(dpotrf)("Upper", (int *)&(m->n), m->v, (int *)&(m->n), info, (FC_LEN_T) 5);
		if (*info != 0) {
	    	if (*info > 0 && DEBUG_COV)
				Rf_warning("the leading minor of order %d is not positive definite", *info);
	    	if (*info < 0)
				Rf_error("argument %d of Lapack routine %s had invalid value", -(*info), "dpotrf");
		}
	} else { /* LDL': */
    	if (piv->size != m->n) 
			Rf_error("CHfactor: 'piv' must have dimension equal to m->n");
		double w, *work;
		/* first query for size of work, then allocate work, then factorize m: */
		int lwork = -1;
		F77_CALL(dsytrf)("Upper", (int *)&(m->n), m->v, (int *)&(m->n), (int *) piv->pe, &w, &lwork, info, (FC_LEN_T) 5);
		lwork = (int) w;
		work = emalloc(lwork * sizeof(double));
		F77_CALL(dsytrf)("Upper", (int *)&(m->n), m->v, (int *)&(m->n), (int *) piv->pe, work, &lwork, info, (FC_LEN_T) 5);
		efree(work);
		if (*info != 0) {
	    	if (*info > 0 && DEBUG_COV)
				Rf_warning("D[%d,%d] is exactly zero, meaning that D is singular", *info, *info);
	    	if (*info < 0)
				Rf_error("argument %d of Lapack routine %s had invalid value", -(*info), "dsytrf");
		}
	}
	return(m);
}

MAT *CHsolve(MAT *m, MAT *b, MAT *out, PERM *piv) { /* solve A X = B after factorizing A */
	int info;
	if (m->m != m->n) 
		Rf_error("CHsolve: 'm' must be a square matrix");
	if (m->m != b->m) 
		Rf_error("CHsolve: b does not match m");
	out = m_copy(b, out); /* column-major */
	if (piv == PNULL) /* Choleski */
		F77_CALL(dpotrs)("Upper", (int *) &(m->m), (int *) &(b->n), m->v, (int *) &(m->m),          out->v, (int *) &(m->m), &info, (FC_LEN_T) 5);
	else /* LDL' */
		F77_CALL(dsytrs)("Upper", (int *) &(m->m), (int *) &(b->n), m->v, (int *) &(m->m), piv->pe, out->v, (int *) &(m->m), &info, (FC_LEN_T) 5);
	if (info < 0)
		Rf_error("CHsolve: argument %d of Lapack routine %s had invalid value", -info, piv == NULL ? "dpotrs" : "dsytrs");
	return(out);
}

VEC *CHsolve1(MAT *m, VEC *b, VEC *out, PERM *piv) { /* solve A x = b after factorizing A */
	int one = 1, info;
	if (m->m != m->n) 
		Rf_error("CHsolve1: 'm' must be a square matrix");
	if (m->m != b->dim) 
		Rf_error("CHsolve1: vector b does not match m");
	out = v_copy(b, out);
	if (piv == PNULL) 
		F77_CALL(dpotrs)("U", (int *) &(m->m), (int *) &one, m->v, (int *) &(m->m),          out->ve, (int *) &(m->m), &info FCONE);
	else
		F77_CALL(dsytrs)("L", (int *) &(m->m), (int *) &one, m->v, (int *) &(m->m), piv->pe, out->ve, (int *) &(m->m), &info FCONE);
	if (info < 0)
		Rf_error("CHsolve1: argument %d of Lapack routine %s had invalid value", -info, piv == NULL ? "dpotrs" : "dsytrs");
	return(out);
}

MAT *m_inverse(MAT *in, int *info) { /* out <- in^{-1} */
	PERM *piv = px_resize(PNULL, in->m);
	in = CHfactor(in, piv, info);
	if (*info != 0) { /* singular */
		px_free(piv);
		return(in);
	}
	MAT *rhs = m_zero(m_resize(MNULL, in->m, in->m));
	for (int i = 0; i < rhs->m; i++)
		ME(rhs, i, i) = 1.0;
	rhs = CHsolve(in, rhs, rhs, piv);
	in = m_copy(rhs, in);
	m_free(rhs);
	px_free(piv);
	return(in);
}
