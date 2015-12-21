#ifndef MTRXH
# define MTRXH
/* interface copied from meschach; implementation rewritten from scratch */
typedef struct {
	size_t m, n, /* #rows, #cols */ 
		max; /* max size, memory allocated */
	double *v;
} MAT; /* dense matrix */
#define ME(X,i,j) X->v[j * X->m + i] /* row i, column j, column-major access */

typedef	struct	{
	size_t dim, max;
	double *ve;
} VEC; /* vector: row or column, whatever matches */

typedef	struct	{
	size_t size, max;
	int *pe;
} PERM;

typedef	struct	{
	size_t size, max;
	int *ive;
} IVEC;

#define PNULL (PERM *) NULL
#define MNULL (MAT *) NULL
#define VNULL (VEC *) NULL
#define IVNULL (IVEC *) NULL

#define M_FREE(x) { if (x != NULL) m_free(x); x = MNULL; }
#define V_FREE(x) { if (x != NULL) v_free(x); x = VNULL; }
#define P_FREE(x) { if (x != NULL) px_free(x); x = PNULL; }
void m_free(MAT *m);
void v_free(VEC *v);
void iv_free(IVEC *v);
void px_free(PERM *p);
#define m_get(i,j) m_resize(MNULL, i, j)
#define v_get(i) v_resize(VNULL, i)

MAT *m_resize(MAT *mat, size_t m, size_t n);
VEC *v_resize(VEC *v, size_t n);
PERM *px_resize(PERM *p, size_t n);
IVEC *iv_resize(IVEC *v, size_t n);
MAT *m_zero(MAT *m);
VEC *v_zero(VEC *v);
MAT *m_inverse(MAT *in, int *info);
VEC *vm_mlt(MAT *m, VEC *v, VEC *out);
VEC *mv_mlt(MAT *m, VEC *v, VEC *out);
MAT *m_mlt(MAT *m1, MAT *m2, MAT *out);
MAT *mtrm_mlt(MAT *m1, MAT *m2, MAT *out);
VEC *v_sub(VEC *v1, VEC *v2, VEC *out);
MAT *m_sub(MAT *m1, MAT *m2, MAT *out);
VEC *v_add(VEC *v1, VEC *v2, VEC *out);
VEC *sv_mlt(double s, VEC *v1, VEC *v2);
MAT *m_add(MAT *m1, MAT *m2, MAT *out);
MAT *m_copy(MAT *in, MAT *out);
VEC *v_copy(VEC *in, VEC *out);
double v_norm2(VEC *v);
MAT *CHsolve(MAT *A, MAT *b, MAT *out, PERM *piv);
VEC *CHsolve1(MAT *A, VEC *b, VEC *out, PERM *piv);
MAT *CHfactor(MAT *A, PERM *piv, int *info);
double in_prod(VEC *a, VEC *b);
MAT *sm_mlt(double s, MAT *m1, MAT *out);
MAT *ms_mltadd(MAT *m1, MAT *m2, double s, MAT *out);
MAT *mmtr_mlt(MAT *m1, MAT *m2, MAT *out);
void m_logoutput(MAT *a);
void v_logoutput(VEC *x);
#endif
