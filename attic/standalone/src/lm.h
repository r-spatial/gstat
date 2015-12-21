#ifndef LM_H
# define LM_H
void pred_lm(DATA **data, int n_vars, DPOINT *where, double *est);
void make_residuals_lm(DATA *d);
double *make_ols(DATA *d);

# ifdef MATRIXH /* MAT, VEC definitions */
MAT *get_X(DATA **d, MAT *X, int nvars);
MAT *get_X0(DATA **d, MAT *X0, DPOINT *where, int nvars);
double calc_mu(const DATA *d, const DPOINT *pt);
VEC *get_y(DATA **d, VEC *y, int nvars);
int is_singular(MAT *X, double epsilon);
void m_logoutput(MAT *a);
void v_logoutput(VEC *x);
 
typedef struct {
	VEC *beta, /* parameter vector */
		*y, /* data vector */
		*Xty, /* X'y */
		*weights; /* weights in a WLS model: V-1, 1/sigma^2_i */
	MAT *X, /* design matrix */
		*Cov, /* covariance matrix of beta */
		*Chol; /* Choleski decomposition of X'X or X'V-1X */
	double MSErr, /* Mean Square Error */
		MSReg, /* Mean Square due to regression */
		SSErr, /* Sum of Squares error */
		SSReg, /* Sum of Squares regression */
		cn_max; /* max. allowed condition number; < 0 => don't check */
	int dfE, /* degrees of freedom error */
		dfReg, /* degrees of freedom regression */
		is_singular, /* flag if X'X is singular */
		has_intercept; /* model has intercept, J is part of X */
} LM ;

LM *calc_lm(LM *lm);
void logprint_lm(DATA *d, LM *lm);

LM *init_lm(LM *lm);
void free_lm(LM *lm);

# endif /* MATRIXH */

#endif /* LM_H */
