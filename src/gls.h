enum GLS_WHAT {
	GLS_BLUE /* generalized least squares best linear unbiased estimate */,
	GLS_BLUP /* gls best linear unbiased predictor */,
	GLS_BLP  /* gls best linear predictor */,
	UPDATE /* update estimate: use previously calculated weights */,
	GLS_INIT /* initial value */
};

void gls(DATA **d, int n_vars, enum GLS_WHAT pred, DPOINT *where, double *est);
double *make_gls(DATA *d, int calc_residuals);
double *make_gls_mv(DATA **d, int n_vars);
void free_glm(void *v_glm);
