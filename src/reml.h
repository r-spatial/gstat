VARIOGRAM *reml_sills(DATA *d, VARIOGRAM *vp);

#ifdef MATRIXH
MAT *XVXt_mlt(MAT *X, MAT *V, MAT *out);
MAT *XtVX_mlt(MAT *X, MAT *V, MAT *out);
MAT *XdXt_mlt(MAT *X, VEC *d, MAT *out);
MAT *XtdX_mlt(MAT *X, VEC *d, MAT *out);
#endif
