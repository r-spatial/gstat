#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

extern SEXP gstat_init(SEXP s_debug_level);
extern SEXP gstat_load_ev(SEXP np, SEXP dist, SEXP gamma);
extern SEXP gstat_fit_variogram(SEXP fit, SEXP fit_sill, SEXP fit_range);
extern SEXP gstat_exit(SEXP x);
extern SEXP gstat_new_data(SEXP sy, SEXP slocs, SEXP sX, SEXP has_intercept, 
			SEXP beta, SEXP nmax, SEXP nmin, SEXP maxdist, SEXP force,
			SEXP vfn, SEXP sw, SEXP grid, SEXP degree, SEXP is_projected,
			SEXP vdist, SEXP lambda, SEXP omax);
extern SEXP gstat_new_dummy_data(SEXP loc_dim, SEXP has_intercept, SEXP beta, 
		SEXP nmax, SEXP nmin, SEXP maxdist, SEXP vfn, SEXP is_projected,
		SEXP vdist);
extern SEXP gstat_debug_level(SEXP level);
extern SEXP gstat_load_variogram(SEXP s_ids, SEXP s_model, SEXP s_sills, SEXP s_ranges, 
		SEXP s_kappas, SEXP s_anis_all, SEXP s_table, SEXP s_max_val); 
extern SEXP gstat_predict(SEXP sn, SEXP slocs, SEXP sX, SEXP block_cols, SEXP block, 
			SEXP weights, SEXP nsim, SEXP blue);
extern SEXP gstat_set_method(SEXP to);
extern SEXP gstat_set_set(SEXP arg, SEXP val);
extern SEXP gstat_set_merge(SEXP a, SEXP b, SEXP c, SEXP d);
extern SEXP gstat_variogram(SEXP s_ids, SEXP cutoff, SEXP width, SEXP direction, 
		SEXP cressie, SEXP dX, SEXP boundaries, SEXP grid, SEXP cov,
		SEXP pseudo);
extern SEXP gstat_variogram_values(SEXP ids, SEXP pars, SEXP covariance, SEXP dist_values);
extern SEXP gstat_get_variogram_models(SEXP dolong);

#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}

const static R_CallMethodDef R_CallDef[] = {
	CALLDEF(gstat_init, 1),
	CALLDEF(gstat_load_ev, 3),
	CALLDEF(gstat_fit_variogram, 3),
	CALLDEF(gstat_exit, 1),
	CALLDEF(gstat_new_data, 17),
	CALLDEF(gstat_new_dummy_data, 9),
	CALLDEF(gstat_debug_level, 1),
	CALLDEF(gstat_load_variogram, 8),
	CALLDEF(gstat_predict, 8),
	CALLDEF(gstat_set_method, 1),
	CALLDEF(gstat_set_set, 2),
	CALLDEF(gstat_set_merge, 4),
	CALLDEF(gstat_variogram, 10),
	CALLDEF(gstat_variogram_values, 4),
	CALLDEF(gstat_get_variogram_models, 1),
	{NULL, NULL, 0}
};

void
// attribute_visible  // optional
R_init_gstat(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, R_CallDef, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}
