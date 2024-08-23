/*
 * all functions exposed to R
 */
#include <R.h>
#include <Rinternals.h>
/* #include <Rdefines.h> */

#include "defs.h"
#include "data.h"
#include "select.h"
#include "utils.h"
#include "userio.h"
#include "vario.h"
#include "fit.h"
#include "sem.h"
#include "glvars.h"
#include "debug.h"
#include "mapio.h"
#include "msim.h"
#include "getest.h"
#include "s.h"

static DATA_GRIDMAP *gstat_S_fillgrid(SEXP gridparams);
static void gstat_set_block(long i, SEXP block, SEXP block_cols, DPOINT *current);
int do_print_progress = 0;
#define NAME_SIZE 20 /* buffer size for name */

extern unsigned int n_pred_locs; /* msim.c */

SEXP gstat_init(SEXP s_debug_level) {

	do_print_progress = 0;
	remove_all();
	init_global_variables();
	init_data_minmax();
	GetRNGstate();
	debug_level = INTEGER(s_debug_level)[0];
	if (debug_level < 0) {
		debug_level = -debug_level;
		do_print_progress = 1;
	}
	return(s_debug_level);
}

SEXP gstat_exit(SEXP x) {

	PutRNGstate(); /* write seed back to R/S engine */
	remove_all();
	return(x);
}

SEXP gstat_new_data(SEXP sy, SEXP slocs, SEXP sX, SEXP has_intercept, 
			SEXP beta, SEXP nmax, SEXP nmin, SEXP maxdist, SEXP force,
			SEXP vfn, SEXP sw, SEXP grid, SEXP degree, SEXP is_projected,
			SEXP vdist, SEXP lambda, SEXP omax) {
	double *y, *locs, *X, *w = NULL;
	long i, j, id, n, dim, n_X, has_int;
	DPOINT current;
	DATA **d;
	char name[NAME_SIZE];

	PROTECT(sy = Rf_coerceVector(sy, REALSXP));
	n = LENGTH(sy);
	y = REAL(sy);
	if (n == 0)
		ErrMsg(ER_IMPOSVAL, "no data read");

	if (LENGTH(slocs) % n != 0)
		Rf_error("dimensions do not match: locations %d and data %ld",
			(int) LENGTH(slocs), n);
	dim = LENGTH(slocs) / n;
	if (dim <= 0)
		Rf_error("too few spatial dimensions: %ld", dim);
	if (dim > 3)
		Rf_error("too many spatial dimensions: %ld", dim);
	locs = REAL(slocs);

	if (LENGTH(sw) == n)
		w = REAL(sw);

	if (LENGTH(sX) % n != 0)
		Rf_error("dimensions do not match: X %d and data %ld: missing values in data?",
			(int) LENGTH(sX), n);
	n_X = LENGTH(sX) / n;
	X = REAL(sX);

	assert(n_X > 0);
	current.z = 0.0;
	current.bitfield = 0;

	id = get_n_vars();
	snprintf(name, NAME_SIZE, "var%ld", id);
	which_identifier(name);
	d = get_gstat_data();
	d[id]->id = id;

	d[id]->n_list = d[id]->n_max = 0;
	d[id]->colnx = d[id]->colny = d[id]->colnvalue = d[id]->colnz = 0;
	d[id]->x_coord = "x";
	d[id]->y_coord = "y";
	d[id]->z_coord = "z";
	d[id]->variable = "R data";
	d[id]->fname = "R data";
	d[id]->lambda = REAL(lambda)[0];
	has_int = INTEGER(has_intercept)[0];
	/* increase d[id]->n_X and set d[id]->colX[i]: */
	for (i = d[id]->n_X = 0; i < n_X; i++) 
		data_add_X(d[id], i + (has_int ? 0 : 1)); 
	assert(d[id]->n_X == n_X);
	for (i = 0; i < LENGTH(beta); i++) /* do nothing if beta is numeric(0) */
		d[id]->beta = push_d_vector(REAL(beta)[i], d[id]->beta);
	if (INTEGER(nmax)[0] > 0) /* leave default (large) if < 0 */
		d[id]->sel_max = INTEGER(nmax)[0];
	if (INTEGER(omax)[0] > 0) /* leave default (0) if <= 0 */
		d[id]->oct_max = INTEGER(omax)[0];
	if (INTEGER(nmin)[0] > 0) /* leave default (0) if <= 0 */
		d[id]->sel_min = INTEGER(nmin)[0];
	if (REAL(maxdist)[0] > 0.0)
		d[id]->sel_rad = REAL(maxdist)[0];
	if (INTEGER(force)[0] > 0)
		d[id]->force = INTEGER(force)[0];
	switch(INTEGER(vfn)[0]) {
		case 1: /* d[id]->variance_fn = v_identity; == leave NULL */ break;
		case 2: d[id]->variance_fn = v_mu; break;
		case 3: d[id]->variance_fn = v_bin; break;
		case 4: d[id]->variance_fn = v_mu2; break;
		case 5: d[id]->variance_fn = v_mu3; break;
		default: Rf_error("unknown variance function %d", INTEGER(vfn)[0]);
	}
	gl_longlat = (INTEGER(is_projected)[0] == 0);
	d[id]->mode = X_BIT_SET | V_BIT_SET;
	if (dim > 1)
		d[id]->mode = d[id]->mode | Y_BIT_SET;
	if (dim > 2)
		d[id]->mode = d[id]->mode | Z_BIT_SET;
	set_norm_fns(d[id]); /* so gstat can calculate distances */
	if (w != NULL)
		d[id]->colnvariance = 1; /* it is set */
	switch(LENGTH(grid)) {
		case 0: case 1: break; /* empty, i.e., numeric(0) */
		case 6: d[id]->grid = gstat_S_fillgrid(grid); break;
		default: Rf_error("length of grid topology %d unrecognized", (int) LENGTH(grid));
	}
	d[id]->polynomial_degree = INTEGER(degree)[0];
	if (d[id]->polynomial_degree < 0 || d[id]->polynomial_degree > 3) {
		Rf_error("polynomial degree should be 0, 1, 2 or 3");
	}
	if (d[id]->polynomial_degree > 0) { 
		/* we're doing polynomials through degree: */
		if (id > 0) {
			Rf_error("polynomial degree will only work for a single variable");
		} if (n_X > 1) {
			Rf_error("polynomial degree only works when no other predictors are given");
		}
		setup_polynomial_X(d[id]); /* standardized coordinate polynomials */
	}
	d[id]->vdist = INTEGER(vdist)[0];
	assert(n_X <= d[id]->n_X);
	current.X = (double *) emalloc(d[id]->n_X * sizeof(double));

	SET_POINT(&current);
	current.u.stratum = 0;
	current.attr = current.x = current.y = current.z = 0.0;
	for (i = 0; i < n; i++) { /* loop over points */
		current.attr = y[i];
		current.x = locs[i];
		if (dim >= 2)
			current.y = locs[n + i];
		if (dim >= 3)
			current.z = locs[2 * n + i];
		/* track min/max coordinates, also for z, for the qtree bbox */
		if (i == 0) {
			d[id]->maxX = d[id]->minX = current.x;
			d[id]->maxY = d[id]->minY = current.y;
			d[id]->maxZ = d[id]->minZ = current.z;
		} else {
			d[id]->minX = MIN(d[id]->minX, current.x);
			d[id]->maxX = MAX(d[id]->maxX, current.x);
			d[id]->minY = MIN(d[id]->minY, current.y);
			d[id]->maxY = MAX(d[id]->maxY, current.y);
			d[id]->minZ = MIN(d[id]->minZ, current.z);
			d[id]->minZ = MIN(d[id]->minZ, current.z);
		}
		for (j = 0; j < n_X; j++)
			current.X[j] = X[j * n + i];
		if (w != NULL)
			current.variance = 1.0/(w[i]);
		push_point(d[id], &current);
	}
	check_global_variables();
	d[id]->n_original = d[id]->n_list;
	efree(current.X);
	UNPROTECT(1); /* sy */
	return(sy);
}

SEXP gstat_new_dummy_data(SEXP loc_dim, SEXP has_intercept, SEXP beta, 
		SEXP nmax, SEXP nmin, SEXP maxdist, SEXP vfn, SEXP is_projected,
		SEXP vdist) {
	int i, id, dim, has_int;
	char name[NAME_SIZE];
	DATA **d = NULL;

	dim = INTEGER(loc_dim)[0];
	if (dim <= 0)
		Rf_error("dimension value impossible: %d", dim);
	if (dim > 3)
		Rf_error("too many dimensions: %d", dim);
	assert(LENGTH(beta) > 0);

	id = get_n_vars();
	snprintf(name, NAME_SIZE, "var%d", id);
	which_identifier(name);
	d = get_gstat_data();
	d[id]->id = id;

	d[id]->n_list = d[id]->n_max = 0;
	d[id]->colnx = d[id]->colny = d[id]->colnvalue = d[id]->colnz = 0;
	d[id]->x_coord = "x";
	d[id]->y_coord = "y";
	d[id]->z_coord = "z";
	d[id]->variable = "R data";
	d[id]->fname = "R data";
	has_int = INTEGER(has_intercept)[0];
	for (i = d[id]->n_X = 0; i < LENGTH(beta); i++)
		data_add_X(d[id], i + (has_int ? 0 : 1));
	assert(d[id]->n_X == LENGTH(beta));
	d[id]->dummy = 1;
	for (i = 0; i < LENGTH(beta); i++)
		d[id]->beta = push_d_vector(REAL(beta)[i], d[id]->beta);
	if (INTEGER(nmax)[0] > 0) /* leave default (large) if < 0 */
		d[id]->sel_max = INTEGER(nmax)[0];
/* I doubt whether using nmin for dummy data _ever_ can have a
 * meaning, but hey, let's add it anyway. */
	if (INTEGER(nmin)[0] > 0) /* leave default (0) if <= 0 */
		d[id]->sel_min = INTEGER(nmin)[0];
	if (REAL(maxdist)[0] > 0.0)
		d[id]->sel_rad = REAL(maxdist)[0];
	switch(INTEGER(vfn)[0]) {
		case 1: /* d[id]->variance_fn = v_identity; -> leave NULL */ break;
		case 2: d[id]->variance_fn = v_mu; break;
		case 3: d[id]->variance_fn = v_bin; break;
		case 4: d[id]->variance_fn = v_mu2; break;
		case 5: d[id]->variance_fn = v_mu3; break;
		default: Rf_error("unknown variance function %d", INTEGER(vfn)[0]);
	}
	gl_longlat = (INTEGER(is_projected)[0] == 0);
	d[id]->vdist = INTEGER(vdist)[0];
	d[id]->mode = X_BIT_SET | V_BIT_SET;
	if (dim > 1)
		d[id]->mode = d[id]->mode | Y_BIT_SET;
	if (dim > 2)
		d[id]->mode = d[id]->mode | Z_BIT_SET;
	set_norm_fns(d[id]); /* so gstat can calculate distances */
	check_global_variables();
	d[id]->n_original = d[id]->n_list;
	return(loc_dim);
}

SEXP gstat_predict(SEXP sn, SEXP slocs, SEXP sX, SEXP block_cols, SEXP block, 
			SEXP weights, SEXP nsim, SEXP blue) {
	double *locs, **est_all, *X;
	long i, j, k, n, nvars, nest, dim, n_X, ncols_block, 
		nrows_block, pos;
	DPOINT current, *bp = NULL;
	DATA **d = NULL, *vd = NULL, *area = NULL;
	SEXP ret;
	SEXP retvector;
	SEXP retvector_dim;
	extern unsigned int n_pred_locs; /* predict.c, used in msim.c */
	float ***msim = NULL;

	nvars = get_n_vars();
	nest = nvars + (nvars * (nvars + 1))/2;
	n = INTEGER(sn)[0];
	if (n <= 0 || LENGTH(slocs) == 0 || LENGTH(sX) == 0)
		ErrMsg(ER_IMPOSVAL, "newdata empty or only NA's");
	if (LENGTH(slocs) % n != 0)
		Rf_error("dimensions do not match: locations %d, nrows in X %ld",
			(int) LENGTH(slocs), n);
	dim = LENGTH(slocs) / n;
	if (dim > 3)
		Rf_error("too many spatial dimensions: %ld", dim);
	if (dim <= 0)
		Rf_error("too few spatial dimensions: %ld", dim);
	locs = REAL(slocs);
	if (LENGTH(sX) % n != 0)
		Rf_error("dimensions do not match: X %d and data %ld", (int) LENGTH(sX), n);
	n_X = LENGTH(sX) / n;

	current.attr = current.x = current.y = current.z = 0.0;
	current.bitfield = 0;
	/* assuming points ... */
	SET_POINT(&current);
	/* and then do the block thing: */
	if (LENGTH(block_cols) == 0) {
		bp = get_block_p();
		bp->x = bp->y = bp->z = 0.0; /* obsolete, I'd guess */
		if (LENGTH(block) >= 1) {
			bp->x = REAL(block)[0];
			SET_BLOCK(&current);
		}
		if (LENGTH(block) >= 2)
			bp->y = REAL(block)[1];
		if (LENGTH(block) >= 3)
			bp->z = REAL(block)[2];
		if (LENGTH(block) > 3)
			pr_warning("block dimension can only be 3; using the first 3");
	} else if (LENGTH(block_cols) == 1) { /* if > 1, block contains multiple 2D blocks */
		ncols_block = INTEGER(block_cols)[0];
		if (ncols_block < 1 || ncols_block > 3)
			ErrMsg(ER_IMPOSVAL, "block dimensions should be in [1..3]");
		nrows_block = LENGTH(block) / ncols_block; /* nr of rows */
		if (nrows_block > 0) {
			area = create_data_area();
			area->colnvariance = 0;
			area->n_list = area->n_max = 0;
			area->id = ID_OF_AREA;
			area->mode = X_BIT_SET;
			if (ncols_block > 1)
				area->mode = area->mode & Y_BIT_SET;
			if (ncols_block > 2)
				area->mode = area->mode & Z_BIT_SET;
			for (i = 0; i < nrows_block; i++) {
				current.x = REAL(block)[i];
				if (ncols_block > 1)
					current.y = REAL(block)[nrows_block + i];
				if (ncols_block > 2)
					current.z = REAL(block)[2 * nrows_block + i];
				if (LENGTH(weights) > 0) {
					area->colnvariance = 1;
					current.variance = REAL(weights)[i];
				}
				push_point(area, &current);
			}
			SET_BLOCK(&current);
		}
		if (DEBUG_FORCE)
			print_data_list(area);
	}

	X = REAL(sX);
	assert(n_X > 0);
	current.X = (double *) emalloc(n_X * sizeof(double));
	current.u.stratum = 0;
	d = get_gstat_data();
	est_all = (double **) emalloc(n * sizeof(double *));
	for (i = 0; i < n; i++)
		est_all[i] = (double *) emalloc(nest * sizeof(double));
	/* 
	 * the following is to fake gstat's default method handling: 
	 * we got to suggest that we'll go through a list of prediction
	 * locations, a la the gstat ``data(): ... ;'' command.
	 */
	vd = get_dataval(); 
	vd->id = ID_OF_VALDATA; 
	vd->mode = d[0]->mode;
	/* set min/max[XYZ] */
	vd->minY = vd->maxY = vd->minZ = vd->maxZ = 0.0;
	vd->minX = vd->maxX = locs[0];
	for (i = 1; i < n; i++) {
		vd->minX = MIN(vd->minX, locs[i]);
		vd->maxX = MAX(vd->maxX, locs[i]);
	}
	if (dim >= 2) {
		vd->minY = vd->maxY = locs[n];
		for (i = 1; i < n; i++) {
			vd->minY = MIN(vd->minY, locs[n + i]);
			vd->maxY = MAX(vd->maxY, locs[n + i]);
		}
	}
	if (dim >= 3) {
		vd->minZ = vd->maxZ = locs[2 * n];
		for (i = 1; i < n; i++) {
			vd->minZ = MIN(vd->minZ, locs[2 * n + i]);
			vd->maxZ = MAX(vd->maxZ, locs[2 * n + i]);
		}
	}

	/* fill, and standardize coordinate predictors from degree = x */
	for (i = 0; i < nvars; i++) 
		setup_data_minmax(d[i]);
	setup_data_minmax(vd);
	for (i = 0; i < nvars; i++) 
		calc_polynomials(d[i]);
	/* calc_polynomials(vd); */ /* still no data in fake vd */

	vd->polynomial_degree = d[0]->polynomial_degree;
	if (vd->polynomial_degree > 0) {
		setup_polynomial_X(vd); /* standardized coordinate polynomials */
		current.X = (double *) erealloc(current.X, vd->n_X * sizeof(double));
	}

	/* so far for the faking; now let's see what gstat makes out of this: */
	if (INTEGER(nsim)[0] == 0) {
		if (INTEGER(blue)[0] == 0) { /* FALSE */
			if (get_method() == NSP) /* choose default */
				set_method(get_default_method());
		} else 
			set_method(LSLM);
	} else {
		if (INTEGER(nsim)[0] < 0) {
			gl_nsim = -(INTEGER(nsim)[0]);
			set_method(ISI);
		} else {
			gl_nsim = INTEGER(nsim)[0];
			set_method(GSI);
		}
		n_pred_locs = n;
		if (gl_nsim > 1)
			init_simulations(d);
		if (get_n_beta_set() != get_n_vars())
			setup_beta(d, get_n_vars(), gl_nsim);
	}
	set_mode();  /* simple, stratified, multivariable? */
	check_global_variables(); /* it's there, better do it now */
	if (debug_level)
		Rprintf("[%s]\n", method_string(get_method()));
#ifdef WIN32
	R_FlushConsole();
	R_ProcessEvents();
#endif
	for (i = 0; i < n; i++) {
		print_progress(i, n);
		if (LENGTH(block_cols) > 1)
			gstat_set_block(i, block, block_cols, &current);
		current.x = locs[i];
		if (dim >= 2)
			current.y = locs[n + i];
		if (dim >= 3)
			current.z = locs[2 * n + i];
		for (j = 0; j < n_X; j++)
			current.X[j] = X[j * n + i];
		/* transform standardized coordinate polynomial here */
		if (vd->polynomial_degree)
			calc_polynomial_point(vd, &current);
		for (j = 0; j < get_n_vars(); j++)
			select_at(d[j], &current);
		get_est(d, get_method(), &current, est_all[i]);
#ifdef WIN32
		R_ProcessEvents(); /* avoid terminal freeze in R/Win */
#endif
		R_CheckUserInterrupt();
	}
	print_progress(100, 100);
	PROTECT(ret = Rf_allocVector(VECSXP, 1));
	PROTECT(retvector_dim = Rf_allocVector(REALSXP, 2));
	REAL(retvector_dim)[0] = n; /* nrows */
	if (gl_nsim > 1) {
		PROTECT(retvector = Rf_allocVector(REALSXP, gl_nsim * nvars * n));
		msim = get_msim();
		for (i = pos = 0; i < nvars; i++)
			for (j = 0; j < gl_nsim; j++) 
				for (k = 0; k < n; k++) {
					if (is_mv_float(&(msim[i][k][j])))
						REAL(retvector)[pos++] = NA_REAL;
					else
						REAL(retvector)[pos++] = msim[i][k][j];
				}
		REAL(retvector_dim)[1] = nvars * gl_nsim; /* ncols */
	} else {
		PROTECT(retvector = Rf_allocVector(REALSXP, n * nest));
		for (j = pos = 0; j < nest; j++) {
			for (i = 0; i < n; i++) {
				if (is_mv_double(&(est_all[i][j])))
					REAL(retvector)[pos] = NA_REAL;
				else
					REAL(retvector)[pos] = est_all[i][j];
				pos++;
			}
		}
		REAL(retvector_dim)[1] = nest; /* ncols */
	}
	if (gl_nsim > 0)
		free_simulations();
	/* SET_DIM(retvector, retvector_dim); */
	Rf_setAttrib(retvector, R_DimSymbol, retvector_dim);
	SET_VECTOR_ELT(ret, 0, retvector);
	for (i = 0; i < n; i++)
		efree(est_all[i]);
	efree(est_all);
	efree(current.X);
	UNPROTECT(3);
	return(ret);
}

static void gstat_set_block(long i, SEXP block, SEXP block_cols, DPOINT *current) {
	DATA *area;
	VARIOGRAM *v;
	long nrows_block, start, end, j;

	if (i >= LENGTH(block_cols) || i < 0)
		ErrMsg(ER_IMPOSVAL, "block_cols length less than nr of prediction locations");
	nrows_block = LENGTH(block) / 2; /* nr of rows */
	start = INTEGER(block_cols)[i];
	if (i == LENGTH(block_cols) - 1)
		end = nrows_block;
	else
		end = INTEGER(block_cols)[i+1] - 1;
	area = get_data_area();
	if (area != NULL)
		free_data(area);
	area = create_data_area();
	area->n_list = area->n_max = 0;
	area->id = ID_OF_AREA;
	area->mode = X_BIT_SET & Y_BIT_SET;
	for (j = start - 1; j < end; j++) {
		current->x = REAL(block)[j];
		current->y = REAL(block)[nrows_block + j];
		push_point(area, current);
	}
	SET_BLOCK(current);
	if (DEBUG_FORCE)
		print_data_list(area); 
	for (j = 0; j < get_n_vgms(); j++) {
		v = get_vgm(j);
		if (v != NULL)
			v->block_semivariance_set = v->block_covariance_set = 0; /* don't store under these circumstances! */
	}
	return;
}

SEXP gstat_variogram(SEXP s_ids, SEXP cutoff, SEXP width, SEXP direction, 
		SEXP cressie, SEXP dX, SEXP boundaries, SEXP grid, SEXP cov,
		SEXP pseudo) {
	SEXP ret;
	SEXP np; 
	SEXP dist;
	SEXP gamma;
	SEXP sx;
	SEXP sy;
	SEXP ev_parameters;
	/* SEXP y; */
	long i, id1, id2, nest;
	VARIOGRAM *vgm;
	DATA **d;

	GRIDMAP *m;
	unsigned int row, col, n;

	id1 = INTEGER(s_ids)[0];
	if (LENGTH(s_ids) > 1)
		id2 = INTEGER(s_ids)[1];
	else
		id2 = id1;
	vgm = get_vgm(LTI(id1,id2));
	vgm->id = LTI(id1,id2);
	vgm->id1 = id1;
	vgm->id2 = id2;
	if (INTEGER(cov)[0] == 0)
		vgm->ev->evt = (id1 == id2 ? SEMIVARIOGRAM : CROSSVARIOGRAM);
	else if (INTEGER(cov)[0] == 1)
		vgm->ev->evt = (id1 == id2 ? COVARIOGRAM : CROSSCOVARIOGRAM);
	else {
		if (id1 != id2)
			ErrMsg(ER_IMPOSVAL,
			"cannot compute pairwise relative cross semivariogram");
		if (INTEGER(cov)[0] == 2)
			vgm->ev->evt = PRSEMIVARIOGRAM;
	}
	/* vgm->ev->is_asym = INTEGER(asym)[0]; */
	vgm->ev->pseudo = INTEGER(pseudo)[0];
	vgm->ev->recalc = 1;
	if (LENGTH(cutoff) > 0)
		gl_cutoff = REAL(cutoff)[0];
	if (LENGTH(width) > 0)
		gl_iwidth = REAL(width)[0];
	gl_alpha = REAL(direction)[0];
	gl_beta = REAL(direction)[1];
	gl_tol_hor = REAL(direction)[2];
	gl_tol_ver = REAL(direction)[3];
	gl_cressie = INTEGER(cressie)[0];
	if (LENGTH(dX) > 0) {
		d = get_gstat_data();
		d[id1]->dX = REAL(dX)[0];
		d[id2]->dX = REAL(dX)[0];
	} 
	for (i = 0; i < LENGTH(boundaries); i++) /* does nothing if LENGTH is 0 */
		push_bound(REAL(boundaries)[i]);
	switch (LENGTH(grid)) {
		case 0: case 1: break;
		case 6: vgm->ev->S_grid = gstat_S_fillgrid(grid); break;
		default: Rf_error("unrecognized grid length in gstat_variogram");
			break;
	}

	calc_variogram(vgm, NULL);

	if (vgm->ev->S_grid != NULL) {
		PROTECT(ret = Rf_allocVector(VECSXP, 4));
		m = vgm->ev->map;
		n = m->rows * m->cols;
		PROTECT(np = Rf_allocVector(REALSXP, n));
		PROTECT(gamma = Rf_allocVector(REALSXP, n));
		PROTECT(sx = Rf_allocVector(REALSXP, n));
		PROTECT(sy = Rf_allocVector(REALSXP, n));

		for (row = i = 0; row < m->rows; row++) {
			for (col = 0; col < m->cols; col++) {
				map_rowcol2xy(m, row, col, &(REAL(sx)[i]), 
								&(REAL(sy)[i]));
				REAL(np)[i] = vgm->ev->nh[i];
				if (vgm->ev->nh[i] > 0)
					REAL(gamma)[i] = vgm->ev->gamma[i];
				else 
					REAL(gamma)[i] = NA_REAL;
				i++;
			}
		}
		SET_VECTOR_ELT(ret, 0, sx);
		SET_VECTOR_ELT(ret, 1, sy);
		SET_VECTOR_ELT(ret, 2, np);
		SET_VECTOR_ELT(ret, 3, gamma);
		free_data_gridmap(vgm->ev->S_grid);
		UNPROTECT(5);
	} else {
		if (vgm->ev->cloud)
			nest = vgm->ev->n_est;
		else {
			if (vgm->ev->zero == ZERO_SPECIAL)
				nest = vgm->ev->n_est;
			else 
				nest = vgm->ev->n_est - 1;
		}
		PROTECT(ret = Rf_allocVector(VECSXP, 4));
		if (nest <= 0) {
			UNPROTECT(1);
			return(ret);
		}
		PROTECT(np = Rf_allocVector(REALSXP, nest));
		PROTECT(dist = Rf_allocVector(REALSXP, nest));
		PROTECT(gamma = Rf_allocVector(REALSXP, nest));
		PROTECT(ev_parameters = Rf_allocVector(REALSXP, 4));
		REAL(ev_parameters)[0] = vgm->ev->cutoff;
		REAL(ev_parameters)[1] = vgm->ev->iwidth;
		REAL(ev_parameters)[2] = vgm->ev->pseudo;
		REAL(ev_parameters)[3] = vgm->ev->is_asym;
		for (i = 0; i < nest; i++) {
			REAL(np)[i] = vgm->ev->nh[i];
			REAL(dist)[i] = vgm->ev->dist[i];
			REAL(gamma)[i] = vgm->ev->gamma[i];
		}
		SET_VECTOR_ELT(ret, 0, np);
		SET_VECTOR_ELT(ret, 1, dist);
		SET_VECTOR_ELT(ret, 2, gamma);
		SET_VECTOR_ELT(ret, 3, ev_parameters);
		UNPROTECT(5);
	}
	return(ret);
}

SEXP gstat_load_variogram(SEXP s_ids, SEXP s_model, SEXP s_sills, SEXP s_ranges, 
		SEXP s_kappas, SEXP s_anis_all, SEXP s_table, SEXP s_max_val) 
{
	VARIOGRAM *vgm;
	long i, n, id1, id2, max_id;
	double anis[5] = {0.0, 0.0, 0.0, 1.0, 1.0}, rpars[2], *sills, *ranges, 
		*kappas, *anis_all;
	const char *model;

	sills = REAL(s_sills);
	ranges = REAL(s_ranges);
	kappas = REAL(s_kappas);
	anis_all = REAL(s_anis_all);

	id1 = INTEGER(s_ids)[0];
	id2 = INTEGER(s_ids)[1];
	max_id = MAX(id1, id2);

	if (get_n_vars() == 0)
		which_identifier("xx"); /* at least "load" one dummy var */
	if (max_id >= get_n_vars())
		ErrMsg(ER_IMPOSVAL,
			"gstat_load_variogram has been called with max_id >= n_vars");

	vgm = get_vgm(LTI(id1,id2));
	assert(vgm != NULL);

	vgm->id = LTI(id1,id2);
	vgm->id1 = id1;
	vgm->id2 = id2;
	vgm->n_models = vgm->n_fit = 0;

	n = LENGTH(s_sills);
	for (i = 0; i < n; i++) { /* loop over sub models */
		model = CHAR(STRING_ELT(s_model, i));
		anis[0] = anis_all[0 * n + i];
		anis[1] = anis_all[1 * n + i];
		anis[2] = anis_all[2 * n + i];
		anis[3] = anis_all[3 * n + i];
		anis[4] = anis_all[4 * n + i];
		rpars[0] = ranges[i];
		rpars[1] = kappas[i];
		if (LENGTH(s_table) > 0)
			push_to_v_table(vgm, rpars[0], 
					LENGTH(s_table), REAL(s_table),
					(anis[3] == 1.0 && anis[4] == 1.0) ? NULL : anis);
		else
			push_to_v(vgm, model, sills[i], rpars, 2,
				(anis[3] == 1.0 && anis[4] == 1.0) ? NULL : anis, 1, 1);
	}
	update_variogram(vgm);
	if (REAL(s_max_val)[0] > 0.0 || REAL(s_max_val)[1] > 0.0 || REAL(s_max_val)[2] > 0.0)
		vgm->max_val = get_semivariance(vgm, 
				REAL(s_max_val)[0], REAL(s_max_val)[1], REAL(s_max_val)[2]);
	if (DEBUG_DUMP)
		logprint_variogram(vgm, 1); 
	return(s_model);
}

SEXP gstat_variogram_values(SEXP ids, SEXP pars, SEXP covariance, SEXP dist_values) {
	double from, to, n, d, x = 1.0, y = 0.0, z = 0.0;
	int i, id1, id2, cov = 0, ndist = 0;
	VARIOGRAM *vgm;
	SEXP dist;
	SEXP gamma;
	SEXP ret;

	if (LENGTH(pars) != 3 && LENGTH(pars) != 6)
		Rf_error("supply three or six distance parameters");
	from = REAL(pars)[0];
	to = REAL(pars)[1];
	n = REAL(pars)[2];
	ndist = LENGTH(dist_values);
	cov = INTEGER(covariance)[0];
	if (LENGTH(pars) == 6) {
		x = REAL(pars)[3];
		y = REAL(pars)[4];
		z = REAL(pars)[5];
	}

	id1 = INTEGER(ids)[0];
	id2 = INTEGER(ids)[1];
	vgm = get_vgm(LTI(id1,id2));

	if (ndist > 0) {
		PROTECT(dist = Rf_allocVector(REALSXP, ndist));
		PROTECT(gamma = Rf_allocVector(REALSXP, ndist));
		for (i = 0; i < ndist; i++) {
			d = REAL(dist_values)[i];
			REAL(dist)[i] = d;
			REAL(gamma)[i] = (cov ? 
				get_covariance(vgm, d * x, d * y, d * z) : 
				get_semivariance(vgm, d * x, d * y, d * z));
		}
	} else {
		PROTECT(dist = Rf_allocVector(REALSXP, n));
		PROTECT(gamma = Rf_allocVector(REALSXP, n));
		for (i = 0; i < n; i++) {
			d = from;
			if (i > 0) /* implies n > 1 */
				d += (i/(n-1))*(to-from);
			REAL(dist)[i] = d;
			REAL(gamma)[i] = (cov ? 
				get_covariance(vgm, d * x, d * y, d * z) : 
				get_semivariance(vgm, d * x, d * y, d * z));
		}
	}
	PROTECT(ret = Rf_allocVector(VECSXP, 2));
	SET_VECTOR_ELT(ret, 0, dist);
	SET_VECTOR_ELT(ret, 1, gamma);
	UNPROTECT(3);
	return(ret);
}

// Added by Paul Hiemstra, 30-06-2008
SEXP get_covariance_list(SEXP ids, SEXP covariance, SEXP dist_list) {
	double d, x = 1.0, y = 0.0, z = 0.0;
	int i, id1, id2, cov = 0;
	VARIOGRAM *vgm;
	SEXP dist;
	SEXP gamma;
	SEXP ret;
	int length_list = LENGTH(dist_list);

	cov = INTEGER(covariance)[0];

	id1 = INTEGER(ids)[0];
	id2 = INTEGER(ids)[1];
	vgm = get_vgm(LTI(id1,id2));

	PROTECT(dist = Rf_allocVector(REALSXP, length_list));
	PROTECT(gamma = Rf_allocVector(REALSXP, length_list));
	for (i = 0; i < length_list; i++) {
		d = REAL(dist_list)[i];
		REAL(dist)[i] = d;
		REAL(gamma)[i] = (cov ? 
			get_covariance(vgm, d * x, d * y, d * z) : 
			get_semivariance(vgm, d * x, d * y, d * z));
	}
	PROTECT(ret = Rf_allocVector(VECSXP, 2));
	SET_VECTOR_ELT(ret, 0, dist);
	SET_VECTOR_ELT(ret, 1, gamma);
	UNPROTECT(3);
	return(ret);
}

SEXP gstat_get_variogram_models(SEXP dolong) {
	SEXP ret;
	int i, n = 0, do_long;
	
	for (i = 1; v_models[i].model != NOT_SP; i++)
		n++;

	do_long = INTEGER(dolong)[0];
	PROTECT(ret = Rf_allocVector(STRSXP, n));
	for (i = 1; v_models[i].model != NOT_SP; i++)
		SET_STRING_ELT(ret, i-1, 
				Rf_mkChar(do_long ? v_models[i].name_long : v_models[i].name));
	UNPROTECT(1);
	return(ret);
}

SEXP gstat_load_ev(SEXP np, SEXP dist, SEXP gamma) {

	int i, cloud = 1;
	VARIOGRAM *vgm;

	which_identifier("xx");
	/*
	 * vgm = get_vgm(LTI(INTEGER(id)[0], INTEGER(id)[1]));
	 * */
	vgm = get_vgm(LTI(0, 0));
	if (vgm->ev == NULL)
		vgm->ev = init_ev();
	vgm->ev->evt = SEMIVARIOGRAM;
	vgm->ev->n_est = LENGTH(np);
	vgm->ev->n_max = LENGTH(np);
	vgm->ev->gamma = (double *) emalloc (sizeof(double) * vgm->ev->n_max);
	vgm->ev->dist = (double *) emalloc (sizeof(double) * vgm->ev->n_max);
	vgm->ev->nh = (unsigned long *) emalloc (sizeof(long) * vgm->ev->n_max);
	for (i = 0; i < vgm->ev->n_est; i++) {
		vgm->ev->nh[i] = REAL(np)[i];
		vgm->ev->dist[i] = REAL(dist)[i];
		vgm->ev->gamma[i] = REAL(gamma)[i];
		if (cloud && vgm->ev->nh[i] > 1)
			cloud = 0;
	}
	vgm->ev->cloud = cloud;
	if (DEBUG_VGMFIT)
		fprint_sample_vgm(vgm->ev);
	return(np);
}

SEXP gstat_fit_variogram(SEXP fit, SEXP fit_sill, SEXP fit_range) {
	int i;
	VARIOGRAM *vgm;
	SEXP ret;
	SEXP sills;
	SEXP ranges;
	SEXP SSErr;
	SEXP fit_is_singular;

	vgm = get_vgm(LTI(0, 0));
	vgm->ev->fit = INTEGER(fit)[0];
	for (i = 0; i < vgm->n_models; i++) {
		vgm->part[i].fit_sill = INTEGER(fit_sill)[i];
		vgm->part[i].fit_range = INTEGER(fit_range)[i];
	}
	update_variogram(vgm);
	if (DEBUG_VGMFIT)
		logprint_variogram(vgm, 1);
	fit_variogram(vgm);
	if (DEBUG_VGMFIT)
		logprint_variogram(vgm, 1);

	PROTECT(sills = Rf_allocVector(REALSXP, vgm->n_models));
	PROTECT(ranges = Rf_allocVector(REALSXP, vgm->n_models));
	for (i = 0; i < vgm->n_models; i++) {
		REAL(sills)[i] = vgm->part[i].sill;
		REAL(ranges)[i] = vgm->part[i].range[0];
	}

	PROTECT(ret = Rf_allocVector(VECSXP, 4));
	SET_VECTOR_ELT(ret, 0, sills);
	SET_VECTOR_ELT(ret, 1, ranges);

	PROTECT(fit_is_singular = Rf_allocVector(REALSXP, 1));
	REAL(fit_is_singular)[0] = vgm->fit_is_singular;
	SET_VECTOR_ELT(ret, 2, fit_is_singular);

	PROTECT(SSErr = Rf_allocVector(REALSXP, 1));
	REAL(SSErr)[0] = vgm->SSErr;
	SET_VECTOR_ELT(ret, 3, SSErr);

	UNPROTECT(5);
	return(ret);
}

SEXP gstat_debug_level(SEXP level) {
	debug_level = INTEGER(level)[0];
	return(level);
}

SEXP gstat_set_method(SEXP to) {
	const char *what;

	what = CHAR(STRING_ELT(to, 0));
	for (int id = 1; methods[id].name != NULL; id++) {
		if (almost_equals(what, methods[id].name)) {
			set_method(methods[id].m);
			break; /* id-loop */
		}
	}
	return(to);
}

SEXP gstat_set_set(SEXP arg, SEXP val) {
	const char *name;
	int i;
	typedef struct {
		const char *name;
		void *ptr;
		enum { UNKNOWN, IS_INT, IS_UINT, IS_REAL, IS_STRING, IS_D_VECTOR, NO_ARG } what;
		enum { NOLIMIT, GEZERO, GTZERO } limit;
	} GSTAT_EXPR;
	const GSTAT_EXPR set_options[] = {
	{ "alpha",          &gl_alpha,        IS_REAL, GEZERO  },
	{ "beta",           &gl_beta,         IS_REAL, GEZERO  },
	{ "blas",           &gl_blas,         IS_INT,  GEZERO  },
	{ "choleski",       &gl_choleski,     IS_INT,  GEZERO  },
	{ "co$incide",      &gl_coincide,     IS_INT,  GEZERO  },
	{ "Cr$essie",       &gl_cressie,      IS_INT,  GEZERO  },
	{ "cutoff",         &gl_cutoff,       IS_REAL, GTZERO  },
	{ "de$bug",         &debug_level,     IS_INT,  GEZERO  },
	{ "fit",            &gl_fit,          IS_INT,  GEZERO  },
	{ "fit_l$imit",     &gl_fit_limit,    IS_REAL, GTZERO  },
	{ "fr$action",      &gl_fraction,     IS_REAL, GTZERO  },
	/* { "display",        &gl_display,      IS_STRING, NOLIMIT }, */
	{ "gls$_residuals", &gl_gls_residuals, IS_INT, GEZERO  },
	{ "id$p",           &gl_idp,          IS_REAL, GEZERO  },
	{ "in$tervals",     &gl_n_intervals,  IS_INT,  GTZERO  },
	{ "it$er",          &gl_iter,         IS_INT,  GEZERO  },
	{ "lhs",            &gl_lhs,          IS_INT,  GEZERO  },
	{ "longlat",        &gl_longlat,      IS_INT, GEZERO },
	{ "sim_beta",  		&gl_sim_beta,     IS_INT, GEZERO },
	{ "n_uk",           &gl_n_uk,         IS_INT,  GEZERO  },
	{ "numbers",        &gl_numbers,      IS_INT,  GEZERO  },
	{ "nb$lockdiscr",   &gl_nblockdiscr,  IS_INT,  GTZERO  },
	{ "no$check",       &gl_nocheck,      IS_INT,  GEZERO  },
	{ "ns$im",          &gl_nsim,         IS_INT,  GTZERO  },
	{ "or$der",         &gl_order,        IS_INT,  GEZERO },
	{ "q$uantile",      &gl_quantile,     IS_REAL, GEZERO  },
	{ "rowwise",        &gl_rowwise,      IS_INT,  GEZERO  },
	{ "rp",             &gl_rp,           IS_INT,  GEZERO  },
	{ "see$d",          &gl_seed,         IS_INT,  GTZERO  },
	{ "useed",          &gl_seed,         IS_UINT,  GEZERO  },
	{ "spa$rse",        &gl_sparse,       IS_INT,  GEZERO  },
	{ "spi$ral",        &gl_spiral,       IS_INT,  GEZERO  },
	{ "spl$it",         &gl_split,        IS_INT,  GTZERO  },
	{ "sy$mmetric",     &gl_sym_ev,       IS_INT,  GEZERO  },
	{ "tol_h$or",       &gl_tol_hor,      IS_REAL, GEZERO  },
	{ "tol_v$er",       &gl_tol_ver,      IS_REAL, GEZERO  },
	{ "v$erbose",       &debug_level,     IS_INT,  GEZERO  },
	{ "w$idth",         &gl_iwidth,       IS_REAL, GEZERO  },
	{ "x$valid",        &gl_xvalid,       IS_INT,  GEZERO  },
	{ "zero_di$st",     &gl_zero_est,     IS_INT,  GEZERO  },
	{ "zero",           &gl_zero,         IS_REAL, GEZERO  },
	{ "zm$ap",          &gl_zmap,         IS_REAL, NOLIMIT },
	{ NULL, NULL, 0, 0 }
	};
	name = CHAR(STRING_ELT(arg, 0));
	for (i = 0; set_options[i].name; i++)
		if (almost_equals(name, set_options[i].name))
			break; /* break out i-loop */
	if (set_options[i].name == NULL)
		ErrMsg(ER_SYNTAX, name);

	if (almost_equals((const char *)name, "nb$lockdiscr"))
		gl_gauss = 0; /* side effect */

	switch (set_options[i].what) {
		case IS_INT: 
			*((int *) set_options[i].ptr) = Rf_asInteger(val);
			/* Rprintf("int arg: %s val %d\n", name, Rf_asInteger(val)); */
			break;
		case IS_UINT: 
			*((unsigned int *) set_options[i].ptr) = (unsigned int) Rf_asInteger(val);
			/* Rprintf("uint arg: %s val %d\n", name, Rf_asInteger(val)); */
			break;
		case IS_REAL: 
			*((double *) set_options[i].ptr) = Rf_asReal(val);
			/* Rprintf("real arg: %s val %d\n", name, asReal(val)); */
			break; 
		case IS_STRING: 
			*((const char **) set_options[i].ptr) = CHAR(STRING_ELT(val, 0));
			break;
		default:
			ErrMsg(ER_SYNTAX, name);
			break;
	}
	return val;
}

SEXP gstat_set_merge(SEXP a, SEXP b, SEXP c, SEXP d) { /* merge a(b) with c(d); */
	DATA **dpp;
	int id, id1, id2, col1, col2;
	id1 = Rf_asInteger(a);
	id2 = Rf_asInteger(c);
	if (id1 >= get_n_vars() || id2 >= get_n_vars() || id1 < 0 || id2 < 0)
		ErrMsg(ER_IMPOSVAL, "id values out of range");
	col1 = Rf_asInteger(b);
	col2 = Rf_asInteger(d);
	if (id1 < id2) { /* swap id and col */
		id = id1; id1 = id2; id2 = id;
		id = col1; col1 = col2; col2 = id;
	}
	dpp = get_gstat_data();
	if (push_to_merge_table(dpp[id1], id2, col1, col2)) 
		ErrMsg(ER_IMPOSVAL, "attempt to merge failed");
	return(a);
}

double r_uniform(void) {
	return(unif_rand());
}

double r_normal(void) {
	return(norm_rand());
}

static DATA_GRIDMAP *gstat_S_fillgrid(SEXP gridparams) {
	double x_ul, y_ul, cellsizex, cellsizey;
	unsigned int rows, cols;

	cellsizex = REAL(gridparams)[2];
	cellsizey = REAL(gridparams)[3];
	rows = (unsigned int) REAL(gridparams)[5];
	cols = (unsigned int) REAL(gridparams)[4];
	x_ul = REAL(gridparams)[0] - 0.5 * cellsizex;
	y_ul = REAL(gridparams)[1] + (rows - 0.5) * cellsizey;
	return gsetup_gridmap(x_ul, y_ul, cellsizex, cellsizey, rows, cols);
}
