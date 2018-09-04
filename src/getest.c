/*
 * getest.c: choose the predicion function at one prediction location
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>

#include "defs.h"
#include "data.h"
#include "utils.h"
#include "mapio.h"
#include "userio.h"
#include "debug.h"
#include "vario.h"
#include "sem.h"
#include "fit.h"
#include "glvars.h"
#include "mtrx.h"
#include "lm.h"
#include "gls.h"
#include "sim.h"
#include "msim.h"
#include "block.h"
#include "getest.h"

static void est_quantile_div(DATA *data, double *est, int div);
static void est_skew_kurt(DATA *data, double *est);
static double inverse_dist(DATA *data, DPOINT *where, double idPow);
static void save_variogram_parameters(VARIOGRAM *v);
static void reset_variogram_parameters(VARIOGRAM *v);
static double sample_mean(double *list, int n);
static double sample_var(double *list, double mean, int n);
static double sample_std(double *list, double mean, int n);
static double est_quant(double *list, double p, int n);
static int CDECL d_cmp(const double *a, const double *b);
static double *vgm_pars = NULL;

void get_est(DATA **data, METHOD method, DPOINT *where, double *est) {
/*
 * given all data[i]->sel, get_est returns one or more values to *est,
 * according to the method of calculation set in *method 
 * sel must be not NULL, where contains the location of the estimate
 * PRE: data, where, est;
 * USES: several estimation routines.
 */
	DPOINT *block = NULL;
	VARIOGRAM *v;
	int i, j, n_vars, n_sel, *is_pt;
	double *X_ori = NULL, *local_sim;
	enum GLS_WHAT gls_mode = GLS_BLUP;
	const double *sim = NULL; /* return value of cond_sim() */

	for (i = 0; i < get_n_outputs(); i++)
		set_mv_double(&est[i]);
	block = get_block_p();
	if (get_mode() == MODE_NSP)
		ErrMsg(ER_IMPOSVAL, "Getest(): mode not specified");
	if (block->x > 0.0 || block->y > 0.0 || block->z > 0.0 || get_data_area())
		SET_BLOCK(where);
	else
		SET_POINT(where);
	n_vars = get_n_vars();
	if (get_mode() == STRATIFY &&
			(where->u.stratum < 0 || where->u.stratum >= n_vars))
		return;
	local_sim = (double *) emalloc(n_vars * sizeof(double));
	is_pt = (int *) emalloc(n_vars * sizeof(int));
	for (i = 0; i < n_vars; i++) {
		set_mv_double(&(local_sim[i]));
		is_pt[i] = 0;
	}
	if (DEBUG_COV) {
		printlog("we're at location X: %g Y: %g Z: %g\n",
			where->x, where->y, where->z);
		if (IS_BLOCK(where)) {
			if (get_data_area())
				printlog("block set in area()\n");
			else
				printlog("block size: dx: %g dy: %g dz: %g\n",
					block->x, block->y, block->z);
		} else
			printlog("zero block size\n");
		if (get_mode() == STRATIFY)
			printlog("stratum: %d\n", where->u.stratum);
	}
	switch (method) {
		case DIV: /* BREAKTHROUGH */
		case MED: 
			if (get_mode() == STRATIFY)
				est_quantile_div(data[where->u.stratum], est, method == DIV);
			else
				for (i = 0; i < n_vars; i++)
					est_quantile_div(data[i], &(est[2*i]), method == DIV);
			break;
		case SKEW: 
			if (get_mode() == STRATIFY)
				est_skew_kurt(data[where->u.stratum], est);
			else
				for (i = 0; i < n_vars; i++)
					est_skew_kurt(data[i], &(est[2*i]));
			break;
		case IDW: 
			if (gl_idp < 0.0)
				ErrMsg(ER_RANGE, "idp must be non-negative");
			if (get_mode() == STRATIFY) {
				if (data[where->u.stratum]->n_sel > 0)
					est[0] = inverse_dist(data[where->u.stratum], where, gl_idp);
			} else {
				for (i = 0; i < n_vars; i++) 
					if (data[i]->n_sel > 0)
						est[2 * i] = inverse_dist(data[i], where, gl_idp);
			}
			break;
		case LSLM:
			if (n_variograms_set()) {
				switch(get_mode()) {
					case SIMPLE: 
						X_ori = where->X; /* remember... */
						for (i = 0; i < n_vars; i++) {
							if (data[i]->n_sel > 0)
								gls(&data[i], 1, GLS_BLUE, where,
									&(est[2 * i])); 
							where->X += data[i]->n_X;
						}
						where->X = X_ori; /* and put back */
						break;
					case STRATIFY:
						X_ori = where->X; /* remember... */
						for (i = 0; i < where->u.stratum; i++)
							where->X += data[i]->n_X; 
						if (data[where->u.stratum]->n_sel > 0) {
							where->X += where->u.stratum; /* strictly unnecessary */
							gls(&data[where->u.stratum], 1, GLS_BLUE, where, est); 
						}
						where->X = X_ori; /* and put back */
						break;
					case MULTIVARIABLE:
						gls(data, n_vars, GLS_BLUE, where, est); 
						break;
					default: ErrMsg(ER_IMPOSVAL, "Getest(): wrong mode");
						break;
				}
			} else {
				switch(get_mode()) {
					case STRATIFY:
						X_ori = where->X; /* remember... */
						for (i = 0; i < where->u.stratum; i++)
							where->X += data[i]->n_X; 
						if (data[where->u.stratum]->n_sel > 0)
							pred_lm(&data[where->u.stratum], 1, where, est);
						where->X = X_ori; /* and put back */
						break;
					case MULTIVARIABLE:
						pred_lm(data, n_vars, where, est);
						break;
					case SIMPLE:
						X_ori = where->X; /* remember... */
						for (i = 0; i < n_vars; i++) {
							if (data[i]->n_sel > 0) {
								pred_lm(&data[i], 1, where, &(est[2*i]));
								where->X += data[i]->n_X;
							}
						}
						where->X = X_ori; /* and put back */
						break;
					default: ErrMsg(ER_IMPOSVAL, "Getest(): wrong mode");
						break;
				}
			} /* else */
			break;
		case SKR: /* BREAK THROUGH (SK handled in gls()) */
		case OKR: /* BREAK THROUGH: (OK is special case of UK) */
		case UKR:
			if (method == SKR)
				gls_mode = GLS_BLP;
			else
				gls_mode = GLS_BLUP;
			switch (get_mode()) {
				case SIMPLE: 
					X_ori = where->X; /* remember... */
					for (i = 0; i < n_vars; i++) {
						gls(&data[i], 1, gls_mode, where, &(est[2 * i])); 
						where->X += data[i]->n_X;
					}
					where->X = X_ori; /* and put back */
					break;
				case STRATIFY:
					X_ori = where->X; /* remember... */
					for (i = 0; i < where->u.stratum; i++)
						where->X += data[i]->n_X; 
					gls(&data[where->u.stratum], 1, gls_mode, where, est); 
					where->X = X_ori; /* and put back */
					break;
				case MULTIVARIABLE:
					gls(data, n_vars, gls_mode, where, est); 
					break;
				default: ErrMsg(ER_IMPOSVAL, "Getest(): wrong mode"); break;
			}
			if (gl_order) /* order relation violations: */
				correct_orv(est, n_vars, gl_order); 
			break;
		case GSI: /* BREAK TRHOUGH: */
		case ISI: 
			switch (get_mode()) {
				case SIMPLE: /* estimates go to est[0+1],est[2+3],... */
					X_ori = where->X; /* remember... */
					for (i = 0; i < gl_nsim; i++) {
						restore_data_sel(data, i, n_vars);
						set_beta(data, i, n_vars, method);
						for (j = 0; j < n_vars; j++) {
							/* estimate or update: */
							if (i == 0) /* estimate */
								gls(&data[j], 1, data[j]->n_sel < gl_n_uk ?
									GLS_BLP : GLS_BLUP, where, &(est[2*j]));
							else /* update */
								gls(&data[j], 1, UPDATE, where, &(est[2*j]));
							/* new simulation: */
							if (gl_order <= 1) {
								sim = cond_sim(&(est[2*j]), 1, method, 
										&(is_pt[j]), gl_order);
								/* put data from simulation i; */
								local_sim[j] = sim[0];
							}
							where->X += data[j]->n_X;
						}
						where->X = X_ori; /* reset */
						if (gl_order > 1) { /* 2,3 if pdf, 4 if cdf: */
							sim = cond_sim(est, n_vars, method, is_pt, gl_order);
							for (j = 0; j < n_vars; j++)
								local_sim[j] = sim[j];
						}
						save_sim(data, where, i, n_vars, local_sim, is_pt);
					}
					X_ori = where->X; /* remember... */
					for (i = 0; i < n_vars; i++) {
						if (! is_pt[i]) {
							where->attr = est[2 * i] = local_sim[i]; 
								/* last simulation */
							push_point(data[i], where);
						} else
							data[i]->nsim_at_data++;
						where->X += data[i]->n_X;
					}
					where->X = X_ori; /* reset */
					break;
				case STRATIFY: /* estimate goes to est[0] and est[1] */
					X_ori = where->X; /* remember... */
					for (i = 0; i < where->u.stratum; i++)
						where->X += data[i]->n_X; 
					for (i = 0; i < gl_nsim; i++) {
						/* load data from simulation i; */
						restore_data_sel(&data[where->u.stratum], i, 0);
						set_beta(&data[where->u.stratum], i, 0, method);
						if (i == 0)  /* estimate: */
							gls(&data[where->u.stratum], 1, data[where->u.stratum]->n_sel <
								gl_n_uk ? GLS_BLP : GLS_BLUP, where, est);
						else /* update */
							gls(&data[where->u.stratum], 1, UPDATE, where, est);
						sim = cond_sim(est, 1, method, is_pt, gl_order);
						/* save data from simulation i; */
						save_sim_strat(data[where->u.stratum], where, i, 
								sim[0], is_pt[0]);
					}
					/* store & push */
					if (! is_pt[0]) {
						where->attr = est[0] = sim[0]; /* last simulation */
						push_point(data[where->u.stratum], where);
					} else
						data[where->u.stratum]->nsim_at_data++;
					where->X = X_ori; /* and put back */
					break;
				case MULTIVARIABLE:
					for (i = 0, n_sel = 0; i < n_vars; i++) 
						n_sel += data[i]->n_sel;
					if (n_sel < gl_n_uk)
						gls_mode = GLS_BLP;
					else
						gls_mode = GLS_BLUP;
					for (i = 0; i < gl_nsim; i++) {
						/* load data from simulation i; */
						restore_data_sel(data, i, n_vars);
						set_beta(data, i, n_vars, method);
						/* estimate or update: */
						if (i == 0)
							gls(data, n_vars, gls_mode, where, est);
						else
							gls(data, n_vars, UPDATE, where, est);
						/* new simulation: */
						sim = cond_sim(est, n_vars, method, is_pt, gl_order);
						/* put data from simulation i; */
						save_sim(data, where, i, n_vars, sim, is_pt);
					}
					X_ori = where->X; /* remember... */
					for (i = 0; i < n_vars; i++) {
						if (! is_pt[i]) {
							est[2 * i] = where->attr = sim[i]; /* last simulation */
							push_point(data[i], where);
						} else
							data[i]->nsim_at_data++;
						where->X += data[i]->n_X;
					}
					where->X = X_ori; /* and put back */
					break;
				default: ErrMsg(ER_IMPOSVAL, "Getest(): wrong mode"); break;
			} /* switch (get_mode()) */
			for (i = 0; i < n_vars; i++)
				set_mv_double(&(est[2 * i + 1])); 
			for (i = 2 * n_vars; i < get_n_outputs(); i++)
				set_mv_double(&(est[i])); 
			/* print_sim(); */
			break;
		case NRS: 
			if (get_mode() != STRATIFY) {
				for (i = 0; i < n_vars; i++) {
					est[2 * i] = (double) data[i]->n_sel;
					est[2 * i + 1] = (double) data[i]->oct_filled;
				}
			} else {
				est[0] = (double) data[where->u.stratum]->n_sel;
				est[1] = (double) data[where->u.stratum]->oct_filled;
			}
			break;
		case SPREAD: 
			if (get_mode() == STRATIFY) {
				if (data[where->u.stratum]->n_sel > 0) {
					if (data[where->u.stratum]->what_is_u != U_ISDIST)
						ErrMsg(ER_IMPOSVAL, "getest() needs distances");
					n_sel = data[where->u.stratum]->n_sel;
					est[0] = sqrt(data[where->u.stratum]->sel[0]->u.dist2);
					est[1] = sqrt(data[where->u.stratum]->sel[n_sel - 1]->u.dist2);
				}
			} else {
				for (i = 0; i < n_vars; i++) {
					if (data[i]->n_sel > 0) {
						if (data[i]->what_is_u != U_ISDIST)
							ErrMsg(ER_IMPOSVAL, "getest() needs distances");
						n_sel = data[i]->n_sel;
						est[2 * i] = sqrt(data[i]->sel[0]->u.dist2);
						est[2 * i + 1] = sqrt(data[i]->sel[n_sel - 1]->u.dist2);
					}
				}
			}
			break;
		case LSEM:
			v = get_vgm(LTI(0,0));
			assert(v);
			v->id1 = v->id2 = v->id = 0;
			v->ev->evt = SEMIVARIOGRAM;
			if (data[0]->is_residual == 0)
				data[0]->is_residual = 1; /* !! */
			v->ev->recalc = 1;
			calc_variogram(v, NULL);
			if (gl_fit) { 
				v->ev->refit = 1;
				v->ev->fit = gl_fit;
				save_variogram_parameters(v);
				fit_variogram(v);
				/* write back locally fitted variogram model parameters: */
				for (i = j = 0; i < MIN(get_n_vars(), v->n_models); i++) {
					est[2 * i + 0] = v->part[i].sill;
					est[2 * i + 1] = v->part[i].range[0];
					j += 2;
				}
				if (j < get_n_vars())
					est[j] = 1.0 * v->fit_is_singular;
				reset_variogram_parameters(v);
				update_variogram(v); /* not sure if this is needed */
			} else {
				/* write back locally estimated sample variogram values: */
				for (i = 0; i < MIN(get_n_vars(), v->ev->n_est); i++) {
					est[2 * i + 0] = v->ev->gamma[i];
					est[2 * i + 1] = 1.0 * v->ev->nh[i];
				}
			}
			break;
		case NSP:  /* FALLTHROUGH: */
		default: 
			ErrMsg(ER_IMPOSVAL, "getest(): method not specified");
	} /* switch */
	efree(local_sim);
	efree(is_pt);
	return;
}

static void est_quantile_div(DATA *data, double *est, int div) {
	static double *list = NULL;
	static int i, size = 0;
	double mod = -9999.0;
	int n, run, longest_run = 0;

	if (data->n_sel > size) 
		list = (double *) erealloc(list, (size = data->n_sel) * sizeof(double));
	for (i = 0; i < data->n_sel; i++)
		list[i] = data->sel[i]->attr;
	qsort(list, (size_t) data->n_sel, sizeof(double),
		(int CDECL (*)(const void *, const void *)) d_cmp);
	if (div) { /* get diversity and modus in sel: */
		n = data->n_sel;
		for (i = 0; i < data->n_sel - 2; i += run) {
			run = 1;
			while (i + run < data->n_sel && list[i] == list[i + run]) {
				run++;
				n--;
			}
			if (run > longest_run) {
				longest_run = run;
				mod = list[i];
			}
		}
		est[0] = 1.0 * n;
		est[1] = mod;
	} else {
		if (data->n_sel < 2)
			return;
		if (gl_quantile == 0.5)
			est[0] = est[1] = est_quant(list, gl_quantile, data->n_sel);
		else {
			est[0] = est_quant(list, gl_quantile, data->n_sel);
			est[1] = est_quant(list, 1.0 - gl_quantile, data->n_sel);
		}
	}
	return;
}

static void est_skew_kurt(DATA *data, double *est) {
	static double *list = NULL;
	double mean, std, skewness = 0.0, kurtosis = 0.0, d;
	static int i, size = 0;

	if (data->n_sel <= 1) /* need at least 2 for variance */ 
		return;
	if (data->n_sel > size) 
		list = (double *) erealloc(list, (size = data->n_sel) * sizeof(double));
	for (i = 0; i < data->n_sel; i++)
		list[i] = data->sel[i]->attr;
	mean = sample_mean(list, data->n_sel);
	std = sample_std(list, mean, data->n_sel);
	for (i = 0; i < data->n_sel; i++) {
			d = (data->sel[i]->attr - mean)/std;
			skewness += pow(d, 3.0);
			kurtosis += pow(d, 4.0);
	}
	est[0] = skewness/data->n_sel;
	est[1] = kurtosis/data->n_sel;
}

static double inverse_dist(DATA *data, DPOINT *where, double idPow) {
	int i, j;
	double sumweights, sumvals, value, weight, dist2;
	static DATA *blockd = NULL;

	if (data->n_sel <= 0)
		ErrMsg(ER_IMPOSVAL, "zero neighbourhood in inverse_dist()");
	if (data->n_sel == 1)
		return data->sel[0]->attr;
	blockd = block_discr(blockd, get_block_p(), where);
	for (i = 0, value = 0.0; i < blockd->n_list; i++) {
		for (j = 0, sumweights = sumvals = 0.0; j < data->n_sel; j++) {
			dist2 = data->pp_norm2(data->sel[j], blockd->list[i]);
			if (dist2 == 0.0) {
				sumvals = data->sel[j]->attr;
				sumweights = 1.0;
				j = data->n_sel; /* break j loop */
			} else {
				if (idPow == 2.0)
					weight = 1.0 / dist2;
				else
					weight = pow(dist2, -0.5 * idPow);
				sumweights += weight;
				sumvals += weight * data->sel[j]->attr;
			}
		}
		value += blockd->list[i]->u.weight * sumvals/sumweights;
	}
	return value;
}

void save_variogram_parameters(VARIOGRAM *v) {
	int i;
	if (vgm_pars == NULL)
		vgm_pars = (double *) emalloc(3 * v->n_models * sizeof(double));
	for (i = 0; i < v->n_models; i++) {
		vgm_pars[3 * i + 0] = v->part[i].sill;
		vgm_pars[3 * i + 1] = v->part[i].range[0];
		vgm_pars[3 * i + 2] = v->part[i].range[1];
	}

}

void reset_variogram_parameters(VARIOGRAM *v) {
	int i;
	for (i = 0; i < v->n_models; i++) {
		v->part[i].sill = vgm_pars[3 * i + 0];
		v->part[i].range[0] = vgm_pars[3 * i + 1];
		v->part[i].range[1] = vgm_pars[3 * i + 2];
	}
	v->fit_is_singular = 0;
}

static double sample_mean(double *list, int n) {
	int i;
	double mn = 0.0;

	if (list == NULL)
		ErrMsg(ER_NULL, "sample_mean()");
	if (n == 0) 
		ErrMsg(ER_RANGE, "sample_mean(): no values");
	for (i = 0; i < n; i++)
		mn += list[i];
	return mn/(1.0 * n);
}

static double sample_var(double *list, double mean, int n) {
	int i;
	double var = 0.0;

	if (list == NULL)
		ErrMsg(ER_NULL, "sample_var()");
	if (n <= 1 || list == NULL) 
		ErrMsg(ER_RANGE, "sample_var(): <= 1 values");
	for (i = 0; i < n; i++)
		var += SQR(list[i] - mean);
	return (var / (n - 1.0));
}

static double sample_std(double *list, double mean, int n) {
	return sqrt(sample_var(list, mean, n));
}

static double est_quant(double *list, double p, int n) {
/*
* function returns the value of the p-th quantile
* of the ordered list *list, row exists from list[0]..list[length-1],
* p is in [0..1];
* a missing value is generated when the quantile lies not within
* the valid data range, and is undetermined therefore
*/
	double order, where;
	int below, above;

	if (n < 2)
		ErrMsg(ER_RANGE, "est_quant(): < 2 obs.");
	if (p < 0.0 || p > 1.0)
		ErrMsg(ER_RANGE, "can't calculate quantile outside [0,1]");
	order = p * (n - 1);
	/* order = n * (p * n)/(n + 1); */
	below = (int) floor(order); /* the index below order */
	if (below < 0) 
		return list[0];
	above = below + 1;			/* the index above order */
	if (above >= n) 
		return list[n - 1];
	where = order - below;
	return (1 - where) * list[below] + where * list[above]; 
}

int CDECL d_cmp(const double *a, const double *b) {
	if (*a < *b)
		return -1;
	if (*a > *b)
		return 1;
	return 0;
}
