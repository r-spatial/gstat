/*
 * fit.c: fit variogram model to experimental variograms; gnuplot fit driver
 */
#include <stdio.h>
#include <stdlib.h> /* getenv() */
#include <string.h> /* strstr() */
#include <math.h> /* fabs(), sqrt() */

#include <R.h> /* Rprintf() */

#include "mtrx.h"
#include "defs.h"
#include "defaults.h"
#include "userio.h"
#include "data.h"
#include "utils.h"
#include "debug.h"
#include "vario.h"
#include "sem.h"
#include "glvars.h"
#include "reml.h"
#include "lm.h"
#include "fit.h"

#define NEARLY_ZERO     1.0e-30

static void wls_fit(VARIOGRAM *vp);
static double getSSErr(const VARIOGRAM *vp, PERM *p, LM *lm);

static int fill_weights(const VARIOGRAM *vp, PERM *p, LM *lm);
static int fit_GaussNewton(VARIOGRAM *vp, PERM *p, LM *lm, 
	int iter, int *bounded);

int fit_variogram(VARIOGRAM *v) {
	DATA **d = NULL;
	int i = 0;
	long n = 0;

	if (v->ev->refit == 0)
		return 0;
	if (v->ev->fit != NO_FIT && v->ev->fit != MIVQUE_FIT) {
		if (! v->ev->cloud) {
			while (n == 0 && i < v->ev->n_est)
				n += v->ev->nh[i++]; /* check if estimates exist */
			if (n == 0) /* bad luck */
				return 1;
		} else if (v->ev->n_est == 0)
			return 1;
	}
	if (v->ev->fit != NO_FIT) {
		if (v->ev->map)
			ErrMsg(ER_IMPOSVAL, "cannot fit model to variogram map");
		for (i = 0; i < v->n_models; i++)
			if (v->part[i].sill == 0.0 && v->part[i].fit_sill != 0)
				v->part[i].sill = 1.0; /* avoid lot'o trouble */
	}

	if (v->ev->fit == WLS_FIT_MOD && !is_variogram(v))
		pr_warning("this fit method is not recommended for covariograms");
	v->ev->direction.x = sin(gl_alpha * PI / 180.0) * cos(gl_beta * PI / 180.0);
	v->ev->direction.y = cos(gl_alpha * PI / 180.0) * cos(gl_beta * PI / 180.0);
	v->ev->direction.z = sin(gl_beta * PI / 180.0);

	switch (v->ev->fit) {
		case NO_FIT: 
			break;
		case OLS_FIT: /* BREAKTHROUGH: */
		case WLS_FIT: /* BREAKTHROUGH: */
		case WLS_FIT_MOD:
		case WLS_NHH:
			wls_fit(v);
			break;
		case MIVQUE_FIT:
			if (v->id1 != v->id2) 
				return 1;
			d = get_gstat_data();
			reml_sills(d[v->id1], v);
			break;
		default:
			Rprintf("%d\n", v->ev->fit);
			ErrMsg(ER_IMPOSVAL, "fit_vgm(): value for fit not recognized or not implemented");
		/* no default: force compile warning on missing option! */
	}
	return 0;
}

static void wls_fit(VARIOGRAM *vp) {
/*
 * non-linear iterative reweighted least squares fitting of variogram model to
 * sample variogram (..covariogram model to sample covariogram, cross, etc.)
 * all information necessary is contained in *vp.
 *
 * uses Marquardt-Levenberg algorithm;
 * the implementation follows gnuplot's fit.c
 */
	static PERM *p = PNULL;
	int i, j, n_iter = 0, bounded = 0, timetostop;
	double SSErr, oldSSErr = DBL_MAX, step;
	LM *lm;

	p = px_resize(p, vp->ev->n_est);
	if (! vp->ev->cloud) {
		for (i = j = 0; i < (vp->ev->zero == ZERO_AVOID ?
				vp->ev->n_est-1 : vp->ev->n_est); i++) {
			if (vp->ev->nh[i] > 0)
				p->pe[j++] = i;
		}
		p->size = j;
	} 
	lm = init_lm(NULL);

	/* oldSSErr = getSSErr(vp, p, lm); */
	do {
		print_progress(n_iter, gl_iter);
		/* if (DEBUG_VGMFIT) 
			printlog("%s: ", vp->descr); */
		if ((vp->fit_is_singular = fit_GaussNewton(vp, p, lm, n_iter, &bounded))) {
			pr_warning("singular model in variogram fit");
			print_progress(gl_iter, gl_iter);
			vp->SSErr = getSSErr(vp, p, lm);
			return;
		} 
		update_variogram(vp);

		SSErr = getSSErr(vp, p, lm);
		/* we can't use lm->SSErr here since that's only in the
		X-filled-with-derivatives, not the true residuals */

		step = oldSSErr - SSErr;
		if (SSErr > gl_zero)
			step /= SSErr;

		n_iter++;

		if (DEBUG_VGMFIT)
			printlog("after it. %d: SSErr %g->%g, step=%g (fit_limit %g%s)\n",
					n_iter, oldSSErr, SSErr, step, gl_fit_limit, 
					bounded ? "; bounded" : "");

		oldSSErr = SSErr;
		timetostop = (step < gl_fit_limit && step >= 0.0 && bounded == 0) || n_iter == gl_iter;
	} while (! timetostop);

	print_progress(gl_iter, gl_iter);

	if (n_iter == gl_iter)
		pr_warning("No convergence after %d iterations: try different initial values?", n_iter);

	if (DEBUG_VGMFIT) {
		printlog("# iterations: %d, SSErr %g, last step %g", n_iter, SSErr, step);
		if (step < 0.0)
			printlog(", last step was in the wrong direction.\n");
		else
			printlog("\n");
	}

	free_lm(lm);
	vp->SSErr = SSErr;
	return;
} /* wls_fit */

static double getSSErr(const VARIOGRAM *vp, PERM *p, LM *lm) {

	int i;
	double x, y, z, dz, SSErr;

	/*
	if (fill_weights(vp, p, lm))
		return 0.0;
	*/
	for (i = 0, SSErr = 0.0; i < p->size; i++) {
		x = vp->ev->direction.x * vp->ev->dist[p->pe[i]];
		y = vp->ev->direction.y * vp->ev->dist[p->pe[i]];
		z = vp->ev->direction.z * vp->ev->dist[p->pe[i]];
		/* fill y with current residuals: */
		dz = vp->ev->gamma[p->pe[i]] - (is_variogram(vp) ? 
			get_semivariance(vp, x, y, z) :
			get_covariance(vp, x, y, z));
		if (lm->weights != NULL)
			SSErr += dz * dz * lm->weights->ve[i];
		else
			SSErr += dz * dz;
	}
	return SSErr;
}

static int fit_GaussNewton(VARIOGRAM *vp, PERM *p, LM *lm, int iter, int *bounded) {
	double s = 0.0, x, y, z;
	int i, j, n_fit, model, fit_ranges = 0;
	IVEC *fit = NULL;
	VEC *start = NULL;

	if (p->size == 0)
		return 1;

	fit = iv_resize(fit, 2 * vp->n_models);
	/* index fit parameters: parameter fit->ive[j] corresponds to model i */
	for (i = n_fit = 0; i < vp->n_models; i++) {
		if (vp->part[i].fit_sill)
			fit->ive[n_fit++] = i;
		if (vp->part[i].fit_range) {
			fit->ive[n_fit++] = i + vp->n_models; /* large -->> ranges */
			fit_ranges = 1;
		}
	}
	if (n_fit == 0) {
		iv_free(fit);
		return 0;
	}

	fit = iv_resize(fit, n_fit); /* shrink to fit */
	lm->X = m_resize(lm->X, p->size, n_fit);
	lm->y = v_resize(lm->y, p->size);
	start = v_resize(start, n_fit);

	for (i = 0; i < n_fit; i++) {
		if (fit->ive[i] < vp->n_models) {
			model = fit->ive[i];
			start->ve[i] = vp->part[model].sill;
		} else {
			model = fit->ive[i] - vp->n_models;
			start->ve[i] = vp->part[model].range[0];
		}
	}

	for (i = 0; i < p->size; i++) {
		x = vp->ev->direction.x * vp->ev->dist[p->pe[i]];
		y = vp->ev->direction.y * vp->ev->dist[p->pe[i]];
		z = vp->ev->direction.z * vp->ev->dist[p->pe[i]];
		/* fill y with current residuals: */
		if (is_variogram(vp))
			s = get_semivariance(vp, x, y, z);
		else
			s = get_covariance(vp, x, y, z);
		lm->y->ve[i] = vp->ev->gamma[p->pe[i]] - s;
		/* fill X: */
		for (j = 0; j < n_fit; j++) { /* cols */
			if (fit->ive[j] < vp->n_models) {
				model = fit->ive[j];
				ME(lm->X, i, j) = (is_variogram(vp) ?
					UnitSemivariance(vp->part[model],x,y,z) :
					UnitCovariance(vp->part[model],x,y,z));
			} else {
				model = fit->ive[j] - vp->n_models;
				ME(lm->X, i, j) = (is_variogram(vp) ?
					da_Semivariance(vp->part[model],x,y,z) :
					-da_Semivariance(vp->part[model],x,y,z));
			}
		}
	}

	if (iter == 0 && fill_weights(vp, p, lm)) {
		iv_free(fit);
		v_free(start);
		return 1;
	}

	lm->has_intercept = 1; /* does not affect the fit */
	lm = calc_lm(lm); /* solve WLS eqs. for beta */

	if (DEBUG_FIT) {
		Rprintf("beta: ");
		v_logoutput(lm->beta);
	}

	if (lm->is_singular) {
		iv_free(fit);
		v_free(start);
		return 1;
	}

	if (fit_ranges) {
		s = v_norm2(lm->beta) / v_norm2(start);
		if (s > 0.2) {
			/* don't allow steps > 20% ---- */
			sv_mlt(0.2 / s, lm->beta, lm->beta); 
			*bounded = 1;
		} else
			*bounded = 0; /* a `free', voluntary step */
	} else /* we're basically doing linear regression here: */
		*bounded = 0;

	for (i = 0; i < n_fit; i++) {
		if (fit->ive[i] < vp->n_models) {
			model = fit->ive[i];
			vp->part[model].sill = start->ve[i] + lm->beta->ve[i];
		} else {
			model = fit->ive[i] - vp->n_models;;
			vp->part[model].range[0] = start->ve[i] + lm->beta->ve[i];
		}
	}
	iv_free(fit);
	v_free(start);
	return 0;
}

static int fill_weights(const VARIOGRAM *vp, PERM *p, LM *lm) {
	double x, y, z, s;
	int i, retval = 0;

	if (vp->ev->fit == OLS_FIT || (vp->ev->cloud && vp->ev->fit == WLS_FIT)) {
		if (lm->weights != NULL)
			v_free(lm->weights);
		lm->weights = NULL;
		return 0;
	}
	lm->weights = v_resize(lm->weights, p->size);
	for (i = 0; i < p->size; i++) {
		if (vp->ev->fit == WLS_NHH)
			s = vp->ev->dist[p->pe[i]];
		else {
			x = vp->ev->direction.x * vp->ev->dist[p->pe[i]];
			y = vp->ev->direction.y * vp->ev->dist[p->pe[i]];
			z = vp->ev->direction.z * vp->ev->dist[p->pe[i]];
			s = get_semivariance(vp, x, y, z);
		}
		if (vp->ev->cloud) {
			lm->weights->ve[i] = 1.0;
			if (fabs(s) > NEARLY_ZERO)
				lm->weights->ve[i] /= s * s;
			else {
				pr_warning("infinite weight during fit");
				retval = 1;
				break;
			}
		} else { /* no cloud: */
			lm->weights->ve[i] = vp->ev->nh[p->pe[i]];
			if (vp->ev->fit != WLS_FIT) {
				if (fabs(s) > NEARLY_ZERO)
					lm->weights->ve[i] /= (s * s);
				else {
					pr_warning("infinite weight during fit");
					retval = 1;
					break;
				}
			}
		} /* else cloud */
	} /* for i */
	return retval;
}
