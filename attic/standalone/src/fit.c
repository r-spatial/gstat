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
 * fit.c: fit variogram model to experimental variograms; gnuplot fit driver
 */
#include <stdio.h>
#include <stdlib.h> /* getenv() */
#include <string.h> /* strstr() */
#include <math.h> /* fabs(), sqrt() */

#include "matrix2.h"
#include "defs.h"
#include "defaults.h"
#include "userio.h"
#include "data.h"
#include "utils.h"
#include "read.h"
#include "debug.h"
#include "vario.h"
#include "sem.h"
#include "plot.h"
#include "glvars.h"
#include "reml.h"
#include "lm.h"
#include "fit.h"

#ifdef USING_R
void Rprintf(const char *, ...);
#endif

#define FIT_LOG     "fit.log"
#define NEARLY_ZERO     1e-30

static void wls_fit(VARIOGRAM *vp);
static double getSSErr(const VARIOGRAM *vp, PERM *p, LM *lm);
static void write_fx(FILE *, VARIOGRAM *v);
static void get_values(const char *fname, VARIOGRAM *v);
#ifndef USING_R
static void gnu_fit(VARIOGRAM *v);
#endif
static void correct_for_anisotropy(VARIOGRAM *v);

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

	if ((v->ev->fit == WLS_FIT_MOD || v->ev->fit == WLS_GNUFIT_MOD) &&
			!is_variogram(v))
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
#ifndef USING_R
		case WLS_GNUFIT: /* BREAKTHROUGH */
		case WLS_GNUFIT_MOD:
			gnu_fit(v);
			break;
#endif
		case MIVQUE_FIT:
			if (v->id1 != v->id2) 
				return 1;
			d = get_gstat_data();
			reml_sills(d[v->id1], v);
			break;
		default:
			ErrMsg(ER_IMPOSVAL, "fit_vgm(): value for fit not recognized");
		/*
		case LMC:
			d = get_gstat_data();
			fit_lmc(d, v, ...
		*/
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

	if (gl_cn_max < 0.0)
		lm->cn_max = 1.0/sqrt(MACHEPS);
	else
	 	lm->cn_max = gl_cn_max;

	/* oldSSErr = getSSErr(vp, p, lm); */
	do {
		print_progress(n_iter, gl_iter);
		if (DEBUG_VGMFIT) 
			printlog("%s: ", vp->descr);
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

		if (step < gl_fit_limit && step >= 0.0 && bounded == 0) 
			timetostop = 1;
		else if (n_iter > gl_iter)
			timetostop = 1;
		else 
			timetostop = 0;

	} while (! timetostop);

	print_progress(gl_iter, gl_iter);

	if (n_iter == gl_iter)
		pr_warning("No convergence after %d iterations", n_iter);

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

static int fit_GaussNewton(VARIOGRAM *vp, PERM *p, LM *lm, int iter,
		int *bounded) {
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
				lm->X->me[i][j] = (is_variogram(vp) ?
					UnitSemivariance(vp->part[model],x,y,z) :
					UnitCovariance(vp->part[model],x,y,z));
			} else {
				model = fit->ive[j] - vp->n_models;
				lm->X->me[i][j] = (is_variogram(vp) ?
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

#ifndef USING_R
	if (DEBUG_FIT) {
		printf("data: ");
		v_foutput(stdout, lm->y);
		printf("weights: ");
		v_foutput(stdout, lm->weights);
		printf("X: ");
		m_foutput(stdout, lm->X);
	}
#endif

	lm->has_intercept = 1; /* does not affect the fit */
	lm = calc_lm(lm); /* solve WLS eqs. for beta */

#ifndef USING_R
	if (DEBUG_FIT) {
		printf("beta: ");
		v_foutput(stdout, lm->beta);
	}
#endif

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

/* 
 * fit a variogram model to a sample variogram,
 * using the gnuplot fit command of gnuplot versions 3.6 and above
 */
#ifndef USING_R
static void gnu_fit(VARIOGRAM *v) {
	int i;
	FILE *f = NULL;
	char *fit_log = NULL, *cmd;

	cmd = NULL; /* to avoid compiler warning */
	if ((fit_log = getenv("FIT_LOG")) == NULL)
		fit_log = FIT_LOG;
	if (file_exists(fit_log))
		eremove(fit_log);
#ifdef HAVE_POPEN
	f = epopen(gl_gnuplot, "w");
#else
	f = efopen(GNUFIT_NAME, "w");
#endif
	if (is_variogram(v)) {
		fprintf(f, "Nug(a,x) = (x == 0.0 ? 0.0 : 1.0) # irresp. a\n");
		for (i = 2; v_models[i].name != NULL; i++)
			fprintf(f, "%s\n", v_models[i].v_gnuplot);
	} else {
		fprintf(f, "Nug(a,x) = (x == 0.0 ? 1.0 : 0.0) # irresp. a\n");
		for (i = 2; v_models[i].name != NULL; i++)
			fprintf(f, "%s\n", v_models[i].c_gnuplot);
	}
	write_fx(f, v);
	fprint_sample_vgm(f, v->ev);
	fprintf(f, "e\n");
#ifdef HAVE_POPEN
	epclose(f);
#else
	fclose(f);
	cmd = (char *) emalloc(strlen(gl_gnuplot) + strlen(GNUFIT_NAME) + 2);
	sprintf(cmd, "%s %s", gl_gnuplot, GNUFIT_NAME);
	esystem(cmd); /* let gnuplot do the fit */
	efree(cmd);
#endif
	get_values(fit_log, v);
	update_variogram(v);
	return;
}

static void write_fx(FILE *f, VARIOGRAM *v) {
/*
 * write:
 * c1 = ..; a1 = ..; ...
 * f(x) = c1 * Nug(0) + c2 * Sph(a2) + ...
 * fit f(x) '-' via c1, c2, ...
 */
 	int i, first = 1;
/*
 * set parateters, c0=..; a0=..;
 */
 	for (i = 0; i < v->n_models; i++) {
 		if (v->part[i].fit_sill)
 			fprintf(f, "c%d = %g; ", i, v->part[i].sill);
 		if (v->part[i].fit_range)
 			fprintf(f, "a%d = %g; ", i, v->part[i].range[0]);
 	}
/*
 * f(x) = ...
 */
 	fprintf(f, "\nf(x) = ");
	fprint_gnuplot_model(f, v, 1);
/*
 * fit f(x) 'file' using a:b:c # c are `uncertainties'
 * after a bit trial and error, they seem to be error standard deviations.
 */
	if (v->ev->fit == WLS_GNUFIT_MOD)
		fprintf(f, "\nfit f(x) '-' using %s via",
			v->ev->cloud ? "3:4:(sqrt(f($3)**2))" : "4:5:(sqrt((f($4)**2)/$3))");
	else
		fprintf(f, "\nfit f(x) '-' using %s via",
			v->ev->cloud ? "3:4" : "4:5:(sqrt(1.0/$3))" );
/*
 * via
 */
 	for (i = 0; i < v->n_models; i++) {
 		if (v->part[i].fit_sill) {
 			fprintf(f, "%sc%d", first ? " " : ", ", i);
 			first = 0;
 		}
 		if (v->part[i].fit_range) {
 			fprintf(f, "%sa%d", first ? " " : ", ", i);
 			first = 0;
 		}
 	}
 	fprintf(f, "\n");
	return;
}

static void get_values(const char *fname, VARIOGRAM *v) {
	FILE *f = NULL;
	char *s = NULL, *cp;
	int size = 0, nr = 0;

	if (! file_exists(fname)) {
		pr_warning("can't find %s -- get gnuplot 3.7 from www.gnuplot.info", 
			fname);
		return;
	}
	f = efopen(fname, "r");
	do {
		if (get_line(&s, &size, f) == NULL) {
			pr_warning("error while reading %s", fname);
			efclose(f);
			return;
		} 
		if (strstr(s, "BREAK:          Singular matrix")) {
			pr_warning("error during variogram fit (singular model)");
			efree(s);
			efclose(f);
			return;
		}

	} while (strstr(s, "Final set of parameters") == NULL);
	get_line(&s, &size, f); /* line with the many ==== */
	get_line(&s, &size, f); /* empty line */
	while (get_line(&s, &size, f) != NULL) { /* values */
		if (s[0] == 'c' || s[0] == 'a') {
			cp = s;
			cp++;
			cp = strtok(cp, " =");
			if (read_int(cp, &nr) || nr < 0 || nr >= v->n_models)
				ErrMsg(ER_IMPOSVAL, fname);
			cp = strtok(NULL, " =");
			if (cp == NULL)
				ErrMsg(ER_IMPOSVAL, fname);
			if (s[0] == 'c')
				read_double(cp, &(v->part[nr].sill));
			else
				read_double(cp, &(v->part[nr].range[0]));
		} else
			break; /* while-loop */
	}
	efclose(f);
	efree(s);
	correct_for_anisotropy(v);
	return;
}
#endif /* USING_R */

static void correct_for_anisotropy(VARIOGRAM *v) {
/*
 * Ok, so now let's get to the hard part. What happened?
 * in case of anisotropy, non-fitted values are correctly adjusted to
 * their directional value in fprint_gnuplot_model().
 * If however a value is fitted, the real fitted sill or range should
 * be the corresponding range/sill in the main direction. So, we have
 * to transfer those.
 */
 	int i;
 
 	for (i = 0; i < v->n_models; i++) {
 		if (v->part[i].fit_range && v->part[i].tm_range)
 			v->part[i].range[0] /= relative_norm(v->part[i].tm_range,
 				v->ev->direction.x, v->ev->direction.y, v->ev->direction.z);
 	}
	return;
}
