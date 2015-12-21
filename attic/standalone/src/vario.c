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
 * vario.c: basic variogram model functions (init, print, update, etc.)
 */

#include <stdio.h>
#include <stdlib.h> /* getenv() */
#include <ctype.h> /* toupper() */
#include <float.h> /* DBL_MIN */
#include <math.h>
#include <string.h>

#include "defs.h"
#include "matrix2.h"

#include "userio.h"
#include "data.h"
#include "utils.h"
#include "debug.h"
#include "vario.h"
#include "vario_fn.h"
#include "glvars.h"
#include "lex.h" /* read_variogram() */
#include "read.h" /* for vario() only */
#include "lm.h"

static int is_valid_cs(const VARIOGRAM *aa, const VARIOGRAM *bb,
	const VARIOGRAM *ab);
static int is_posdef(MAT *m);
static void strcat_tm(char *cp, ANIS_TM *tm);
static ANIS_TM *get_tm(double anis[5]);
static void init_variogram_part(VGM_MODEL *v);

const V_MODEL v_models[] = { /* the variogram model ``data base'': */
	{	 NOT_SP, "Nsp", "Nsp (not specified)",  /* DON'T CHANGE THIS ONE!! */
		"# NSP: should never occur",
		"# NSP: should never occur",
		NULL, NULL },
	{	NUGGET, "Nug", "Nug (nugget)", 
		"Nug(a,x) = 1 # I don't let gnuplot draw what happens at x=0",
		"Nug(a,x) = 0 # I don't let gnuplot draw what happens at x=0",
		fn_nugget, da_is_zero },
	{	EXPONENTIAL, "Exp", "Exp (exponential)", 
		"Exp(a,x) = 1 - exp(-x/a)",
		"Exp(a,x) = exp(-x/a)",
		fn_exponential, da_fn_exponential },
	{ 	SPHERICAL, "Sph", "Sph (spherical)", 
		"Sph(a,x) = (x < a ? (1.5 * x/a) - 0.5*((x/a)**3) : 1)",
		"Sph(a,x) = (x < a ? (1 - ((1.5 * x/a) - 0.5*((x/a)**3))) : 0)",
		fn_spherical, da_fn_spherical },
	{	GAUSSIAN, "Gau", "Gau (gaussian)", 
		"Gau(a,x) = 1 - exp(-((x/a)**2))",
		"Gau(a,x) = exp(-((x/a)**2))",
		fn_gaussian, da_fn_gaussian },
	{	EXCLASS, "Exc", "Exclass (Exponential class)", 
		"# Exponential class model not supported by gnuplot",
		"# Exponential class model not supported by gnuplot",
		fn_exclass, NULL },
#ifdef USING_R
	{	MATERN, "Mat", "Mat (Matern)", 
		"# Matern model not supported by gnuplot",
		"# Matern model not supported by gnuplot",
		fn_matern, NULL },
	{	STEIN, "Ste", "Mat (Matern, M. Stein's parameterization)", 
		"# Matern model not supported by gnuplot",
		"# Matern model not supported by gnuplot",
		fn_matern2, NULL },
#endif
	{ 	CIRCULAR, "Cir", "Cir (circular)", 
		"Cir(a,x) = (x < a ? ((2*x)/(pi*a))*sqrt(1-(x/a)**2)+(2/pi)*asin(x/a) : 1)",
		"Cir(a,x) = (x < a ? (1-(((2*x)/(pi*a))*sqrt(1-(x/a)**2)+(2/pi)*asin(x/a))) : 0)",
		fn_circular, NULL },
	{	LINEAR, "Lin", "Lin (linear)",  /* one-parameter (a = 0), or two-parameter with sill */
		"Lin(a,x) = (a > 0 ? (x < a ? x/a : 1) : x)",
		"Lin(a,x) = (a > 0 ? (x < a ? (1 - x/a) : 0) : 1 - x)",
		fn_linear, da_fn_linear },
	{	BESSEL, "Bes", "Bes (bessel)", 
		"# Bessel model not suppurted by gnuplot",
		"# Bessel model not suppurted by gnuplot",
		fn_bessel, NULL },
	{	PENTASPHERICAL, "Pen", "Pen (pentaspherical)", 
		"Pen(a,x) = (x < a ? ((15.0/8.0)*(x/a)+(-5.0/4.0)*((x/a)**3)+\
(3.0/8.0)*((x/a)**5)) : 1)",
		"Pen(a,x) = (x < a ? (1 - ((15.0/8.0)*(x/a)+(-5.0/4.0)*((x/a)**3)+\
(3.0/8.0)*((x/a)**5))) : 0)",
		fn_pentaspherical, da_fn_pentaspherical },
	{	PERIODIC, "Per", "Per (periodic)", 
		"Per(a,x) = 1 - cos(2*pi*x/a)",
		"Per(a,x) = cos(2*pi*x/a)",
		fn_periodic, da_fn_periodic },
	{ 	HOLE, "Hol", "Hol (hole)",
		"Hol(a,x) = 1 - sin(x/a)/(x/a)",
		"Hol(a,x) = sin(x/a)/(x/a)",
		fn_hole, da_fn_hole },
	{	LOGARITHMIC, "Log", "Log (logarithmic)", 
		"Log(a,x) = log(x + a)",
		"Log(a,x) = 1 - log(x + a)",
		fn_logarithmic, da_fn_logarithmic },
	{	POWER, "Pow", "Pow (power)",
		"Pow(a,x) = x ** a",
		"Pow(a,x) = 1 - x ** a",
		fn_power, da_fn_power },
	{	SPLINE, "Spl", "Spl (spline)", 
		/* Wackernagel 2nd ed., p. 225 -- not working yet */
		"Spl(a,x) = x == 0 ? 0 : x * x * log(x)",
		"Spl(a,x) = x == 0 ? 1 : 1 - x * x * log(x)",
		fn_spline, NULL },
	{	LEGENDRE, "Leg", "Leg (Legendre)", 
		"",
		"",
		fn_legendre, NULL },
	{	MERROR, "Err", "Err (Measurement error)", 
		"Err(a,x) = 1 # I don't let gnuplot draw what happens at x=0",
		"Err(a,x) = 0 # I don't let gnuplot draw what happens at x=0",
		fn_nugget, da_is_zero },
	/* the folowing two should always be the last ``valid'' one: */
	{	INTERCEPT, "Int", "Int (Intercept)",
		"Int(a,x)   = 1",
		"Int(a,x)   = 1",
		fn_intercept, da_is_zero },
	{	NOT_SP, NULL, NULL, NULL, NULL, NULL, NULL } /* THIS SHOULD BE LAST */
};

const char *vgm_type_str[] = { 
	"not specified",
	"semivariogram",
	"cross variogram",
	"covariogram",
	"cross covariogram"
};

VARIOGRAM *init_variogram(VARIOGRAM *v) {
/*
 * initializes one variogram structure
 * if v is NULL, memory is allocated for the structure
 */
	int i;

	if (v == NULL)
		v = (VARIOGRAM *) emalloc(sizeof(VARIOGRAM));

	v->id = v->id1 = v->id2 = -1;
	v->n_models = 0;
	v->is_valid_covariance = 1;
	v->isotropic = 1;
	v->n_fit = 0;
	v->fit_is_singular = 0;
	v->descr = v->fname = v->fname2 = (char *) NULL;
	v->max_range = (double) DBL_MIN;
	v->sum_sills = 0.0;
	v->measurement_error = 0.0;
	v->max_val = 0.0;
	v->min_val = 0.0;
	vgm_init_block_values(v);
	v->part = (VGM_MODEL *) emalloc(INIT_N_VGMM * sizeof(VGM_MODEL));
	v->table = NULL;
	for (i = 0; i < INIT_N_VGMM; i++)
		init_variogram_part(&(v->part[i]));
	v->max_n_models = INIT_N_VGMM;
	v->SSErr = 0.0;
	v->ev = init_ev();

	return v;
}

void vgm_init_block_values(VARIOGRAM *v) {
	v->block_semivariance_set = 0;
	v->block_covariance_set = 0;
	v->block_covariance = -999999.0;
	v->block_semivariance = -999999.0;
}

static void init_variogram_part(VGM_MODEL *p) {
	int i;

	p->sill = 0.0;
	for (i = 0; i < NRANGEPARS; i++)
		set_mv_double(&(p->range[i])); /* trigger errors if misused */
	p->model = NOT_SP;
	p->fit_sill = p->fit_range = 1;
	p->fnct = p->da_fnct = NULL;
	p->tm_range = NULL;
	p->id = -1;
}

SAMPLE_VGM *init_ev(void) {
	SAMPLE_VGM *ev = NULL;

	ev = (SAMPLE_VGM *) emalloc(sizeof(SAMPLE_VGM));
	set_mv_double(&(ev->cutoff));
	set_mv_double(&(ev->iwidth));
	ev->gamma = NULL;
	ev->dist = NULL;
	ev->nh = NULL;
	ev->pairs = NULL;
	ev->n_max = 0;
	ev->n_est = 0;
	ev->zero = ZERO_DEFAULT;
	ev->plot_numbers = 1;
	ev->is_directional = 0;
	ev->evt = NOTSPECIFIED;
	ev->fit = NO_FIT;
	ev->recalc = 1;
	ev->refit = 1;
	ev->pseudo = 0;
	ev->is_asym = -1;
	ev->map = NULL;
	ev->S_grid = NULL;
	ev->direction.x = 1.0;
	ev->direction.y = ev->direction.z = 0.0;
	return ev;
}


void free_variogram(VARIOGRAM *v) {
	int i;

	assert(v != NULL);
	if (v->ev) {
		if (v->ev->n_max > 0) {
			efree(v->ev->gamma);
			efree(v->ev->dist);
			efree(v->ev->nh);
			if (v->ev->pairs)
				efree(v->ev->pairs);
		}
		efree(v->ev);
	}
	for (i = 0; i < v->max_n_models; i++)
		if (v->part[i].tm_range != NULL)
			efree(v->part[i].tm_range);
		
	if (v->descr != NULL)
		efree(v->descr);

	efree(v->part);
	if (v->table) {
		efree(v->table->values);
		efree(v->table);
	}
	efree(v);
}

void logprint_variogram(const VARIOGRAM *v, int verbose) {
	printlog("%s", sprint_variogram(v, verbose));
}

#ifndef USING_R
void fprint_variogram(FILE *f, const VARIOGRAM *v, int verbose) {
	fprintf(f, "%s", sprint_variogram(v, verbose));
}
#endif

const char *sprint_variogram(const VARIOGRAM *v, int verbose) {
/* prints contents of VARIOGRAM v on string */
	static char tmp[ERROR_BUFFER_SIZE], s[ERROR_BUFFER_SIZE];
	int i, j, k;

	tmp[0] = s[0] = '\0';

	if (v->id1 < 0 && v->id2 < 0)
		return s; /* never set */

	if ((v->descr == NULL || v->n_models == 0) && (v->fname == NULL))
		return s; /* nothing to print */

	if (v->id1 == v->id2)
		sprintf(s, "variogram(%s)", name_identifier(v->id1)); 
	else
		sprintf(s, "variogram(%s,%s)", name_identifier(v->id1),
			name_identifier(v->id2)); 

	if (v->fname) {
		sprintf(tmp, ": '%s'", v->fname);
		strcat(s, tmp);
	}

	if (v->descr && v->n_models > 0) {
		sprintf(tmp, ": %s", v->descr);
		strcat(s, tmp);
	}

	strcat(s, ";\n");

	if (verbose == 0)
		return s;

	for (i = 0; i < v->n_models; i++) {
		sprintf(tmp, "# model: %d type: %s sill: %g range: %g\n", 
			i, v_models[v->part[i].model].name_long, 
			v->part[i].sill, v->part[i].range[0]);
		strcat(s, tmp);
		if (v->part[i].tm_range != NULL) {
			sprintf(tmp, "# range anisotropy, rotation matrix:\n");
			strcat(s, tmp);
			for (j = 0; j < 3; j++) {
				for (k = 0; k < 3; k++) {
					sprintf(tmp, "%s%8.4f", k == 0 ? "# " : " ",
						v->part[i].tm_range->tm[j][k]);
					strcat(s, tmp);
				}
				strcat(s, "\n");
			}
		}
	}
	sprintf(tmp, "# sum sills %g, max %g, min %g, flat at distance %g\n",
		v->sum_sills, v->max_val, v->min_val, v->max_range);
	strcat(s, tmp);
	return s;
}

void update_variogram(VARIOGRAM *vp) {
/*
 * update min/max, n_fit, descr
 * assumes that models are not changed: they can only be changed through
 * read_variogram();
 */
	char s[LENGTH_OF_MODEL], *cp;
	VGM_MODEL *p;
	int i;

	vp->descr = (char *) erealloc(vp->descr,
		vp->max_n_models * LENGTH_OF_MODEL * sizeof(char));
	cp = vp->descr;
	*cp = '\0';
	/* update sum_sills: */
	vp->sum_sills = vp->min_val = vp->max_val = 0.0;
	vp->n_fit = 0;
	vp->max_range = DBL_MIN;
	for (i = 0; i < vp->n_models; i++) {
		p = &(vp->part[i]);
		vp->sum_sills += p->sill;
		if (p->sill < 0.0)
			vp->min_val += p->sill;
		else
			vp->max_val += p->sill;
		vp->max_range = MAX(p->range[0], vp->max_range);

		if (p->model == BESSEL || p->model == GAUSSIAN ||
				p->model == EXPONENTIAL || p->model == LOGARITHMIC ||
				p->model == POWER || p->model == PERIODIC ||
				p->model == EXCLASS || p->model == LEGENDRE ||
				p->model == HOLE || /* more??? */
#ifdef USING_R
				p->model == MATERN ||
				p->model == STEIN ||
#endif
				(p->model == LINEAR && p->range[0] == 0)) 
					/* sill is reached asymptotically or oscillates */
			vp->max_range = DBL_MAX;
		else  /* transitive model: */
			vp->max_range = MAX(p->range[0], vp->max_range);

		if (p->fit_sill == 0)
			strcat(cp, "@ ");
		sprintf(s, gl_format, i == 0 ? p->sill : fabs(p->sill));
		strcat(cp, s);
		strcat(cp, " ");
		sprintf(s, "%s(", v_models[p->model].name);
		strcat(cp, s);
		if ((p->model == LINEAR && p->range[0] == 0.0) || p->model == NUGGET || p->model == INTERCEPT)
			p->fit_range = 0; /* 1 would lead to singularity */
		else if (p->fit_range == 0)
			strcat(cp, "@ ");
		sprintf(s, gl_format, p->range[0]);
		strcat(cp, s);
		if (p->tm_range != NULL) 
			strcat_tm(cp, p->tm_range);
		strcat(cp, ")");
		if (i != vp->n_models - 1)
			strcat(cp, vp->part[i+1].sill < 0.0 ? " - " : " + ");
		if (p->model == LOGARITHMIC || p->model == POWER || p->model == INTERCEPT
				|| (p->model == LINEAR && p->range[0] == 0))
		 	vp->is_valid_covariance = 0;
		if (p->fit_sill)
			vp->n_fit++;
		if (p->fit_range)
			vp->n_fit++;
		if (p->model == MERROR)
			vp->measurement_error += p->sill;
	}
	if (vp->table != NULL) {
		vp->sum_sills = vp->table->values[0];
		vp->max_val = vp->table->values[0];
		vp->min_val = vp->table->values[0];
		for (i = 1; i < vp->table->n; i++) {
			vp->max_val = MAX(vp->max_val, vp->table->values[i]);
			vp->min_val = MIN(vp->min_val, vp->table->values[i]);
		}
	}
	return;
}

static void strcat_tm(char *cp, ANIS_TM *tm) {
	char s[100];

	strcat(cp, ",");
	sprintf(s, gl_format, tm->angle[0]);
	if (TM_IS3D(tm)) {
		strcat(cp, s); strcat(cp, ",");
		sprintf(s, gl_format, tm->angle[1]);
		strcat(cp, s); strcat(cp, ",");
		sprintf(s, gl_format, tm->angle[2]);
	}
	strcat(cp, s); strcat(cp, ",");
	sprintf(s, gl_format, tm->ratio[0]);
	if (TM_IS3D(tm)) {
		strcat(cp, s); strcat(cp, ",");
		sprintf(s, gl_format, tm->ratio[1]);
	}
	strcat(cp, s);
}

double get_max_sill(int n) {
	int i, j;
	VARIOGRAM *vp;
	static double max_sill;

	vp = get_vgm(0);
	max_sill = vp->max_val;
	for (i = 0; i < n; i++) {
		for (j = 0; j <= i; j++) {
			vp = get_vgm(LTI(i,j));
			max_sill = MAX(max_sill, vp->max_val);
		}
	}
	return max_sill;
}

double get_semivariance(const VARIOGRAM *vp, double dx, double dy, double dz) {
/* returns gamma(dx,dy,dz) for variogram v: gamma(h) = cov(0) - cov(h) */

	int i;
	double sv = 0.0, dist = 0.0;

	if (vp->table != NULL)
		return(SEM_TABLE_VALUE(vp->table, 
				transform_norm(vp->table->tm_range, dx, dy, dz)));
	if (! vp->isotropic) {
		for (i = 0; i < vp->n_models; i++)
			sv += vp->part[i].sill * vp->part[i].fnct(
				transform_norm(vp->part[i].tm_range, dx, dy, dz),
				vp->part[i].range);
	} else {
		dist = transform_norm(NULL, dx, dy, dz);
		if (dist > vp->max_range)
			return vp->sum_sills;
		for (i = 0; i < vp->n_models; i++)
			sv += vp->part[i].sill * vp->part[i].fnct(dist, vp->part[i].range);
	}
	return sv;
}

double get_covariance(const VARIOGRAM *vp, double dx, double dy, double dz) {
/* returns cov(dx,dy,dz) for variogram v */

	int i;
	static int warning = 0;
	double ctmp = 0.0, dist;

	if (vp == NULL) {
		warning = 0;
		return 0.0;
	}

	if (! vp->is_valid_covariance && !warning) {
		pr_warning(
			"%s: non-transitive variogram model not allowed as covariance function",
			vp->descr);
		warning = 1;
	}
	if (!vp->is_valid_covariance && !DEBUG_FORCE)
		ErrMsg(ER_IMPOSVAL, "covariance from non-transitive variogram not allowed ");
	if (vp->table != NULL)
		return(COV_TABLE_VALUE(vp->table, 
				transform_norm(vp->table->tm_range, dx, dy, dz)));
	if (! vp->isotropic) {
		for (i = 0; i < vp->n_models; i++)
			ctmp += vp->part[i].sill * (1.0 - vp->part[i].fnct(
				transform_norm(vp->part[i].tm_range, dx, dy, dz),
				vp->part[i].range));
	} else {
		dist = transform_norm(NULL, dx, dy, dz);
		for (i = 0; i < vp->n_models; i++)
			ctmp += vp->part[i].sill * (1.0 - vp->part[i].fnct(dist,
				vp->part[i].range));
	}
	return ctmp;
}

static int is_valid_cs(const VARIOGRAM *aa, const VARIOGRAM *bb,
	const VARIOGRAM *ab)
/*
 * Purpose       : Check Cauchy-Schwartz inequality on cross/variograms    
 * Created by    : Edzer J. Pebesma                                       
 * Date          : may 6th 1992                                          
 * Prerequisites :                                                     
 * Returns       : return nonzero if |g_ab(h)| > sqrt(g_aa(h)g_bb(h))
 * Side effects  : none                                             
 */
{
	int i, check_failed = 0;
	double maxrange = 0, dist, dx, dy, dz;

	for (i = 0; i < aa->n_models; i++)
		if (aa->part[i].range[0] > maxrange)
			maxrange = aa->part[i].range[0];
	for (i = 0; i < ab->n_models; i++)
		if (ab->part[i].range[0] > maxrange)
			maxrange = ab->part[i].range[0];
	for (i = 0; i < bb->n_models; i++)
		if (bb->part[i].range[0] > maxrange)
			maxrange = bb->part[i].range[0];
	for (i = 0; i < 101 && !check_failed; i++) {
		dist = (i * maxrange)/100;
		dx = dy = dz = 0.0;
		if (i % 3 == 0) dx = dist;
		if (i % 3 == 1) dy = dist;
		if (i % 3 == 2) dz = dist;
		if (fabs(get_semivariance(ab, dx, dy, dz)) >
			sqrt(get_semivariance(aa, dx, dy, dz) * 
			get_semivariance(bb, dx, dy, dz))) {
			check_failed = 1; /* yes, the check failed */
			pr_warning("%s %d %s %d %s %d\n%s\n%s\n%s\n%s\n%s %g %g %g",
				"Cauchy-Schwartz violation: variogram",
				aa->id,",",bb->id, "and cross variogram", ab->id,
				"descriptors: ", aa->descr, bb->descr, ab->descr,
				"first failure on dx, dy and dz:", dx, dy, dz);
		}
	} /* for */
	if (check_failed)
		return 0;
	else
		return 1;
}

void check_variography(const VARIOGRAM **v, int n_vars)
/*
 * check for intrinsic correlation, linear model of coregionalisation
 * or else (with warning) Cauchy Swartz
 */
{
	int i, j, k, ic = 0, lmc, posdef = 1;
	MAT **a = NULL;
	double b;
	char *reason = NULL;

	if (n_vars <= 1)
		return;
/* 
 * find out if lmc (linear model of coregionalization) hold: 
 * all models must have equal base models (sequence and range)
 */
	for (i = 1, lmc = 1; lmc && i < get_n_vgms(); i++) {
		if (v[0]->n_models != v[i]->n_models) {
			reason = "number of models differ";
			lmc = 0;
		}
		for (k = 0; lmc && k < v[0]->n_models; k++) {
			if (v[0]->part[k].model != v[i]->part[k].model) {
				reason = "model types differ";
				lmc = 0;
			}
			if (v[0]->part[k].range[0] != v[i]->part[k].range[0]) {
				reason = "ranges differ";
				lmc = 0;
			}
		}
		for (k = 0; lmc && k < v[0]->n_models; k++)
			if (v[0]->part[k].tm_range != NULL) {
				if (v[i]->part[k].tm_range == NULL) {
					reason = "anisotropy for part of models";
					lmc = 0;
				} else if (
		v[0]->part[k].tm_range->ratio[0] != v[i]->part[k].tm_range->ratio[0] ||
		v[0]->part[k].tm_range->ratio[1] != v[i]->part[k].tm_range->ratio[1] ||
		v[0]->part[k].tm_range->angle[0] != v[i]->part[k].tm_range->angle[0] ||
		v[0]->part[k].tm_range->angle[1] != v[i]->part[k].tm_range->angle[1] ||
		v[0]->part[k].tm_range->angle[2] != v[i]->part[k].tm_range->angle[2]
				) {
					reason = "anisotropy parameters are not equal";
					lmc = 0;
				}
			} else if (v[i]->part[k].tm_range != NULL) {
				reason = "anisotropy for part of models";
				lmc = 0;
			}
	}
	if (lmc) {
/*
 * check for ic:
 */
		a = (MAT **) emalloc(v[0]->n_models * sizeof(MAT *));
		for (k = 0; k < v[0]->n_models; k++)
			a[k] = m_get(n_vars, n_vars);
		for (i = 0; i < n_vars; i++) {
			for (j = 0; j < n_vars; j++) { /* for all variogram triplets: */
				for (k = 0; k < v[0]->n_models; k++)
					a[k]->me[i][j] = v[LTI(i,j)]->part[k].sill;
			}
		}
		/* for ic: a's must be scaled versions of each other: */
		ic = 1;
		for (k = 1, ic = 1; ic && k < v[0]->n_models; k++) {
			b = a[0]->me[0][0]/a[k]->me[0][0];
			for (i = 0; ic && i < n_vars; i++)
				for (j = 0; ic && j < n_vars; j++)
					if (fabs(a[0]->me[i][j] / a[k]->me[i][j] - b) > EPSILON)
						ic = 0;	
		}
		/* check posdef matrices */
		for (i = 0, lmc = 1, posdef = 1; i < v[0]->n_models; i++) {
			posdef = is_posdef(a[i]);
			if (posdef == 0) {
				reason = "coefficient matrix not positive definite";
				if (DEBUG_COV) {
					printlog("non-positive definite coefficient matrix %d:\n", 
						i);
					m_logoutput(a[i]);
				}
				ic = lmc = 0;
			}
			if (! posdef)
				printlog(
				"non-positive definite coefficient matrix in structure %d", 
				i+1);
		}
		for (k = 0; k < v[0]->n_models; k++)
			m_free(a[k]);
		efree(a);

		if (ic) {
			printlog("Intrinsic Correlation found. Good.\n");
			return;
		} else if (lmc) {
			printlog("Linear Model of Coregionalization found. Good.\n");
			return;
		}
	}
/*
 * lmc does not hold: check on Cauchy Swartz
 */
	pr_warning("No Intrinsic Correlation or Linear Model of Coregionalization found\nReason: %s", reason ? reason : "unknown");
	if (gl_nocheck == 0) {
		pr_warning("[add `set nocheck = 1;' to the command file to ignore the following error]\n");
		ErrMsg(ER_IMPOSVAL, "variograms do not satisfy a legal model");
	}
	printlog("Now checking for Cauchy-Schwartz inequalities:\n");
	for (i = 0; i < n_vars; i++)
		for (j = 0; j < i; j++)
			if (is_valid_cs(v[LTI(i,i)], v[LTI(j,j)], v[LTI(i,j)])) {
				printlog("variogram(%s,%s) passed Cauchy-Schwartz\n",
					name_identifier(j), name_identifier(i));
			} else
				pr_warning("Cauchy-Schwartz inequality found for variogram(%s,%s)",
						name_identifier(j), name_identifier(i) );
	return;
}

/* from meschach matrix library (c) , see matrix[2].h
 try CHfactor -- Cholesky L.L' factorisation of A in-situ */
static int is_posdef(MAT *A) {
	u_int	i, j, k;
	Real	sum, tmp;

	for (k = 0; k < A->n; k++)
	{	
		/* do diagonal element */
		sum = A->me[k][k];
		for (j = 0; j < k; j++)
		{
			tmp = A->me[k][j];
			sum -= tmp*tmp;
		}
		/*
		if (sum <= 0.0)
			return 0;
		A->me[k][k] = sqrt(sum);
		*/
		if (sum < -gl_zero)
			return 0;
		if (sum < gl_zero)
			A->me[k][k] = sqrt(gl_zero);
		else
			A->me[k][k] = sqrt(sum);

		/* set values of column k */
		for (i = k + 1; i < A->n; i++)
		{
			sum = A->me[i][k];
			sum -= __ip__(A->me[i],A->me[k],(int)k);
			A->me[j][i] = A->me[i][j] = sum/A->me[k][k];
		}
	}
	return 1;
}

double transform_norm(const ANIS_TM *tm, double dx, double dy, double dz) {
/* returns variogram distance given dx, dy, dz and VARIOGRAM v */

	double dist = 0.0, tmp;
	int i;

	if (dx == 0.0 && dy == 0.0 && dz == 0.0)
		return 0.0;
	if (tm != NULL) {
		for (i = 0, tmp = 0.0; i < 3; i++) {
			tmp = tm->tm[i][0] * dx + tm->tm[i][1] * dy + tm->tm[i][2] * dz;
			dist += tmp * tmp;
		}
		return sqrt(dist);
	} 
	return sqrt((dx * dx) + (dy * dy) + (dz * dz));
}

double da_general(VGM_MODEL *part, double h) {
	int i;
	double low, high, range, r[NRANGEPARS];

	for (i = 0; i < NRANGEPARS; i++) {
		if (is_mv_double(&(part->range[i])))
			set_mv_double(&(r[i]));
		else
			r[i] = part->range[i];
	}
	range = MAX(1e-20, part->range[0]);
	r[0] = range * (1.0 + DA_DELTA);
	low = part->fnct(h, r);
	r[0] = range * (1.0 - DA_DELTA);
	high = part->fnct(h, r);
	return part->sill * (low - high) / (2.0 * range * DA_DELTA);
}

int push_variogram_model(VARIOGRAM *v, VGM_MODEL part) {
	int i, max_id, where = -1;
/*
 * add the part submodel to v (if part.id < 0) or else
 * modify the appropriate part of v, having the id of part.id.
 * do a lot of checks, and set .fn and .da_fn functions.
 */

	if (v->n_models == v->max_n_models) {
		v->part = (VGM_MODEL *) erealloc(v->part, 
				(v->max_n_models + INIT_N_VGMM) * sizeof(VGM_MODEL));
		for (i = v->max_n_models; i < v->max_n_models + INIT_N_VGMM; i++)
			init_variogram_part(&(v->part[i]));
		v->max_n_models += INIT_N_VGMM;
#ifndef USING_R
		printf("enlarging v->max_n_models\n");
#endif
	}
	/*
	 * check some things: 
	 */
	if (part.model == NOT_SP)
		ErrMsg(ER_IMPOSVAL, "model NSP not allowed in variogram structure");
	if (part.range[0] < 0.0)
		ErrMsg(ER_RANGE, "variogram range cannot be negative");
	if (part.model == LINEAR) {
		if (part.range[0] == 0.0)
			part.fit_range = 0;
	} else if (part.model == NUGGET || part.model == INTERCEPT || 
			part.model == MERROR) {
		part.fit_range = 0;
		if (part.range[0] > 0.0) 
			ErrMsg(ER_RANGE, "range must be zero");
	} else if (part.range[0] == 0.0) 
		ErrMsg(ER_RANGE, "range must be positive");
	if (part.model == POWER && part.range[0] > 2.0)
		ErrMsg(ER_RANGE, "power model range (parameter) cannot exceed 2.0");
	if (part.model == EXCLASS && part.range[1] > 2.0)
		ErrMsg(ER_RANGE, "exponentical class model shape parameter cannot exceed 2.0");

	if (part.id < 0) {
		where = v->n_models;
		v->n_models++;
		for (i = max_id = 0; i < v->n_models; i++)
			max_id = MAX(v->part[i].id, max_id);
		part.id = max_id + 1;
	} else { /* search in list: */
		for (i = 0; where < 0 && i < v->n_models; i++)
			if (v->part[i].id == part.id)
				where = i;
		assert(where >= 0); /* i.e., it should really be in the list */
	}

	if (v->isotropic)
		v->isotropic = (part.tm_range == NULL);

	/* 
	 * check that the .fn and .da_fn functions in v_models 
	 * will indeed be the correct ones: 
	 */
	assert(part.model == v_models[part.model].model);

	v->part[where] = part;
	v->part[where].fnct = v_models[part.model].fn;
	v->part[where].da_fnct = v_models[part.model].da_fn;

	return part.id;
}

VGM_MODEL_TYPE which_variogram_model(const char *m) {
	char s[4];
	int i;

	strncpy(s, m, 3);
	s[0] = toupper(s[0]);
	s[1] = tolower(s[1]);
	s[2] = tolower(s[2]);
	s[3] = '\0';
	for (i = 1; v_models[i].name != NULL; i++)
		if (almost_equals(s, v_models[i].name))
			return v_models[i].model;
	return NOT_SP;
}

double relative_nugget(VARIOGRAM *v) {
	int i;
	double nug = 0.0, sill = 0.0;
	
	assert(v->n_models != 0);

	if (v->n_models == 1)
		return (v->part[0].model == NUGGET ? 1.0 : 0.0);

	for (i = 0; i < v->n_models; i++) {
		if (v->part[i].model == NUGGET)
			nug += v->part[i].sill;
		else
			sill += v->part[i].sill;
	}
	assert(nug + sill > 0.0);
	return (nug/(nug+sill));
}

#ifndef USING_R
int vario(int argc, char **argv) {
/* model from to nsteps */
	double dist, from, to;
	int i, is_vgm, nsteps = 0;
	VARIOGRAM vgm;

	is_vgm = almost_equals(argv[0], "se$mivariance");
	if (argc < 3) {
		printlog("usage: %s variogram_model dist [to_dist [n_intervals]]\n", argv[0]);
		exit(0);
	}
	init_variogram(&vgm);
	vgm.id = 0;
	if (read_variogram(&vgm, string_dup(argv[1])))
		ErrMsg(ER_SYNTAX, argv[1]);
	if (read_double(argv[2], &from))
		ErrMsg(ER_RDFLT, argv[2]);
	if (argc >= 4) {
		if (read_double(argv[3], &to))
			ErrMsg(ER_RDFLT, argv[3]);
		nsteps = 1;
	} else
		to = from;
	if (argc >= 5)
		if (read_int(argv[4], &nsteps))
			ErrMsg(ER_RDINT, argv[4]);
	if (DEBUG_DUMP)
		logprint_variogram(&vgm, 1);
	if (nsteps < 0)
		ErrMsg(ER_RANGE, "n_steps must be >= 0");
	dist = from;
	for (i = 0; i <= nsteps; i++) {
		printlog("%g %g\n", dist, (is_vgm ? 
			get_semivariance(&vgm, dist, 0, 0) : 
			get_covariance(&vgm, dist, 0, 0)));
		if (i < nsteps) /* nsteps > 0 */
			dist += (to - from)/(1.0*nsteps);
	}	
	return 0;
}
#endif

FIT_TYPE fit_int2enum(int fit) {
	switch (fit) {
		case 0: return NO_FIT;
		case 1: return WLS_FIT; 
		case 2: return WLS_FIT_MOD;
		case 3: return WLS_GNUFIT;
		case 4: return WLS_GNUFIT_MOD; 
		case 5: return MIVQUE_FIT;
		case 6: return OLS_FIT;
		case 7: return WLS_NHH;
	}
	ErrMsg(ER_IMPOSVAL, "invalid value for fit");
	return OLS_FIT; /* never reached */
}

FIT_TYPE fit_shift(FIT_TYPE now, int next) {
	switch (now) {
		case NO_FIT: return (next ? WLS_FIT : WLS_NHH);
		case WLS_FIT: return (next ? WLS_FIT_MOD : NO_FIT); 
		case WLS_FIT_MOD: return (next ? WLS_GNUFIT : WLS_FIT);
		case WLS_GNUFIT: return (next ? WLS_GNUFIT_MOD : WLS_FIT_MOD);
		case WLS_GNUFIT_MOD: return (next ? MIVQUE_FIT : WLS_GNUFIT); 
		case MIVQUE_FIT: return (next ? OLS_FIT : WLS_GNUFIT_MOD);
		case OLS_FIT: return (next ? WLS_NHH : MIVQUE_FIT);
		case WLS_NHH: return (next ? NO_FIT : OLS_FIT);
	}
	return NO_FIT; /* never reached */
}

DO_AT_ZERO zero_int2enum(int zero) {
	switch(zero) {
		case 0: return ZERO_DEFAULT;
		case 1: return ZERO_INCLUDE;
		case 2: return ZERO_AVOID;
		case 3: return ZERO_SPECIAL;
	}
	ErrMsg(ER_IMPOSVAL, "invalid value for zero");
	return ZERO_DEFAULT; /* never reached */
}

DO_AT_ZERO zero_shift(DO_AT_ZERO now, int next) {
	if (next) {
		switch(now) {
			case ZERO_DEFAULT: return ZERO_INCLUDE;
			case ZERO_INCLUDE: return ZERO_AVOID;
			case ZERO_AVOID: return ZERO_SPECIAL;
			case ZERO_SPECIAL: return ZERO_DEFAULT;
		}
	} else {
		switch(now) {
			case ZERO_DEFAULT: return ZERO_SPECIAL;
			case ZERO_INCLUDE: return ZERO_DEFAULT;
			case ZERO_AVOID: return ZERO_INCLUDE;
			case ZERO_SPECIAL: return ZERO_SPECIAL;
		}
	}
	return ZERO_DEFAULT; /* never reached */
}

VGM_MODEL_TYPE model_shift(VGM_MODEL_TYPE now, int next) {
	int i;

	if (now == NOT_SP)
		return NUGGET;

	for (i = 1; v_models[i].model != NOT_SP; i++) {
		if (v_models[i].model == now) {
			if (next) {
				if (v_models[i+1].model == NOT_SP)
					return now;
				else
					return v_models[i+1].model;
			} else {
				if (v_models[i-1].model == NOT_SP)
					return now;
				else
					return v_models[i-1].model;
			}
		}
	}
	return NUGGET; /* never reached */
}

double effective_range(const VARIOGRAM *v) {
	int i;
	double er = 0.0, range;
	for (i = 0; i < v->n_models; i++) {
		switch (v->part[i].model) {
			case EXPONENTIAL: range = 3 * v->part[i].range[0]; break;
			case GAUSSIAN: range = sqrt(3) * v->part[i].range[0]; break;
			default: range = v->part[i].range[0]; break;
		}
		er = MAX(er, range);
	}
	return er;
}

int get_n_variogram_models(void) {
	int i, n = 0;
	for (i = 1; v_models[i].model != NOT_SP; i++)
		n++;
	return(n);
}

void push_to_v(VARIOGRAM *v, const char *mod, double sill, double *range, 
		int nrangepars, double *d, int fit_sill, int fit_range) {
	VGM_MODEL vm;
	int i;

	init_variogram_part(&vm);
	vm.model = which_variogram_model(mod);
	if (nrangepars > NRANGEPARS)
		ErrMsg(ER_IMPOSVAL, "too many range parameters");
	for (i = 0; i < nrangepars; i++)
		vm.range[i] = range[i];
	vm.sill = sill;
	vm.fit_sill = fit_sill;
	vm.fit_range = fit_range;
	if (d != NULL && d[0] != -9999.0)
		vm.tm_range = get_tm(d);
#ifdef USING_R
	if (vm.model == STEIN && range[1] > 100.0) {
		vm.model = GAUSSIAN;
		vm.range[1] = 0.0;
		pr_warning("kappa values over 100 overflow gammafn: taking Gaussian approximation");
	}
#endif
	push_variogram_model(v, vm);
}

void push_to_v_table(VARIOGRAM *v, double maxdist, int length, double *values,
			double *anis) {
	int i;

	v->table = (COV_TABLE *) emalloc(sizeof(COV_TABLE));
	v->table->n = length;
	v->table->maxdist = maxdist;
	v->table->values = (double *) emalloc(length * sizeof(double));
	for (i = 0; i < length; i++)
		v->table->values[i] = values[i];
	if (anis != NULL)
		v->table->tm_range = get_tm(anis);
	else
		v->table->tm_range = NULL;
}

static ANIS_TM *get_tm(double anis[5]) {
/* Part of this routine was taken from GSLIB, first edition:
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C                                                                      %
C Copyright (C) 1992 Stanford Center for Reservoir Forecasting.  All   %
C rights reserved.  Distributed with: C.V. Deutsch and A.G. Journel.   %
C ``GSLIB: Geostatistical Software Library and User's Guide,'' Oxford  %
C University Press, New York, 1992.                                    %
C                                                                      %
C The programs in GSLIB are distributed in the hope that they will be  %
C useful, but WITHOUT ANY WARRANTY.  No author or distributor accepts  %
C responsibility to anyone for the consequences of using them or for   %
C whether they serve any particular purpose or work at all, unless he  %
C says so in writing.  Everyone is granted permission to copy, modify  %
C and redistribute the programs in GSLIB, but only under the condition %
C that this notice and the above copyright notice remain intact.       %
C                                                                      %
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/
	int i;
	double alpha, beta, theta, sina, sinb, sint, cosa, cosb, cost, afac1, afac2;
	ANIS_TM *t = NULL;

/* 
	About naming convention:

	gstat     GSLIB
	===============
	anis[0]    ang1 (first anis. par. for 2D)
	anis[1]    ang2
	anis[2]    ang3
	anis[3]   anis1 (second anis. par. for 2D)
	anis[4]   anis2
*/

#define ANIS_ERR(x) message("parsing anis. pars. %g,%g,%g,%g,%g -- error on %g\n", \
	anis[0],anis[1],anis[2],anis[3],anis[4],x)
#define DEG2RAD (PI/180.0)

	for (i = 0; i < 3; i++) {
		if (anis[i] < 0 || anis[i] >= 360) {
			ANIS_ERR(anis[i]);
			ErrMsg(ER_RANGE, "this value should be in [0..360>");
		}
	}
	for (i = 3; i < 5; i++) {
		if (anis[i] <= 0.0 || anis[i] > 1.0) {
			ANIS_ERR(anis[i]);
			ErrMsg(ER_RANGE, "this value should be in <0..1]");
		}
	}

	/* from GSLIB: */
	if (anis[0] >= 0.0 && anis[0] < 270)
		alpha = (double) (90.0 - anis[0]) * DEG2RAD;
	else
		alpha = (double) (450.0 - anis[0]) * DEG2RAD;
	beta = -1.0 * (double) anis[1] * DEG2RAD;
	theta =       (double) anis[2] * DEG2RAD;

	sina = sin(alpha);
	sinb = sin(beta);
	sint = sin(theta);
	cosa = cos(alpha);
	cosb = cos(beta);
	cost = cos(theta);

	afac1 = 1.0 / MAX((double) anis[3], (double) EPSILON);
	afac2 = 1.0 / MAX((double) anis[4], (double) EPSILON);

	t = emalloc(sizeof(ANIS_TM));

	t->angle[0] = anis[0];
	t->angle[1] = anis[1];
	t->angle[2] = anis[2];
	t->ratio[0] = anis[3];
	t->ratio[1] = anis[4];
	t->tm[0][0] =       (cosb * cosa);
	t->tm[0][1] =       (cosb * sina);
	t->tm[0][2] =       (-sinb);
	t->tm[1][0] = afac1*(-cost*sina + sint*sinb*cosa);
	t->tm[1][1] = afac1*(cost*cosa + sint*sinb*sina);
	t->tm[1][2] = afac1*( sint * cosb);
	t->tm[2][0] = afac2*(sint*sina + cost*sinb*cosa);
	t->tm[2][1] = afac2*(-sint*cosa + cost*sinb*sina);
	t->tm[2][2] = afac2*(cost * cosb);
	return t;
}
