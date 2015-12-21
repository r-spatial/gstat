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
 * glvars.c: has global variables, choose defaults, check setting integrity
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "defs.h"
#include "userio.h"
#include "debug.h"
#include "data.h"
#include "utils.h"
#include "vario.h"
#include "lex.h" /* read_variogram() */
#include "random.h"
#include "predict.h"
#include "defaults.h" /* default values for gl_* variables */

#ifndef USING_R
# include "writecmd.h" /* printlog_glvars() */
#endif

#include "glvars.h"
#include "gls.h"
#include "block.h"
#include "polygon.h"
/* 
 * this module contains all global variables, along with some 
 * routines for initialisation (of structures), resizing during runtime,
 * checks after reading command file and before starting real calculations,
 * obtaining some "local" global structures, printing all that's held,
 * and so on. get_n_vars() etc. are for clarity that these variables are not
 * changed outside the area where they're set.
 */

static void init_gstat_data(int n);
static void clean_up(void);

/*
 * global variables (glvars.h, defautls.h):
 */
int debug_level; /* debug level */
int gl_coincide; /* do the variable locations coincide? */
int gl_cressie; /* use cressie's estimator ? */
int gl_dots; /* treshold nr for points -> dots, gnuplot */
int gl_fit; /* do not fit a variogram */
int gl_gauss; /* gaussian quadr. block covariances ? */
int gl_iter; /* max. n. iter for mivque estimates */
int gl_jgraph; /* do jgraph plot in batch mode ? */
int gl_lhs; /* apply Latin hypercube sampling to Gaussian simulations */
int gl_longlat; /* indicates whether coordinates are longitude/latitude, or Euclidian */
int gl_nblockdiscr; /* block discrimination in each dimension */
int gl_n_intervals; /* n variogram intervals */
int gl_n_marginals; /* the n marginal distributions */
int gl_nsim; /* number of simultanious simulations */
int gl_n_uk; /* min. # of cs points to use ok */
int gl_numbers; /* plot numbers on variogram plot? */
int gl_nocheck; /* do not check LMC/IC ? */
int gl_order; /* do order relation correction */
int gl_plotweights; /* plot kriging weights? */
int gl_register_pairs; /* register sample variogram pairs? */
int gl_rowwise; /* deal with raster maps row-wise, or as complete blocks? */
int gl_rp; /* follow random path for gs/is? */
int gl_secure; /* disallow system() and popen()? */
int gl_seed; /* seed is set? */
int gl_sim_beta; /* simulation mode for beta: 0 multiv GLS, 1 univ GLS, 2 OLS */
int gl_spiral; /* do spiral search if possible? */
int gl_split; /* see nsearch.c: was Q_SPLIT_AT */
int gl_sym_ev; /* default symmetric ps.cr.v./cr.cv. ? */
int gl_gls_residuals; /* calc. gls residuals? */
int gl_sparse; /* use sparse covariance matrices? */
int gl_xvalid; /* do cross validation on first variable */
int gl_zero_est; /* est. variogram at h 0 seperately? */
double *gl_bounds; /* boundaries semivariogram intervals */
double *gl_marginal_values;
double gl_cutoff; /* variogram cutoff */
double gl_fit_limit; /* convergence criterion on fit */
double gl_fraction; /* fraction of max dist for cutoff */
double gl_gcv; /* generalized covariance constant */
double gl_idp; /* default inverse distance power */
double gl_iwidth; /* variogram class width */
double gl_quantile; /* sample quantile */
double gl_zmap; /* height of the map */
double gl_cn_max; /* check maximum condition number. Default: don't */
double gl_alpha; /* alpha, beta, tol_[hor|ver]: anisotropy parameters */
double gl_beta;
double gl_tol_hor;
double gl_tol_ver;
double gl_zero; /* zero tolerance; 2-squared */

char **gl_marginal_names;
char *gl_display;
char *gl_mv_string;
char *gl_plotfile;  /* gnuplot file name */
char *gl_gnuplot;
char *gl_pager;
char *gl_format;
char *gl_gnuplot35;
char *gl_gpterm;
char *command_line;
char *command_file_name;
char *logfile_name;
char *o_filename;
char *argv0;
FILE *plotfile;

const METHODS methods[] = { /* methods and codes */
	{ NSP,      0, "nsp" }, /* do nothing */
	{ UIF,      0, "ui" },    /* variogram modelling user interface */
	{ OKR,      0, "ok" },   /* ordinary kriging */
	{ UKR,      0, "uk" },   /* universal kriging */
	{ SKR,      0, "sk" },  /* simple kriging */
	{ IDW,      0,  "id" }, /* inverse distance interpolation */
	{ MED,      0, "med" }, /* (local) sample median or quantiles */
	{ NRS,      0, "n$r" }, /* neighbourhood size */
	{ LSLM,     0, "tr$end" },  /* uncorrelated (or weighted) linear model */
	{ GSI,      1,  "gs" }, /* gaussian (conditional) simulation */
	{ ISI,      1,  "is" }, /* indicator (conditional) simulation */
	{ MAPVALUE, 0, "map" },  /* mask map value at data location */
	{ SEM,      0, "se$mivariogram" }, /* sample (cross) semivariance */
	{ COV,      0, "co$variogram" },  /* sample (cross) covariance */
	{ SPREAD,   0, "di$stance" }, /* distance to nearest sample */
	{ XYP,      0, "xy" },  /* x and y coordinate of location */
	{ POLY,     0, "point-in-polygon" }, /* point-in-polygon */
	{ DIV,      0, "div" }, /* diversity and modus */
	{ SKEW,     0, "skew" }, /* skewness and kurtosis */
	{ LSEM,     0, "lsem" }, /* locally estimated/fitted variogram parameters */
	{ TEST,     0, "test" },  /* do-nothing? */
	{ NSP,      0, NULL } /* terminating field */
};

/*
 * "local" database, accesible through get_* functions
 */
static VARIOGRAM **vgm = NULL;
static DATA **data = NULL;
static char **outfile_names = NULL, **ids = NULL, **mask_names = NULL;
static DATA *valdata = NULL;
static DATA *data_area = NULL; /* area that discretises block */
static DPOINT block;
static METHOD method = NSP;
static int n_vars = 0, n_masks = 0;
static int n_last = 0, n_v_last = 0, n_o_last = 0;
static MODE mode = MODE_NSP; /* MODE_NSP, SIMPLE, STRATIFY or MULTIVARIABLE */
static int n_edges = 0; /*really n of edges files */
static char **edges_names = NULL;  /* names of edges files, n_edges long */
static int *n_edges_polys = NULL;  /* number of edges in every edge file */
static POLYGON **edges = NULL;            /* per-file edges */

int init_global_variables(void) {
/*
 * global variables init. (glvars.h; defautls.h):
 */
 	method             = NSP;
	mode               = MODE_NSP;
	debug_level        = DB_NORMAL;
	gl_coincide        = DEF_coincide;
	gl_cressie         = DEF_cressie;
	gl_dots            = DEF_dots;
	gl_fit             = DEF_fit;
	gl_gauss           = DEF_gauss;
	gl_iter            = DEF_iter;
	gl_jgraph          = DEF_jgraph;
	gl_lhs             = DEF_lhs;
	gl_longlat         = DEF_longlat;
	gl_nblockdiscr     = DEF_nblockdiscr;
	gl_n_intervals     = DEF_intervals;
	gl_n_marginals     = DEF_n_marginals;
	gl_nsim            = DEF_nsim;
	gl_n_uk            = DEF_n_uk;
	gl_numbers         = DEF_numbers;
	gl_nocheck         = DEF_nocheck;
	gl_order           = DEF_order;
	gl_plotweights     = DEF_plotweights;
	gl_register_pairs  = DEF_pairs;
	gl_rowwise         = DEF_rowwise;
	gl_rp              = DEF_rp;
	gl_secure          = DEF_secure;
	gl_seed            = DEF_seed;
	gl_sim_beta        = DEF_sim_beta;
	gl_spiral          = DEF_spiral;
	gl_split           = DEF_split;
	gl_sym_ev          = DEF_sym_ev;
	gl_gls_residuals   = DEF_gls_residuals;
	gl_sparse          = DEF_sparse;
	gl_xvalid          = DEF_xvalid;
	gl_zero_est        = DEF_zero_est;
	gl_bounds          = DEF_bounds;
	gl_marginal_values = DEF_marginal_values;
	gl_cutoff       = DEF_cutoff;
	gl_fit_limit    = DEF_fit_limit;
	gl_fraction     = DEF_fraction;
	gl_gcv          = DEF_gcv;
	gl_idp          = DEF_idp;
	gl_iwidth       = DEF_iwidth;
	gl_quantile     = DEF_quantile;
	gl_zmap         = DEF_zmap;
	gl_cn_max       = DEF_cn_max;
	gl_alpha        = DEF_alpha;
	gl_beta         = DEF_beta;
	gl_tol_hor      = DEF_tol_hor;
	gl_tol_ver      = DEF_tol_ver;
	gl_zero         = DEF_zero;
	
	gl_marginal_names = DEF_marginal_names;
	gl_display       = DEF_display;
	gl_mv_string     = DEF_mv_string;
	gl_plotfile       = DEF_plotfile; 
	gl_gnuplot       = DEF_gnuplot;
	gl_pager         = DEF_pager;
	gl_format        = DEF_format;
	gl_gnuplot35     = NULL;
	gl_gpterm        = DEF_gpterm;
	command_line     = NULL;
	command_file_name = NULL;
	logfile_name     = NULL;
	o_filename       = NULL;
	argv0            = NULL;
	plotfile         = NULL;
	
	init_gstat_data(0);
	/* EJPXX 
	 * 	if (valdata == NULL)
	 * 	*/
	valdata = init_one_data(valdata);
	block.x = block.y = block.z = 0.0;
	set_mv_double(&gl_zmap);
	get_covariance(NULL, 0, 0, 0);
	return 0;
}

void push_mask_name(const char *name) {
	if (mask_names == NULL)
		mask_names = (char **) emalloc((n_masks + 1) * sizeof(char *));
	else
		mask_names = (char **) erealloc(mask_names,
				(n_masks + 1) * sizeof(char *));
	mask_names[n_masks] = string_dup(name);
	n_masks++;
}

void push_bound(double value) {
	static int n_bound;

	if (gl_bounds == NULL) {
		n_bound = 0;
		gl_bounds = (double *) emalloc((n_bound+2) * sizeof(double));
	} else
		gl_bounds = (double *) erealloc(gl_bounds,(n_bound+2) * sizeof(double));
	gl_bounds[n_bound] = value;
	gl_bounds[n_bound + 1] = -1.0;
	if (n_bound > 0 && gl_bounds[n_bound] <= gl_bounds[n_bound-1])
		ErrMsg(ER_IMPOSVAL, "bounds must be strictly increasing");
	n_bound++;
}

const char *get_mask_name(int i) {
	if (i >= n_masks || i < 0)
		return NULL;
	return mask_names[i];
}

int get_n_masks(void) {
	return n_masks;
}

int n_variograms_set(void) {
	int i, n;
	for (i = 0, n = 0; i < get_n_vgms(); i++)
		if (vgm[i] != NULL && vgm[i]->id >= 0)
			n++;
	return n;
}

static void init_gstat_data(int n) {
	int i, n_vgms, n_outfl;

	n_vgms = (n * (n + 1))/2;
	n_outfl = n + n_vgms;
	if (n <= n_last)
		return;
	data = (DATA **) erealloc(data, n * sizeof(DATA *));
	for (i = n_last; i < n; i++)
		data[i] = init_one_data(NULL);
	vgm = (VARIOGRAM **) erealloc(vgm, n_vgms * sizeof(VARIOGRAM *));
	for (i = n_v_last; i < n_vgms; i++)
		vgm[i] = NULL;
	outfile_names = (char **) erealloc (outfile_names, n_outfl * sizeof(char *));
	for (i = n_o_last; i < n_outfl; i++)
		outfile_names[i] = NULL;
	n_last = n;
	n_o_last = n_outfl;
	n_v_last = n_vgms;
	n_vars = n;
	return;
}

int which_identifier(const char *id) {
	int i; 

	assert(id);

	for (i = 0; i < n_vars; i++) {
		if (ids[i] == NULL) {
			ids[i] = string_dup(id);
			return i;
		} else /* ids[i] == NULL: */ {
			if (strcmp(ids[i], id) == 0)
				return i;
		}
	}
	/* else: extend data space */
	if (strlen(id) > MAX_ID_LENGTH) {
		pr_warning("maximum length for an identifier is %d", MAX_ID_LENGTH);
		ErrMsg(ER_IMPOSVAL, "which_identifier()");
	}
	n_vars++;
	ids = (char **) erealloc(ids, n_vars * sizeof(char *));
	ids[n_vars - 1] = string_dup(id);
	init_gstat_data(n_vars);
	return n_vars - 1;
}

const char *name_identifier(int i) {
	static const char *cp_val = "data()", *cp_area = "area";
	switch (i) {
		case ID_OF_VALDATA:
			return cp_val;
		case ID_OF_AREA:
			return cp_area;
		default:
			if (i >= get_n_vars() || i < 0) {
				pr_warning("i = %d", i);
				ErrMsg(ER_RANGE, "name_identifier(i): i outside range");
			}
			return ids[i];
	}
}

const char *what_is_outfile(int i) {
	static char **what = NULL;
	static int sizeof_what = 0;
	int j, k;

	if (i < 0) {
		if (what != NULL) {
			for (j = 0; j < sizeof_what; j++)
				efree(what[j]);
			efree(what);
			what = NULL;
			sizeof_what = 0;
		}
		return NULL;
	}
	if (i > get_n_outfile())
		ErrMsg(ER_RANGE, "what_is_outfile(i): i outside range");
	if (what == NULL) {
		sizeof_what = get_n_outfile();
		what = (char **) emalloc(get_n_outfile() * sizeof(char *));
		for (j = 0; j < get_n_outfile(); j++) {
			what[j] = (char *) emalloc (100 * sizeof(char)); 
			what[j][0] = '\0';
		}
		if (get_mode() == STRATIFY) {
			sprintf(what[0], "[predicted value, per stratum]");
			sprintf(what[1], "[prediction variance, per stratum]");
		} else if (get_method() == MAPVALUE) {
			for (j = 0; j < get_n_outfile(); j++) {
				if (j < get_n_masks()) {
					strncpy(what[j], get_mask_name(j), 39);
					what[j][39] = '\0'; /* Konstantin Malakhanov, 30/3/99 */
				} else
					sprintf(what[j], "%s", "[ignore]");
			}
		} else {
			for (j = 0; j < 2 * n_vars; j++) {
				if (j % 2 == 0)
					sprintf(what[j], "pred(%s)", name_identifier(j / 2));
						else
					sprintf(what[j], "var(%s)", name_identifier(j / 2));
			}
			for (j = 0; j < get_n_vars(); j++) {
				for (k = 0; k < j; k++) {
					what[2 * get_n_vars() + LTI2(j,k)] =
						(char *) emalloc((40 + strlen(name_identifier(j)) +
							strlen(name_identifier(k))) *sizeof(char));
					sprintf(what[2 * get_n_vars() + LTI2(j,k)],
						"cov(%s,%s)", name_identifier(j),
						name_identifier(k));
				}
			}
		}
	}
	return what[i];
}

const char *method_string(METHOD i) {
	static char mstr[100];
	char *str, *co, *un, *gsum = "";

	if ((i == ISI || i == GSI) && gl_n_uk == DEF_n_uk &&
			get_n_beta_set() != get_n_vars())
		gsum = " with unknown means";

	str = (get_mode() == STRATIFY ? "stratified " : "");
	un = (get_n_vars() > 0 && data[0]->dummy ? "un" : "");
	co = (get_mode() == MULTIVARIABLE ? "co" : "");

	switch (i) {
		case NSP:
			sprintf(mstr, "exit");
			break;
		case TEST:
			sprintf(mstr, "Test Option");
			break;
		case UIF:
			sprintf(mstr, "starting interactive mode");
			break;
		case MAPVALUE:
			sprintf(mstr, "mask map values on data() locations");
			break;
		case SEM:
			sprintf(mstr, "calculating sample variogram");
			break;
		case COV:
			sprintf(mstr, "calculating sample covariogram");
			break;
		case SPREAD:
			sprintf(mstr, "spread value (distance to nearest observation) on output");
			break;
		case XYP:
			sprintf(mstr, "mask map coordinates (x,y) on first and second output vars");
			break;
		case IDW:
			sprintf(mstr, "%sinverse distance weighted interpolation", str);
			break;
		case MED:
			if (gl_quantile == 0.5)
				sprintf(mstr, "%smedian estimation", str);
			else
				sprintf(mstr, "%s%g-quantile estimation", str, gl_quantile);
			break;
		case NRS:
			sprintf(mstr, "(%s:) neighbourhood size on first output variable", str);
			break;
		case LSLM:
			if (n_variograms_set())
				sprintf(mstr, "%sgeneralized least squares trend estimation", str);
			else
				sprintf(mstr, "%sordinary or weighted least squares prediction", str);
			break;
		case OKR:
			sprintf(mstr, "using %sordinary %skriging", str, co);
			break;
		case SKR:
			sprintf(mstr, "using %ssimple %skriging", str, co);
			break;
		case UKR:
			sprintf(mstr, "using %suniversal %skriging", str, co);
			break;
		case GSI:
			sprintf(mstr, "using %s%sconditional Gaussian %ssimulation%s",
				str, un, co, gsum);
			break;
		case ISI:
			sprintf(mstr, "using %s%sconditional indicator %ssimulation",
				str, un, co);
			break;
		case POLY:
			sprintf(mstr, "point-in-polygon");
			break;
		case DIV:
			sprintf(mstr, "within-neighbourhood diversity and modus");
			break;
		case SKEW:
			sprintf(mstr, "skewness and kurtosis");
			break;
		case LSEM:
			sprintf(mstr, "local semivariance or locally fitted semivariogram parameters");
			break;
	}
	return mstr;
}

int get_n_vars(void) {
	return n_vars;
}

int get_n_beta_set(void) {
	int i, nbeta;

	for (i = nbeta = 0; i < get_n_vars(); i++)
		if (data[i]->beta != NULL)
			nbeta++;
	return nbeta;
}

MODE get_mode(void) {
	return mode;
}

int get_n_vgms(void) {
	int n;
	n = get_n_vars();
	return (n * (n + 1))/2;
}

int get_n_outfile(void) {
	return get_n_vars() + get_n_vgms();
}

const char **get_outfile_name(void) {
	return (const char **) outfile_names;
}

VARIOGRAM *get_vgm(int i) {
	assert(i < get_n_vgms());
	assert(i >= 0);
	if (vgm[i] == NULL)
		vgm[i] = init_variogram(NULL);
	return vgm[i];
}

DATA **get_gstat_data(void) {
	return data;
}

DATA *get_dataval(void) {
	return valdata;
}

DATA *get_data_area(void) {
	return data_area;
}

DATA *create_data_area(void) {
	data_area = init_one_data(NULL);
	return data_area;
}

DPOINT *get_block_p(void) {
	return &block;
}

const char *get_outfile_namei(int i) {
	if (i >= get_n_outfile() || i < 0)
		ErrMsg(ER_RANGE, "get_outfile_namei(i) i outside range");
	return outfile_names[i];
}

METHOD get_method(void) {
	return method;
}

void set_method(METHOD m) {
	method = m;
	return;
}

int is_simulation(METHOD m) {
	assert(methods[m].m == m);
	return methods[m].is_simulation;
}

double max_block_dimension(int reset) {
	static double dim = -1.0;

	if (reset)
		dim = -1.0;
	if (dim < 0.0) {
		if (data_area != NULL)
			dim = data_block_diagonal(data_area);
		else
			dim = sqrt(SQR(block.x)+SQR(block.y)+SQR(block.z));
	} 
	return dim;
}

#ifndef USING_R
int dump_all(void) {
	int i;
	char *cp = NULL;

	printlog("# global variables:\n");
	printlog("%s", sprint_glvars(1));
	if (n_masks > 0) {
		printlog("masks: ");
		for (i = 0; i < n_masks; i++)
			printlog("'%s'%s", get_mask_name(i), i < n_masks-1 ? ", " : ";\n");
	}
	if (gl_marginal_names) {
		printlog("marginals: ");
		for (i = 0; i < gl_n_marginals; i++)
			printlog("'%s'%s", gl_marginal_names[i],
				i < gl_n_marginals - 1 ? ", " : ";\n");
	}
	for (i = 0; i < get_n_outfile(); i++)
		if (outfile_names[i] != NULL)
			printlog("outfile[%d]: '%s'; # %s\n", i, outfile_names[i],
				what_is_outfile(i));
	for (i = 0; i < get_n_vgms(); i++)
		if (vgm[i] != NULL && vgm[i]->descr)
			logprint_variogram(vgm[i], 1);
	printlog("blocksize: dx = %g, dy = %g, dz = %g;\n", 
		block.x, block.y, block.z);
	for (i = 0; i < get_n_vars(); i++) {
		print_data_line(data[i], &cp);
		printlog("data(%s): %s\n", name_identifier(i), cp);
	}
	if (valdata->id > -1) {
		print_data_line(valdata, &cp);
		printlog("data(): %s\n", cp);
	}
	if (get_n_edges())
		report_edges();
	printlog("method: %s\n", method_string(get_method()));
	return 0;
}
#endif

void setup_valdata_X(DATA *d) {
/* 
 * fills '0'-X columns (intercept) at the right place
 * e.g. data(a): .. X=1&2;data(b): .. X=3&4; data(): .. X=1&2&3&4;
 * this leads to "0 1 2 3 4" for colX, but should be "0 1 2 0 3 4";
 */
	int i = 0, j = 0, n_d, n_all;
/*
 * # positive X's in all variables should equal # positive X's in this
 */
	for (i = 0, n_all = 0; i < get_n_vars(); i++)
		for (j = 0; j < data[i]->n_X; j++)
			if (data[i]->colX[j] > 0)
				n_all++;
	for (i = 0, n_d = 0; i < d->n_X; i++)
		if (d->colX[i] > 0)
			n_d++;
	if (n_all != n_d) {
		pr_warning(
		"nr of X's in data: (%d) should match X's in other data(...) (%d)",
		n_d, n_all);
		ErrMsg(ER_IMPOSVAL, "X column definition mismatch");
	}
/*
 * now correct for 0's
 */
	for (i = 0, n_all = 0; i < get_n_vars(); i++)
		n_all += data[i]->n_X;
	if (n_all == d->n_X)
		return; /* we're done */
	n_d = d->n_X;
	d->n_X = n_all;
	d->colX = (int *) realloc(d->colX, d->n_X * sizeof(int));
	/* fill backwards */
	for (i = get_n_vars() - 1; i >= 0; i--) {
		for (j = data[i]->n_X - 1; j >= 0; j--) {
			n_all--; /* position of current X in d */
			if (data[i]->colX[j] <= 0) /* intercept, x, xy, x2, etc. */
				d->colX[n_all] = data[i]->colX[j];
			else {
				n_d--;
				if (n_d < 0)
					ErrMsg(ER_IMPOSVAL, "setup_X(): n_d < 0");
				if (d->colX[n_d] == 0)
					ErrMsg(ER_IMPOSVAL, "setup_X(): zero error");
				d->colX[n_all] = d->colX[n_d];
			}
			if (n_all < 0)
				ErrMsg(ER_IMPOSVAL, "setup_X(): n_all < 0");
		}
	}
	return;
}

METHOD get_default_method(void) {
	int i, Xset, Vgm_set;
/* 
 * no no prediction locations or no data:
 */
 	if (get_n_vars() == 0)
 		return NSP;
	if (get_n_masks() <= 0 && valdata->id < 0 &&
			gl_xvalid == 0 && data_area == NULL) {
		return UIF;
	}

/*
 * check on X variables
 */
	for (i = Xset = 0; i < get_n_vars(); i++)
		if (!(data[i]->n_X == 1 && data[i]->colX[0] == 0))
			Xset++;
/*
 * check on variograms
 */
	for (i = 0, Vgm_set = 0; i < get_n_vars(); i++)
		if (vgm[LTI(i,i)] != NULL && 
				(vgm[LTI(i,i)]->n_models > 0 || vgm[LTI(i,i)]->table != NULL)) 
					/* was: ->id >= 0*/
			Vgm_set++;
	if (!(Vgm_set == 0 || Vgm_set == get_n_vars()))
		ErrMsg(ER_SYNTAX, "set either all or no variograms");

	if (Vgm_set > 0) {
		if (get_n_beta_set() > 0)
			return SKR;
		else 
			return (Xset > 0 ? UKR : OKR);
	} else 
		return (Xset > 0 ? LSLM : IDW);
}

void set_mode(void) {
	int i, j, nm, check_failed = 0;

	if (method == NSP)
		return;
/*
 * simple, univariate:
 */
	if (get_n_vars() <= 1) {
		mode = SIMPLE;
		return;
	}
/*
 * (get_n_vars() > 1): 
 * multivariable prediction if all cross variograms set parameters merge
 */
	for (i = check_failed = 0; i < get_n_vars(); i++)
		for (j = 0; j < i; j++)
			if (vgm[LTI(i,j)] == NULL || vgm[LTI(i,j)]->id < 0)
				check_failed = 1;
	if (check_failed == 0) {
		mode = MULTIVARIABLE;
		return;
	}
	if (n_variograms_set() == 0) {
		for (i = 0; i < get_n_vars(); i++)
			if (data[i]->n_merge > 0) {
				mode = MULTIVARIABLE;
				return;
			}
	}

/*
 * stratify? ONLY if: 
 * 0. get_n_vars() > 1; no cross variograms set ==>> has been checked.
 * 1. no pred(): or var(): except for first variable;
 * 2. No masks and valdata->what_is_u == U_ISSTRATUM
 * 3. mask is a valid strata map, n categories > 1
 */
	if (n_masks == 0) {
		mode = (valdata->what_is_u == U_ISSTRATUM) ? STRATIFY : SIMPLE;
		return;
	} 
	for (i = 2, check_failed = 0; i < get_n_outfile(); i++)
		if (get_outfile_namei(i) != NULL)
			check_failed = 1;
	if (! check_failed) {
		for (i = nm = 0; i < get_n_vars(); i++)
			for (j = 0; j < data[i]->n_X; j++)
				/* check for user-defined base function values in mask map */
				if (data[i]->colX[j] > 0)
					nm++; /* UK or lm */
		/* if STRATIFIED: first mask categories, then base functions */
		check_failed = (1 + nm != n_masks);
	}
#ifndef USING_R
	if (! check_failed)
		check_failed = !is_valid_strata_map(get_mask_name(0), get_n_vars());
#endif
	if (check_failed) {
		mode = SIMPLE;
		return;
	} else {
		mode = STRATIFY;
		return;
	}
}

int decide_on_coincide(void) {
	int i, j;

	if (get_n_vars() <= 1)
		return 0;
	if (get_mode() == STRATIFY)
		return 0; /* data may coincide, but prediction locations won't */
	for (i = 1; i < get_n_vars(); i++) {
		if (data[i]->n_list != data[0]->n_list ||
				data[i]->colnx != data[0]->colnx ||
				data[i]->colny != data[0]->colny ||
				data[i]->colnz != data[0]->colnz ||
				data[i]->sel_min != data[0]->sel_min ||
				data[i]->sel_max != data[0]->sel_max ||
				data[i]->force != data[0]->force ||
				data[i]->sel_rad != data[0]->sel_rad)
			return 0; /* don't check filename: 'file.dat' and './file.dat' */
/*
 * Consider the data file
 * 1 2 NA 3
 * 2 3 4 NA
 * with x=1, y=2, v=3 (for var 1), v=4 (for var 2).
 * This is only distinguishable by doing it the ``hard way'':
 */
		for (j = 0; j < data[0]->n_list; j++) {
			if (data[0]->list[j]->x != data[i]->list[j]->x ||
				data[0]->list[j]->y != data[i]->list[j]->y ||
				data[0]->list[j]->z != data[i]->list[j]->z)
				return 0;
		}
	}
	if (DEBUG_DUMP)
		printlog("(identical search conditions found for all variables)\n");
	return 1;
}

void check_global_variables(void) {
/*
 * Purpose       : check internal variable consistency, add some parameters 
 * Created by    : Edzer J. Pebesma
 * Date          : april 13, 1992
 * Prerequisites : none
 * Returns       : -
 * Side effects  : none
 * also check Cauchy-Schwartz unequality on cross/variograms.
 */
	int i, j, nposX, n_merge = 0;
	METHOD m;
	VARIOGRAM *v_tmp;

	if (!almost_equals(gl_plotfile, DEF_plotfile))
		plotfile = efopen(gl_plotfile, "w");
	/* UK: check if n_masks equals total nr of unbiasedness cond. */
	if (gl_nblockdiscr < 2)
		ErrMsg(ER_RANGE, "nblockdiscr must be >= 2");
	if (method == UKR || method == LSLM) {
		nposX = 0;
		for (i = 0; i < get_n_vars(); i++)
			for (j = 0; j < data[i]->n_X; j++) {
				if (data[i]->colX[j] > 0)
					nposX++;
			}
		if (n_masks > 0 && nposX > 0 &&
				nposX + (get_mode() == STRATIFY ? 1 : 0) != n_masks) {
			message("nposX %d n_masks %d stratify %d\n", nposX, n_masks,
				(get_mode() == STRATIFY ? 1 : 0));

			ErrMsg(ER_IMPOSVAL,
			"number of masks does not equal number of (positive) X's");
		}
	}
    if (method == POLY) {
        if (valdata->id > -1) {
            message("Do not use dummy data() with point-in-polygon method,\n");
            message("Use the only data(...) command.\n");
        }
        
        if (data_area) message("area: ... command is ignored with point-in-polygon method\n"); 
        if ( get_n_vars()>2) {
            message("More then one data(...) command cannot be used with point-in-polygon method,\n");
            message("use the only data(...) command.\n");
            message("Only first data(...) specification will be used.\n");
        } else if ( get_n_vars() == 0) {
            message("Specify the data you want to test with data(...) statement.\n");
            ErrMsg(ER_SYNTAX, "data specification error");
        }
        if (get_n_edges() == 0) {
            message("Point-in-polygon method required, but no edges data is given.\n");
            ErrMsg(ER_SYNTAX, "data specification error");
        }
        if (!(data[0]->mode & X_BIT_SET))
			ErrMsg(ER_VARNOTSET, "x coordinate not set");
		if (!(data[0]->mode & Y_BIT_SET))
			ErrMsg(ER_VARNOTSET, "y coordinate not set");
		
        /* return; */
    }
    
	if (method == SPREAD) {
		for (i = 0; i < get_n_vars(); i++)
			if (data[i]->sel_rad == DBL_MAX)
				data[i]->sel_rad *= 0.99; /* force distance calculation */
	}

	if (get_n_beta_set() != 0 && get_n_beta_set() != get_n_vars())
		ErrMsg(ER_SYNTAX, 
			"set sk_mean or beta either for all or for no variables");

	if (!(method == ISI || method == GSI)) {
		if (gl_nsim > 1)
			ErrMsg(ER_IMPOSVAL, "nsim only allowed for simulation");
	}

	if (method == ISI && max_block_dimension(0) > 0.0)
		ErrMsg(ER_IMPOSVAL, "indicator simulation only for points");
	/*
	 * check if both block and area are set
	 */
	if (data_area != NULL && (block.x > 0.0 || block.y > 0.0 || block.z > 0.0))
		ErrMsg(ER_IMPOSVAL, "both block and area set: choose one");
	/*
	 * check for equality of coordinate dimensions:
	 */
	for (i = 1; i < get_n_vars(); i++) {
		if ((data[i]->mode & V_BIT_SET) != (data[0]->mode & V_BIT_SET))  {
			message("data(%s) and data(%s):\n", name_identifier(0), 
				name_identifier(i));
			ErrMsg(ER_IMPOSVAL, "data have different coordinate dimensions");
		}
	}
	if (valdata->id > 0 && data[0]->dummy == 0 &&
			((data[0]->mode | (V_BIT_SET | S_BIT_SET)) !=
			 (valdata->mode | (V_BIT_SET | S_BIT_SET)))) {
		message("data() and data(%s):\n", name_identifier(0));
		ErrMsg(ER_IMPOSVAL, "data have different coordinate dimensions");
		for (i = 0; i < get_n_vars(); i++) {
			if (data[i]->dummy) {
				data[i]->mode = (valdata->mode | V_BIT_SET);
				data[i]->minX = valdata->minX;
				data[i]->minY = valdata->minY;
				data[i]->minZ = valdata->minZ;
				data[i]->maxX = valdata->maxX;
				data[i]->maxY = valdata->maxY;
				data[i]->maxZ = valdata->maxZ;
				set_norm_fns(data[i]);
			}
		}
	}
	if (get_n_masks()) {
		if (! get_n_vars() > 0)
			ErrMsg(ER_VARNOTSET, "no data");
		if (!(data[0]->mode & X_BIT_SET))
			ErrMsg(ER_VARNOTSET, "x coordinate not set");
		if (!(data[0]->mode & Y_BIT_SET))
			ErrMsg(ER_VARNOTSET, "y coordinate not set");
		if ((data[0]->mode & Z_BIT_SET) && is_mv_double(&gl_zmap))
			ErrMsg(ER_VARNOTSET, "parameter zmap must be set with 3D data");
	}

	for (i = 0; i < get_n_vars(); i++) {
		if (data[i]->fname == NULL && !data[i]->dummy) {
			message("file name for data(%s) not set\n", name_identifier(i));
			ErrMsg(ER_VARNOTSET, (const char *) command_file_name);
		}
		if (data[i]->id < 0) {
			message("data(%s) not set\n", name_identifier(i));
			ErrMsg(ER_VARNOTSET, (const char *) command_file_name);
		}
		if (data[i]->beta && data[i]->beta->size != data[i]->n_X) {
			pr_warning("beta dimension (%d) should equal n_X (%d)", 
				data[i]->beta->size, data[i]->n_X);
			ErrMsg(ER_IMPOSVAL, "sizes of beta and X don't match");
		}
		if (data[i]->sel_rad == DBL_MAX && data[i]->oct_max > 0)
			ErrMsg(ER_IMPOSVAL, 
				"define maximum search radius (rad) for octant search");
#ifndef USING_R
		if (data[i]->vdist && (vgm[LTI(i,i)] == NULL 
				|| vgm[LTI(i,i)]->n_models <= 0))
			ErrMsg(ER_IMPOSVAL, "define direct variogram models to use vdist");
#endif
		if (data[i]->vdist && data[i]->sel_rad == DBL_MAX)
			ErrMsg(ER_IMPOSVAL, "when using vdist, radius should be set");
		if (! data[i]->dummy && ! (data[i]->mode & V_BIT_SET)) {
			message("no v attribute set for data(%s)\n", 
				name_identifier(data[i]->id));
			ErrMsg(ER_VARNOTSET, (const char *) command_file_name);
		}
		if (method != SEM && method != COV && method != MAPVALUE) {
			/* check neighbourhood settings */
			if (data[i]->sel_rad < 0.0 || data[i]->sel_min < 0 || 
				data[i]->sel_max < 0 || (data[i]->sel_min > data[i]->sel_max)) {
				message(
				"invalid neighbourhood selection: radius %g max %d min %d\n", 
				data[i]->sel_rad, data[i]->sel_max, data[i]->sel_min);
				ErrMsg(ER_IMPOSVAL, (const char *) command_file_name);
			}
		}
		if (data[i]->id > -1 && (method == OKR || method == SKR || 
				is_simulation(method) || method == UKR)) {
			if (vgm[LTI(i,i)] == NULL || vgm[LTI(i,i)]->id < 0) {
				message("variogram(%s) not set\n", name_identifier(i));
				ErrMsg(ER_VARNOTSET, "variogram()");
			}
		}
		n_merge += data[i]->n_merge;
	}
	if (n_merge && get_mode() != MULTIVARIABLE)
		ErrMsg(ER_IMPOSVAL, "merge only works in multivariable mode");
	if (mode == SIMPLE && get_method() != UIF) {  /* check if it's clean: */
		for (i = 0; i < get_n_vars(); i++)
			for (j = 0; j < i; j++)
				if (vgm[LTI(i,j)] != NULL && vgm[LTI(i,j)]->id > 0) {
					message("variogram(%s, %s) %s\n", name_identifier(i),
						name_identifier(j),
					"can only be set for ck, cs, uk, sk, ok, sem or cov");
					ErrMsg(ER_IMPOSVAL, "variogram()");
				}
	} 
	if ((m = get_default_method()) != get_method()) {
		if (m == UKR && (get_method() == OKR || get_method() == SKR))
			ErrMsg(ER_IMPOSVAL,
				"\nremove X=... settings for ordinary or simple kriging");
		if (m == OKR && get_method() == SKR)
			ErrMsg(ER_IMPOSVAL, "method: something's terribly wrong!");
		if (m == OKR && get_method() == UKR) {
			message("I would recommend:\n");
			message("Do not specify uk if ok is all you'll get\n");
		}
	}
	if (mode == MULTIVARIABLE && get_method() != UIF && get_method() != SEM
			&& get_method() != COV && n_variograms_set() > 0)
		check_variography((const VARIOGRAM **) vgm, get_n_vars());
	if (gl_n_marginals > 0) {
		if (gl_marginal_names) {
			if (get_mode() == STRATIFY && gl_n_marginals != 2)
				ErrMsg(ER_IMPOSVAL, "# marginal maps should be 2 (stratify)");
			if (get_mode() != STRATIFY && gl_n_marginals != 2 * get_n_vars())
				ErrMsg(ER_IMPOSVAL, "# marginals should be 2 x # of variables");
		} else if (gl_n_marginals != 2 * get_n_vars())
			ErrMsg(ER_IMPOSVAL, "# marginals should be 2 x # of variables");
	}
	v_tmp = init_variogram(NULL);
#ifndef USING_R
	for (i = 0; i < get_n_vgms(); i++) {
		if (vgm[i]) {
			if (vgm[i]->fname && read_variogram(v_tmp, vgm[i]->fname) == 0)
				pr_warning("NOTE that `%s' is a file name, variogram models never have quotes", vgm[i]->fname);
			if (vgm[i]->fname2 && read_variogram(v_tmp, vgm[i]->fname2) == 0)
				pr_warning("NOTE that `%s' is a file name, variogram models never have quotes", vgm[i]->fname2);
		}
	}
#endif
	free_variogram(v_tmp);
} 

void remove_all(void) {

	while (n_vars)
		remove_id(0); /* hard way */
	/* for (i = n_vars-1; i >= 0; i--)
		remove_id(i);
	*/
	/* remove_id(n_vars - 1);  */
	/* the hard way; remove_id(n_vars-1) would be the ``easy'' alternative */
	gls(NULL, 0, GLS_INIT, NULL, NULL); /* cleans up static arrays */
	what_is_outfile(-1); /* cleans up static array */
	reset_block_discr(); /* resets block settings */
	max_block_dimension(1); /* reset */
	if (gl_bounds != NULL) {
		efree(gl_bounds);
		gl_bounds = NULL;
	}
	if (valdata != NULL)
		free_data(valdata);
	valdata = NULL;
}

int remove_id(const int id) {
/*
 * remove id id, and reset data, vgm, ids, outfile_names
 */
	int i, j, id_new, id_old;
	VARIOGRAM *vp;

	assert(id >= 0 && id < n_vars);


	/* reset data */
	free_data(data[id]);
	data[id] = NULL;
	for (i = id; i < n_vars - 1; i++) {
		data[i] = data[i+1];
		data[i]->id = i;
	}

	for (i = 0; i < n_vars; i++) {
		j = LTI(i,id);
		if (vgm[j]) {
			free_variogram(vgm[j]);
			vgm[j] = NULL;
		}
	}

	/* copy variograms: */
	for (i = id; i < n_vars - 1; i++) {
		for (j = id; j <= i; j++) {
			id_new = LTI(i,j);
			id_old = LTI(i+1,j+1);
			vp = vgm[id_new] = vgm[id_old];
			if ((vp != NULL) && (vp->id1 >= 0 || vp->id2 >= 0)) {
				vp->id1 = i;
				vp->id2 = j;
				vp->id = id_new;
			}
		}
	}

	/* reset identifiers: */
	efree(ids[id]);
	for (i = id; i < n_vars - 1; i++)
		ids[i] = ids[i+1];

	/* free outfilenames */
	if (outfile_names[2 * id]) {
		efree(outfile_names[2 * id]);
		outfile_names[2 * id] = NULL;
	}
	if (outfile_names[2 * id + 1]) {
		efree(outfile_names[2 * id + 1]);
		outfile_names[2 * id + 1] = NULL;
	}

	/* shift pred(xx)/variances(xx) names: */
	for (i = id; i < n_vars - 1; i++) {
		outfile_names[2 * i] = outfile_names[2 * (i + 1)];
		outfile_names[2 * i + 1] = outfile_names[2 * (i + 1) + 1];
	}

	/* shift covariances(xx): */
        
	for (i = id; i < n_vars - 1; i++) {
		id_old = 2 *  n_vars + LTI2(i,id);
		if (outfile_names[id_old]) {
			efree(outfile_names[id_old]);
			outfile_names[id_old] = NULL;
		}
		for (j = id; j < i; j++) {
			id_new = 2 * (n_vars - 1) + LTI2(i,j);
			id_old = 2 * n_vars + LTI2(i+1,j+1);
			outfile_names[id_new] = outfile_names[id_old];
		}
	}

	n_vars -= 1;

	if (n_vars == 0)
		clean_up();

	init_gstat_data(n_vars); /* reset sizes */
	return n_vars;
}

static void clean_up(void) {
	int i;

	/* free variograms */
	if (vgm) {
		efree(vgm);
		vgm = NULL;
	}

	/* free data */
	if (data) {
		efree(data);
		data = NULL;
	}

	/* free valdata */
	if (valdata) {
		free_data(valdata);
		valdata = NULL;
	}

	/* free data_area */
	if (data_area) {
		free_data(data_area);
		data_area = NULL;
	}

	/* free outfile names */
	if (outfile_names) {
		efree(outfile_names);
		outfile_names = NULL;
	}

	/* free identifiers */
	if (ids) {
		efree(ids);
		ids = NULL;
	}

	/* free mask names */
	if (mask_names) {
		for (i = 0; i < n_masks; i++) {
			if (mask_names[i])
				efree(mask_names[i]);
		}
		efree(mask_names);
		mask_names = NULL;
	}
	n_vars = 0;
	n_masks = 0;
	n_last = 0;
	n_v_last = 0;
	n_o_last = 0;
	mode = MODE_NSP;
}

void push_edges_name(const char *name) {
	if (edges_names == NULL)
		edges_names = (char **) emalloc((n_edges + 1) * sizeof(char *));
	else
		edges_names = (char **) erealloc(edges_names,
				(n_edges + 1) * sizeof(char *));
	edges_names[n_edges] = string_dup(name);
	n_edges++;
}

const char *get_edges_name(const int i) {
	if (i >= n_edges || i < 0)
		return NULL;
	return edges_names[i];
}

int *get_n_edges_polys(void) {

	assert(n_edges_polys != NULL);
	return n_edges_polys;
}

int *set_n_edges_polys(int *n) {
	return n_edges_polys=n;
}

int get_n_edges(void) {
	return n_edges;
}

POLYGON **get_edges(void) {
	return edges;
}

POLYGON **set_edges(POLYGON **new_edges) {
    return edges=new_edges;
}

DATA *setup_poly_method(void) {
/* assuming data[0] exists! */
	assert(data != NULL && data[0] != NULL);

    
	if (get_n_masks() == 0) {
		if (valdata) 
			free_data(valdata);
		valdata = data[0];
	}

	return valdata; /* why? */
}
