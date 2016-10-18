/*
 * msim.c: multiple simulation database + output
 * Written during my Post Doc at UvA, end of 1996.
 * Rewrite started on Sat Apr 10 20:54:05 WET DST 1999
 */
#include <stdio.h>
#include <stdlib.h> /* qsort() */
#include <string.h> /* memmove() */
#include <math.h>   /* sqrt(), ... */

#include "defs.h"

#include <R.h>

#include "debug.h"
#include "data.h"
#include "utils.h"
#include "vario.h"
#include "glvars.h" /* gl_nsim, gl_format */
#include "userio.h"
#include "mapio.h"
#include "mtrx.h"
#include "lm.h"
#include "gls.h"
#include "sim.h"
#include "msim.h"

static DPOINT *which_point(DATA *d, DPOINT *where);
static unsigned int *get_n_sim_locs_table(unsigned int *size);

/* global variables formerly in predict.c; set in s.c */
unsigned int n_pred_locs = 0;

#ifdef SIM_DOUBLE
typedef double Float; /* doubles the memory requirement -> may be pretty much */
#else
typedef float Float;
#endif

static Float 
	***msim = NULL, 
	/* 
	 * The msim table entry for variable i, for simulation node j, 
   	 * replicate k is msim[i][j][k] 
   	 */
	**msim_base = NULL; /* base structure for blocked allocation */
static double
	***beta = NULL;
	/* 
	 * the beta realisation for variable i, draw j is in
	 * beta[i][j] -- and has dimension data[i]->beta->size
	 */
static unsigned int 
	*n_sim_locs = NULL, /* n simulation locations per data variable */
	table_size = 0, /* offset strata table size */
	**s2d = NULL, 
	/*
	 * s2d: <data list entries> -- find (POINT *) from msim:
	 * msim[i][j][...] -->> data[i]->list[s2d[i][j]]
	 */
	**d2s = NULL;
	/* 
	 * d2s: <msim entries> find msim entry from (POINT *):
	 * data[i]->list[j] -->> msim[i][d2s[i][j]][...]
	 * ((the latter two are necessary because simple counting fails
	 * when point simulation locations coincide with data locations. In
	 * this case, a data DPOINT is not added to the data list, and so we
	 * need to keep track _where_ simulated values were located in order
	 * to write them back to output (file/maps) at the end))
	 */

void print_sim(void) {
/* print complete contents of sim_values -- for debug purposes only */
	int i, j, k;
	for (i = 0; i < get_n_vars(); i++) {
		printlog("variable %d:\n", i);
		for (j = 0; j < n_sim_locs[i]; j++) {
			for (k = 0; k < gl_nsim; k++)
				printlog(" %g", msim[i][j][k]);
			printlog("\n");
		}
	}
}

void init_simulations(DATA **d) {
	int i, j, size;

	assert(n_pred_locs > 0); /* should be set by now... */

	if (msim != NULL) 
		free_simulations();

	n_sim_locs = get_n_sim_locs_table(&table_size);

	if (DEBUG_DUMP) {
		printlog("n_sim_locs_table: ");
		for (i = 0; i < table_size; i++)
			printlog("[%d] ", n_sim_locs[i]);
		printlog("\n");
	}

	msim = (Float ***) emalloc(get_n_vars() * sizeof(Float **));
	msim_base = (Float **) emalloc(get_n_vars() * sizeof(Float *));
	s2d = (unsigned int **) emalloc(get_n_vars() * sizeof(unsigned int *));
	d2s = (unsigned int **) emalloc(get_n_vars() * sizeof(unsigned int *));
	for (i = 0; i < get_n_vars(); i++) {
		/* msim stuff: */
		size = n_sim_locs[i] * gl_nsim;
		msim_base[i] = (Float *) emalloc(size * sizeof(Float));
		memset(msim_base[i], 0xFF, size * sizeof(Float));
		msim[i] = (Float **) emalloc(n_sim_locs[i] * sizeof(Float *));
		for (j = 0; j < n_sim_locs[i]; j++)
			msim[i][j] = &(msim_base[i][j * gl_nsim]);

		/* index stuff: */
		s2d[i] = (unsigned int *) emalloc(n_sim_locs[i] * sizeof(unsigned int));
		d2s[i] = (unsigned int *) emalloc(n_sim_locs[i] * sizeof(unsigned int));
		/* now let's trigger some Seg Faults if on error: */
		memset(s2d[i], 0xFF, n_sim_locs[i] * sizeof(unsigned int));
		memset(d2s[i], 0xFF, n_sim_locs[i] * sizeof(unsigned int));
	}
}

void save_sim(DATA **data, DPOINT *where, int sim, int n_vars,
	const double *value, int *is_pt) {
/*
 * save the last simulated value(s) in msim;
 * data[0]->n_list and data[0]->n_original (mode != STRATIFY) or
 * else n_vars (..) denote where it should go.
 */
	int i, row;
	DPOINT *which = NULL;

	if (gl_nsim <= 1)
		return;

	for (i = 0; i < n_vars; i++) { /* store current simulation */
		row = data[i]->n_list - data[i]->n_original + data[i]->nsim_at_data;
		if (sim == 0) { /* fill d2s and s2d entries: */
			assert(row >= 0 && row < n_sim_locs[i]); 
			if (is_pt[i]) {
				which = which_point(data[i], where);
				s2d[i][row] = GET_INDEX(which);
				/* d2s remains MV */
			} else { /* newly simulated */
				s2d[i][row] = data[i]->n_list; 
				d2s[i][data[i]->n_list - data[i]->n_original] = row;
				/* please go and check it -- this last line took me
				4 hours to get right */
			}
		} 
		msim[i][row][sim] = value[i];
	}
}

void save_sim_strat(DATA *d, DPOINT *where, int sim, double value, int is_pt) {
/* the same, but for stratified mode */
	int row;
	DPOINT *which = NULL;

	if (gl_nsim <= 1)
		return;

	row = d->n_list - d->n_original + d->nsim_at_data;
	if (sim == 0) { /* fill d2s and s2d entries: */
		assert(row >= 0 && row < n_sim_locs[d->id]);
		if (is_pt) {
			which = which_point(d, where);
			s2d[d->id][row] = GET_INDEX(which);
			/* the entry for d2s does not exist */
		} else { /* not an original point, but a simulated one: */
			s2d[d->id][row] = d->n_list; /* newly simulated */
			d2s[d->id][d->n_list - d->n_original] = row;
		}
	}
	msim[d->id][row][sim] = value;
}

static DPOINT *which_point(DATA *d, DPOINT *where) {
	int i;
	double dzero2;

#define WPWARNING "if you are simulating with a Gaussian variogram model without nugget\n\
then try to add a small nugget variance to avoid the following error message"

	dzero2 = gl_zero * gl_zero;
	for (i = 0; i < d->n_sel; i++)
		if (fabs(d->pp_norm2(d->sel[i], where)) <= dzero2)
			return d->sel[i];
	pr_warning(WPWARNING);
	ErrMsg(ER_NULL, "which_point(): point not found");
	return where; /* never reached */
}

void restore_data_sel(DATA **data, int sim, int n_vars) {
	unsigned int i, j;
	int id, idx;
	DATA *d0 = NULL;

	if (gl_nsim <= 1)
		return;

	if (n_vars == 0) {
		assert(get_mode() == STRATIFY);
		d0 = data[0];
		for (j = 0; j < d0->n_sel; j++) {
			id = d0->id;
			idx = GET_INDEX(d0->sel[j]) - d0->n_original;
			/* current sim: */
			if (idx >= 0 && d2s[id][idx] != 0xFFFFFFFF)
				d0->sel[j]->attr = msim[id][d2s[id][idx]][sim]; 
		}
	} else {
		for (i = 0; i < n_vars; i++) {
			for (j = 0; j < data[i]->n_sel; j++) { 
				idx = GET_INDEX(data[i]->sel[j]) - data[i]->n_original;
				/* idx < 0 -->> original data, don't restore */
				if (idx >= 0 && d2s[i][idx] != 0xFFFFFFFF)
					data[i]->sel[j]->attr = msim[i][d2s[i][idx]][sim];
	 				/* data[i]->list[j] -->> msim[i][d2s[i][j]][...] */
			}
		}
	}
}

void free_simulations(void) {
	int i, j;

	if (msim != NULL) {
		for (i = 0; i < get_n_vars(); i++) {
			efree(msim[i]);
			efree(msim_base[i]);
			efree(s2d[i]);
			efree(d2s[i]);
		}
		efree(msim);
		msim = NULL;
		efree(msim_base);
		msim_base = NULL;
	}
	if (s2d != NULL) {
		efree(s2d);
		s2d = NULL;
	}
	if (d2s != NULL) {
		efree(d2s);
		d2s = NULL;
	}
	if (beta != NULL) {
		for (i = 0; i < get_n_vars(); i++) {
			for (j = 0; j < gl_nsim; j++)
				efree(beta[i][j]);
			efree(beta[i]);
		}
		efree(beta);
		beta = NULL;
	}
	if (n_sim_locs != NULL)
		free(n_sim_locs);
	n_sim_locs = NULL;
}

void setup_beta(DATA **d, int n_vars, int n_sim) {
	double *est;
	const double *sim = NULL;
	int i, j, k, sum_n_X = 0, offset, *is_pt = NULL;

	assert(beta == NULL); /* entered only once */

	/* allocate beta */
	beta = (double ***) emalloc(n_vars * sizeof(double **));
	for (i = 0; i < n_vars; i++) {
		assert(d[i]->n_list > 0);
		beta[i] = (double **) emalloc(n_sim * sizeof(double *));
		for (j = 0; j < n_sim; j++)
			beta[i][j] = (double *) emalloc(d[i]->n_X * sizeof(double));
	}
	for (i = 0; i < n_vars; i++) {
		if (d[i]->beta == NULL) /* push bogus values */
			for (j = 0; j < d[i]->n_X; j++)
				d[i]->beta = push_d_vector(-9999.9, d[i]->beta);
		sum_n_X += d[i]->n_X;
	}
	printlog("drawing %d %s%s realisation%s of beta...\n", n_sim,
		n_vars > 1 ? (gl_sim_beta == 0 ? "multivariate " : "univariate ") : "",
		gl_sim_beta == 2 ? "OLS" : "GLS",
		n_sim > 1 ? "s" : "");
	is_pt = (int *) emalloc(sum_n_X * sizeof(int));
	if (gl_sim_beta == 0) {
		est = make_gls_mv(d, n_vars);
		for (j = 0; j < n_sim; j++) {
			sim = cond_sim(est, sum_n_X, GSI, is_pt, 0); /* length sum_n_X */
			for (i = offset = 0; i < n_vars; i++) {
				for (k = 0; k < d[i]->n_X; k++)
					beta[i][j][k] = sim[offset + k];
				offset += d[i]->n_X;
				if (DEBUG_DUMP || DEBUG_COV) {
					printlog("var=%d, sim=%d, beta=[ ", i, j);
					for (k = 0; k < d[i]->n_X; k++)
						printlog("%g ", beta[i][j][k]);
					printlog("]\n");
				}
			}
		}
		efree(est);
	} else {
		for (i = 0; i < n_vars; i++) {
			if (gl_sim_beta == 1)
				est = make_gls(d[i], 0);
			else /* gl_sim_beta == 2 */
				est = make_ols(d[i]);
			for (j = 0; j < n_sim; j++) {
				sim = cond_sim(est, d[i]->n_X, GSI, is_pt, 0);
				for (k = 0; k < d[i]->n_X; k++)
					beta[i][j][k] = sim[k];
				if (DEBUG_DUMP || DEBUG_COV) {
					printlog("var=%d, sim=%d, beta=[ ", i, j);
					for (k = 0; k < d[i]->n_X; k++)
						printlog("%g ", beta[i][j][k]);
					printlog("]\n");
				}
			}
			efree(est);
		}
	}
	efree(is_pt);
	return;
}

void set_beta(DATA **d, int sim, int n_vars, METHOD method) {
/* n_vars == 0 --> stratified mode, check d[0]->id to find out which
 * data we've got here
 */
 	int i;

 	assert(d[0]->beta);
	if (beta == NULL) /* use the values entered by the user */
		return;
 	
 	if (get_mode() != STRATIFY) {
 		for (i = 0; i < n_vars; i++)
 			d[i]->beta->val = beta[i][sim];
 	} else
 		d[0]->beta->val = beta[d[0]->id][sim];
 	return;
}

float ***get_msim(void) {
	return msim;
}

static unsigned int *get_n_sim_locs_table(unsigned int *size) {
	unsigned int i, *table;

	*size = (int) get_n_vars();
	table = (unsigned int *) emalloc(*size * sizeof(int));
	for (i = 0; i < *size; i++)
		table[i] = n_pred_locs;
	return table;
}
