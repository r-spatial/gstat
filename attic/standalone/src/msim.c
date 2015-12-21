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
 * msim.c: multiple simulation database + output
 * Written during my Post Doc at UvA, end of 1996.
 * Rewrite started on Sat Apr 10 20:54:05 WET DST 1999
 */
#include <stdio.h>
#include <stdlib.h> /* qsort() */
#include <string.h> /* memmove() */
#include <math.h>   /* sqrt(), ... */

#include "defs.h"

#ifdef USING_R
# include <R.h>
#endif

#include "debug.h"
#include "random.h"
#include "data.h"
#include "utils.h"
#include "vario.h"
#include "glvars.h" /* gl_nsim, gl_format */
#include "predict.h" /* n_pred_locs */
#include "userio.h"
#include "mapio.h"
#include "lm.h"
#include "gls.h"
#include "sim.h"
#include "msim.h"

static DPOINT *which_point(DATA *d, DPOINT *where);
static void restore_data_list(DATA **data, int sim, int n_vars);
static void lhs_one(Double_index *list, unsigned int dim, double mean, double var);

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
		printf("\n");
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

static void restore_data_list(DATA **data, int sim, int n_vars) {
/*
 * for data[0]..data[n_vars-1] restore the values of the `sim' simulation
 * back into the attribute values
 */
	int i, j, id, idx;
	DATA *d0;

	assert(gl_nsim > 1);

	if (get_mode() == STRATIFY) {
	/*
	 * n_vars is in this case n_pred_locs, so the loop is over all locations:
	 */
		d0 = data[0];
		for (j = d0->n_original; j < d0->n_list; j++) {
			id = d0->id;
			idx = j - d0->n_original;
			if (d2s[id][idx] != 0xFFFFFFFF)
				d0->list[j]->attr = msim[id][d2s[id][idx]][sim];
		}
	} else {
		for (i = 0; i < n_vars; i++)
			for (j = data[i]->n_original; j < data[i]->n_list; j++) {
				idx = j - data[i]->n_original;
				if (d2s[i][idx] != 0xFFFFFFFF)
					data[i]->list[j]->attr = msim[i][d2s[i][idx]][sim];
			}
	}
}

#ifndef USING_R
void save_simulations_to_ascii(const char *fname) {
	DATA **d, *dv;
	FILE *f;
	int i, j, k, n_vars, n_cols = 0;

	if (msim == NULL || gl_nsim <= 1)
		return;

#ifdef DEBUGSIM
	print_sim();
#endif
	f = efopen(fname, "w");
	d = get_gstat_data();
	dv = get_dataval();
	n_vars = get_n_vars();
	if (d[0]->mode & X_BIT_SET) n_cols++;
	if (d[0]->mode & Y_BIT_SET) n_cols++;
	if (d[0]->mode & Z_BIT_SET) n_cols++;

	if (get_mode() == STRATIFY)
		n_cols += (1 + gl_nsim);
	else
		n_cols += n_vars * gl_nsim;

	if (dv->type.type == DATA_EAS) {
		fprintf(f, "#%s simulation results (seed %lu)\n%d\n",
			get_mode() == STRATIFY ?  "stratified" : "multiple", get_seed(), n_cols);
		if (dv->mode & X_BIT_SET) fprintf(f, "%s\n", NULS(dv->x_coord));
		if (dv->mode & Y_BIT_SET) fprintf(f, "%s\n", NULS(dv->y_coord));
		if (dv->mode & Z_BIT_SET) fprintf(f, "%s\n", NULS(dv->z_coord));
		if (get_mode() == STRATIFY)
			fprintf(f, "stratum\n");
	}

	if (get_mode() != STRATIFY) {
		for (i = 0; dv->type.type == DATA_EAS && i < n_vars; i++) /* header */
			for (j = 0; j < gl_nsim; j++)
				fprintf(f, "%s--sim#%d\n", name_identifier(i), j);
		for (i = 0; i < n_sim_locs[0]; i++) {
            if (d[0]->mode & X_BIT_SET) {
                fprintf(f, " ");
                fprintf(f, gl_format, d[0]->list[s2d[0][i]]->x);
            }

			if (d[0]->mode & Y_BIT_SET) {
                fprintf(f, " ");
                fprintf(f, gl_format, d[0]->list[s2d[0][i]]->y);
            }

			if (d[0]->mode & Z_BIT_SET) {
                fprintf(f, " ");
                fprintf(f, gl_format, d[0]->list[s2d[0][i]]->z);
            }

			for (j = 0; j < n_vars; j++) {
				for (k = 0; k < gl_nsim; k++) {
                    fprintf(f, " ");
    				fprintf(f, gl_format, msim[j][i][k]);
                }
			}

			fprintf(f, "\n");
		}
	} else { /* STRATIFIED */
		for (j = 0; dv->type.type == DATA_EAS && j < gl_nsim; j++)
			fprintf(f, "stratified sim#%d\n", j);
		for (i = 0; i < n_vars; i++) { 
			for (j = 0; j < n_sim_locs[i]; j++) { 
				if (d[i]->mode & X_BIT_SET) {
                    fprintf(f, " ");
                    fprintf(f, gl_format, d[i]->list[s2d[i][j]]->x);
                }

				if (d[i]->mode & Y_BIT_SET) {
                    fprintf(f, " ");
                    fprintf(f, gl_format, d[i]->list[s2d[i][j]]->y);
                }

				if (d[i]->mode & Z_BIT_SET) {
                    fprintf(f, " ");
                    fprintf(f, gl_format, d[i]->list[s2d[i][j]]->z);
                }

				if (d[i]->mode & S_BIT_SET)
					fprintf(f, " %d", d[i]->list[s2d[i][j]]->u.stratum + strata_min);
				else
					fprintf(f, " %d", d[i]->id + strata_min);

				for (k = 0; k < gl_nsim; k++) {
                    fprintf(f, " ");
					fprintf(f, gl_format, msim[i][j][k]);
                }

				fprintf(f, "\n");
			}
		}
	}
	free_simulations();
	efclose(f);
}

void save_simulations_to_maps(GRIDMAP *mask) {
	unsigned int i, j, k, r, c;
	char *name = NULL, *what = NULL;
	const char *cp;
	DATA **d;
	GRIDMAP *m;
	void map_sign(GRIDMAP *m, const char *what);

	if (msim == NULL || gl_nsim <= 1)
		return;

#ifdef DEBUGSIM
	print_sim();
#endif
	d = get_gstat_data();
	/* remove standard output files: they equal the last one  ********
	for (i = 0; i < get_n_vars(); i++)
		if (get_outfile_namei(2 * i) && remove(get_outfile_namei(2 * i)))
			pr_warning("could not remove %s", get_outfile_namei(2 * i));
	********************/
	for (i = 0; i < gl_nsim; i++) {
		if (get_mode() != STRATIFY) {
			restore_data_list(d, i, get_n_vars());
			for (j = 0; j < get_n_vars(); j++) {
				if ((cp = get_outfile_namei(2 * j)) != NULL) {
					name = (char *) erealloc(name, 
							(strlen(cp) + 12) * sizeof(char));
					map_name_nr(mask, cp, name, i , gl_nsim);
					what = (char *) erealloc(what, 
							(strlen(name_identifier(j)) + 20) * sizeof(char));
					sprintf(what, "%s : simulation %d", name_identifier(j), i);
					m = map_dup(name, mask);

					for (k = 0; k < n_sim_locs[j]; k ++) {
						assert(s2d[j][k] < d[j]->n_list);
						if (map_xy2rowcol(m, 
								d[j]->list[s2d[j][k]]->x,
								d[j]->list[s2d[j][k]]->y, 
								&r, &c)) {
							pr_warning("x %g y %g", 
								d[j]->list[s2d[j][k]]->x,
								d[j]->list[s2d[j][k]]->y);
							ErrMsg(ER_IMPOSVAL, "map_xy2rowcol");
						}
						map_put_cell(m, r, c, d[j]->list[s2d[j][k]]->attr);
					}
					map_sign(m, what);
					m->write(m);
					map_free(m);
				}
			}
		} else { /* stratified; will not occur if n_pred_locs == 0: */
			cp = get_outfile_namei(0);
			name = (char *) erealloc(name, (strlen(cp) + 12) * sizeof(char));
			map_name_nr(mask, cp, name, i , gl_nsim);
			what = (char *) erealloc(what, 40 * sizeof(char));
			sprintf(what, "stratified simulation %d", i);
			m = map_dup(name, mask);
			for (j = 0; j < get_n_vars(); j++) {
				restore_data_list(&d[j], i, n_pred_locs);
				for (k = 0; k < n_sim_locs[j]; k++) {
					if (map_xy2rowcol(m, 
							d[j]->list[s2d[j][k]]->x, 
							d[j]->list[s2d[j][k]]->y, 
							&r, &c)) {
						pr_warning("x %g y %g", 
							d[j]->list[s2d[j][k]]->x, 
							d[j]->list[s2d[j][k]]->y);
						ErrMsg(ER_IMPOSVAL, "map_xy2rowcol");
					}
					map_put_cell(m, r, c, d[j]->list[s2d[j][k]]->attr);
				}
			}
			map_sign(m, what);
			m->write(m);
			map_free(m);
		} /* else */
	} /* for i */
	efree(name);
	efree(what);
	free_simulations();
}
#endif

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
}

#ifndef USING_R
void lhs(DATA **d, int n_vars, int stratify) {
/*
 * lhs(): Latin hypercube sampling of multi-Gaussian random fields
 * written: Sept 1996
 * Dec 1997: added warning for small sample size
 *
 * given the table sim_values[][], create a Latin Hypercube Sample (LHS),
 * assuming that the marginal distribution of each variable is equal to
 * the unconditional distribution (N(mu(s_0), sill)) or the marginal distrib.
 * as defined by the `marginals' command.
 * 
 * references:
 * Stein, M.L. (1987) Large Sample Properties of Simulations
 *  Using Latin Hypercube Sampling. Technometrics 29 (2): 143-151.
 * Pebesma, E.J., and Heuvelink, G.B.M.,
 *  Latin Hypercube Sampling of Gaussian Random Fields
 *  (accepted for publication in Technometrics, planned for Nov 1999 issue)
 *
 * marginal distribution of variable k, 0 <= k <= K: 
 *   if !stratify: get from data[k % n_vars].
 *   else: get from data[strata[k]].
 */
	int i, j, k, map;
	unsigned int r, c;
	double mean = 0.0, variance = 0.0;
	Double_index *row_Y = NULL;
	VARIOGRAM *vp = NULL;
	GRIDMAP **lhsmean = NULL, **lhsvar = NULL;
	char *fname;

	if (DEBUG_DUMP || DEBUG_COV)
		printlog("\nLatin hypercube sampling...");
	if (gl_nsim <= 1) {
		if (DEBUG_DUMP || DEBUG_COV)
			printlog("\nno data to do so");
		return;
	}
	if (gl_nsim <= 15)
		pr_warning("small sample size (%d) may cause disturbed joint distributions", gl_nsim);
 	if (get_method() != GSI)
 		pr_warning("lhs(): will produce nonsense if not used for Gaussian simulation");
	if (gl_marginal_names) {
		lhsmean = (GRIDMAP **) emalloc((gl_n_marginals / 2) * sizeof(GRIDMAP *));
		lhsvar = (GRIDMAP **) emalloc((gl_n_marginals / 2) * sizeof(GRIDMAP *));
		for (i = 0; i < gl_n_marginals / 2; i++) {
			fname = string_dup(gl_marginal_names[i * 2]);
			lhsmean[i] = new_map(READ_ONLY);
			lhsmean[i]->filename = fname;
			if ((lhsmean[i] = map_read(lhsmean[i])) == NULL)
				ErrMsg(ER_READ, fname);
			if (DEBUG_COV)
				printlog("mean: %s, ", fname);
			fname = string_dup(gl_marginal_names[i * 2 + 1]);
			lhsvar[i] = new_map(READ_ONLY);
			lhsvar[i]->filename = fname;
			if ((lhsvar[i] = map_read(lhsvar[i])) == NULL)
				ErrMsg(ER_READ, fname);
			if (! map_equal(lhsmean[i], lhsvar[i]))
				ErrMsg(ER_IMPOSVAL, "lhsmean and lhsvar: map topology differs");
			if (DEBUG_COV)
				printlog("variance: %s\n", fname);
		}
	}
/*
 * initial check:
 */
 	row_Y = (Double_index *) emalloc(gl_nsim * sizeof(Double_index));
	for (i = 0; i < n_vars; i++) {
		if (d[i]->n_original > 0 && gl_marginal_names == NULL)
			pr_warning("lhs(): define conditional distributions maps");
		else if (gl_n_marginals == 0) { /* will derive marginal from sk_mean/vgm: */
			vp = get_vgm(LTI(i, i));
			if (! vp->is_valid_covariance)
				pr_warning("lhs(): variogram without sill will lead to invalid marginal distribution");
			if (d[i]->beta == NULL)
				ErrMsg(ER_VARNOTSET, "beta/sk_mean not defined for data");
			if (d[i]->beta->size > 1)
				pr_warning("marginal distributions wrt beta will be incorrect");
		}
	}
/*
 * for all variables do:
 */
 	for (i = 0; i < n_vars; i++) { /* loop over variables */
 		for (j = 0; j < n_sim_locs[i]; j++) { /* loop over sim. locations */

 			for (k = 0; k < gl_nsim; k++) {
 				row_Y[k].d = msim[i][j][k];
 				row_Y[k].index = k;
 			}
/*
 * get mean:
 */
			if (gl_marginal_names) { /* conditional mean & variance maps */
				if (stratify)
					map = 0;
				else
					map = i;
				if (map_xy2rowcol(lhsmean[map], 
						d[i]->list[s2d[i][j]]->x, d[i]->list[s2d[i][j]]->y, 
						&r, &c)) {
					pr_warning("x %g y %g", 
						d[i]->list[j]->x, d[i]->list[j]->y);
					ErrMsg(ER_IMPOSVAL, "map_xy2rowcol");
				} 
				mean = map_get_cell(lhsmean[map], r, c);
				variance = map_get_cell(lhsvar[map], r, c);
			} else if (gl_marginal_values) { 
				/* unconditional: prespecified marginal */
				mean = gl_marginal_values[i * 2];
				variance = gl_marginal_values[i * 2 + 1];
			} else { /* unconditional: get marginals from data */
 				/* mean = d[i]->sk_mean; */
 				mean = d[i]->beta->val[0]; /* indeed, inconsistent!! */
 				/* mean = calc_mu(d[i], d[i]->??? ...) */
 				vp = get_vgm(LTI(i,i));
 				variance = vp->sum_sills;
 			}

 			if (DEBUG_COV)
 				printlog("marginals for %s: mean %g, variance %g\n",
 					name_identifier(i), mean, variance);
 			if (variance < 0.0)
 				ErrMsg(ER_IMPOSVAL, "lhs(): negative conditional variance");

			lhs_one(row_Y, gl_nsim, mean, variance);
			for (k = 0; k < gl_nsim; k++)
				msim[i][j][k] = row_Y[k].d; /* copy back */

 		} /* for j */
 	} /* for i */

	if (gl_marginal_names) {
		for (i = 0; i < gl_n_marginals / 2; i++) {
			map_free(lhsmean[i]);
			map_free(lhsvar[i]);
		}
	}

	efree(row_Y);

	if (DEBUG_DUMP || DEBUG_COV)
		printlog("done");

	return;
}
#endif

static void lhs_one(Double_index *list, unsigned int dim, 
		double mean, double var) {
 	static int *rank = NULL, sizeof_rank = -1;
 	int i;
 	double xi;

 	if (rank == NULL || dim > sizeof_rank) {
 		rank = (int *) emalloc(dim * sizeof(int));
 		sizeof_rank = dim;
 	}

	qsort((void *) list, (size_t) dim, sizeof(Double_index),
		(int CDECL (*)(const void *, const void *)) double_index_cmp);
/* 
 * Get ranks -- after sorting row_Y on the ->d field,
 * row_Y[0].index,...,row_Y[gl_nsim-1].index contains the original
 * locations of consecutive elements in sim_values[j].
 * Now we get ranks by:
 */

	for (i = 0; i < dim; i++)
 		rank[list[i].index] = i;
/*
 * produce Z : use eq. (10) in Stein (1987)
 * (store in-place: replace .d values)
 */
	for (i = 0; i < dim; i++) {

 		if (gl_lhs == 2)
 			xi = 0.5; /* Owen's suggestion */
 		else
 			xi = r_uniform();

 		list[i].d = mean + sqrt(var) * q_normal( (rank[i] + xi) / dim);

 		if (DEBUG_COV)
 			printlog("%g [%d]%s", list[i].d, rank[i], (i + 1) % 5 ? " " : "\n");
 	}
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

#ifndef SIM_DOUBLE
float ***get_msim(void) {
	return msim;
}
#endif
