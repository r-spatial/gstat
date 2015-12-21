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
 * sem.c: calculate sample (cross, co-) variogram from data
 * K.M. refers to changes by Konstantin Malakhanov, see mapio.c
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <limits.h>

#include "config.h"

#ifdef USING_R
# include <R.h>
# include <Rdefines.h>
#endif

#include "defs.h"
#include "read.h"
#include "mapio.h"
#include "userio.h"
#include "data.h"
#include "utils.h"
#include "debug.h"
#include "vario.h"
#include "glvars.h"
#include "select.h"
#include "gls.h"
#include "lm.h"
#include "defaults.h"
#include "direct.h"
#include "version.h"
#include "sem.h"

#define SEM_INCREMENT 1000

static double valid_distance(DPOINT *a, DPOINT *b, double max, 
		int symmetric, DATA *d1, DATA *d2, GRIDMAP *m);
static void divide(SAMPLE_VGM *ev);
static SAMPLE_VGM *alloc_exp_variogram(DATA *a, DATA *b, SAMPLE_VGM *ev);
/* variograms: */
static SAMPLE_VGM *semivariogram(DATA *a, SAMPLE_VGM *ev);
static SAMPLE_VGM *cross_variogram(DATA *a, DATA *b, SAMPLE_VGM *ev);
/* covariograms: */
static SAMPLE_VGM *covariogram(DATA *a, SAMPLE_VGM *ev);
static SAMPLE_VGM *cross_covariogram(DATA *a, DATA *b, SAMPLE_VGM *ev);
static int get_index(double dist, SAMPLE_VGM *ev);
static void ev2map(VARIOGRAM *v);
static SAMPLE_VGM *load_ev(SAMPLE_VGM *ev, const char *fname);
static SAMPLE_VGM *semivariogram_list(DATA *d, SAMPLE_VGM *ev);
static SAMPLE_VGM *semivariogram_grid(DATA *d, SAMPLE_VGM *ev);
static void push_to_cloud(SAMPLE_VGM *ev, double gamma, double dist,
	unsigned long index);
static void resize_ev(SAMPLE_VGM *ev, unsigned int size);
static void *register_pairs(void *p, unsigned long nh, DPOINT *a, DPOINT *b);

/*
 * gl_cressie: use Cressie's sqrt(absdiff) estimator;
 * ev->zero: 
 *   case ZERO_INCLUDE: use zero distances in first interval (omit);
 *   case ZERO_AVOID:   avoid zero distances;
 *   case ZERO_SPECIAL: make special estimate for distance zero;
 */

/*
 * calculate sample variogram from data
 * calculates variogram, crossvariogram, covariogram or crosscovariogram
 * from sample data. Data are obtained from the central data base (glvars)
 * using get_gstat_data(), and the variogram requested is that of data id
 * v->id1 and v->id2 -- a direct (co) variogram when id1 == id2, a cross
 * (co) variogram when id1 != id2. 
 * 
 * if v->fname is set and (one of) id1 or id2 is a dummy data, the
 * actual sample variogram is not calculated but rather read from the
 * file v->fname. This is done to enable separate sample variogram
 * calculation (in batch or on a fast remote computer) and model fitting
 * (e.g. on the desk top).
 *
 * returns: non-zero if writing sample variogram to file failed.
 */
int calc_variogram(VARIOGRAM *v /* pointer to VARIOGRAM structure */,
		const char *fname /* pointer to output file name, or NULL if
		no output has to be written to file */ )
{
	DATA **d = NULL, *d1 = NULL, *d2 = NULL;
#ifndef USING_R
	FILE *f = NULL;
#endif

	assert(v);

	d = get_gstat_data();
	d1 = d[v->id1];
	d2 = d[v->id2];
#ifndef USING_R
	if (v->fname && (d1->dummy || d2->dummy)) {
		if ((v->ev = load_ev(v->ev, v->fname)) == NULL)
			ErrMsg(ER_READ, "could not read sample variogram");
		v->ev->cloud = 0;
		v->ev->recalc = 0;
		return 0;
	}
#endif
	if (d1->sel == NULL) 
		select_at(d1, NULL); /* global selection (sel = list) */
	if (d2->sel == NULL)
		select_at(d2, NULL);

	if (v->ev->evt == CROSSVARIOGRAM && 
			(v->ev->pseudo == -1 || v->ev->is_asym == -1)) {
		/* v's first time */
		if (coordinates_are_equal(d[v->id1], d[v->id2]))
			v->ev->pseudo = 0;
		else
			v->ev->pseudo = 1;
		if (gl_sym_ev == 0)
			v->ev->is_asym = v->ev->pseudo;
			/* pseudo: always, else: only if set */
		else
			v->ev->is_asym = 0;
	}
	if (gl_zero_est == ZERO_DEFAULT) { /* choose a suitable default */
		if (is_covariogram(v))
			v->ev->zero = ZERO_SPECIAL;
		else { /* v is variogram */
			if (v->ev->pseudo)
				v->ev->zero = ZERO_SPECIAL;
			else
				v->ev->zero = ZERO_INCLUDE;
		}
	} else
		v->ev->zero = zero_int2enum(gl_zero_est);

	assert(v->ev->zero != ZERO_DEFAULT);

	fill_cutoff_width(d1, v);

	if (v->ev->map && v->fname == NULL && v->ev->S_grid == NULL)
		return -1;

	v->ev->cloud = (v->ev->iwidth <= 0.0);
	if (v->ev->cloud &&
			 (d[v->id1]->n_sel >= MAX_NH ||  d[v->id2]->n_sel >= MAX_NH))
		pr_warning("observation numbers in cloud will be wrong");
	set_direction_values(gl_alpha, gl_beta, gl_tol_hor, gl_tol_ver);

	v->ev->is_directional = is_directional(v);
	if (v->ev->recalc) {
		switch (v->ev->evt) {
			case PRSEMIVARIOGRAM:
			case SEMIVARIOGRAM:
				semivariogram(d[v->id1], v->ev);
				break;
			case CROSSVARIOGRAM:
				cross_variogram(d[v->id1], d[v->id2], v->ev);
				break;
			case COVARIOGRAM:
				v->ev->is_asym = gl_sym_ev;
				covariogram(d[v->id1], v->ev);
				break;
			case CROSSCOVARIOGRAM:
				cross_covariogram(d[v->id1], d[v->id2], v->ev);
				break;
			case NOTSPECIFIED:
			default:
				assert(0); /* aborts */
				break;
		}
	}
#ifndef USING_R
	if (v->ev->map && !v->ev->S_grid)
		ev2map(v);
	else if (fname != NULL) {
		f = efopen(fname, "w");
		fprint_header_vgm(f, d[v->id1], d[v->id2], v->ev);
		fprint_sample_vgm(f, v->ev);
	}
	if (f && f != stdout)
		return efclose(f);
#endif
	return 0;
}

static SAMPLE_VGM *semivariogram(DATA *d, SAMPLE_VGM *ev) {
/*
 *  calculate sample variogram of 0.5 E[(Z(x)-Z(x+h))2]
 */
	if (ev->evt == PRSEMIVARIOGRAM)
		d->calc_residuals = 0;
	ev = alloc_exp_variogram(d, NULL, ev);
	if (d->grid != NULL && d->prob > 0.5 && d->every == 1)
		ev = semivariogram_grid(d, ev);
	else 
		ev = semivariogram_list(d, ev);
	divide(ev);
	ev->recalc = 0;
	return ev;
} /* semivariogram() */

static SAMPLE_VGM *semivariogram_list(DATA *d, SAMPLE_VGM *ev) {
	unsigned long uli, ulj;
	int i, j, index = 0, divide_by = 1;
	unsigned int total_steps;
	double gamma, ddist, head, tail, gam;

	while (d->n_sel / divide_by > 0.5 * sqrt(INT_MAX))
		divide_by <<= 1; /* prevent overflow on calculating total_steps */

	total_steps = (d->n_sel / divide_by) * (d->n_sel - 1) / 2;
	print_progress(0, total_steps);
	
	if (DEBUG_DUMP)
		printlog("Calculating semivariogram from %d points...\n", d->n_sel);

	for (i = 0; i < d->n_sel; i++) {
		print_progress((i / divide_by) * (i - 1) / 2, total_steps);
#ifdef USING_R
		R_CheckUserInterrupt();
#endif
		/*
		printlog("step: %u of %u\n", (i /divide_by) * (i - 1) / 2, total_steps);
		*/
		for (j = 0; j < (ev->map != NULL ? d->n_sel : i); j++) {
			ddist = valid_distance(d->sel[i], d->sel[j], ev->cutoff, 1,
				d, d, (GRIDMAP *) ev->map);
			if (ddist >= 0.0 && i != j) {
				head = d->sel[i]->attr;
				tail = d->sel[j]->attr;
				if (! ev->cloud) {
					index = get_index(ddist, ev);
					if (gl_cressie)  /* sqrt abs diff */
						ev->gamma[index] += sqrt(fabs(head - tail));
					else if (ev->evt == PRSEMIVARIOGRAM) {
						gam = 2.0 * (head - tail)/(head + tail);
						ev->gamma[index] += SQR(gam);
					} else { /* SEMIVARIOGRAM: */
						ev->gamma[index] += SQR(head - tail);
#ifdef ADJUST_VARIANCE
						if (d->colnvariance)
							ev->gamma[index] -= d->sel[i]->variance +
									 d->sel[j]->variance;
#endif
					}
					ev->dist[index] += ddist;
					ev->pairs[index] = register_pairs(ev->pairs[index],
						ev->nh[index], d->sel[i], d->sel[j]);
					ev->nh[index]++;
				} else { /* cloud: */
					if (! (ev->zero == ZERO_AVOID && ddist == 0.0)) {
						if (gl_cressie)
							gamma = sqrt(fabs(head - tail));
						else if (ev->evt == PRSEMIVARIOGRAM) {
							gam = 2.0 * (head - tail)/(head + tail);
							gamma = gam * gam;
						} else {
							gamma = SQR(head - tail);
#ifdef ADJUST_VARIANCE
							if (d->colnvariance)
								gamma -= d->sel[i]->variance + d->sel[j]->variance;
#endif
						}
						uli = i; ulj = j;
						push_to_cloud(ev, gamma / 2.0, ddist, TO_NH(uli,ulj));
					}
				}
			}/* if ddist >= 0 */
		}  /* for j */
	} /* for i */
	print_progress(total_steps, total_steps);
	if (DEBUG_DUMP)
		printlog("ready\n");
	return ev;
}

static SAMPLE_VGM *semivariogram_grid(DATA *d, SAMPLE_VGM *ev) {

	typedef struct {
		int row, col, ev_index;
		double dist;
	} grid_index;
	struct {
		int n;
		grid_index *gi;
	} grid_ev;
	int row, col, irow, icol, i, max_index, index;
	unsigned long ula, ulb;
	double gamma, ddist, head, tail, gam;
	DPOINT a, b, *dpa = NULL, *dpb = NULL;

	max_index = (int) floor(ev->cutoff / SQUARECELLSIZE(d->grid));
	grid_ev.gi = (grid_index *) emalloc(2 * (max_index + 1) * (max_index + 1)
		* sizeof(grid_index));
	grid_ev.n = 0;
	a.x = a.y = a.z = b.z = 0.0;
	/* setup the grid: */
	for (row = 0; row <= max_index; row++) {
		for (col = (row == 0 ? 1 : -max_index); col <= max_index; col++) {
			b.x = col * SQUARECELLSIZE(d->grid);
			b.y = row * SQUARECELLSIZE(d->grid);
			ddist = valid_distance(&a, &b, ev->cutoff, 1,
				d, d, (GRIDMAP *) ev->map);
			if (ddist > 0.0) {
				grid_ev.gi[grid_ev.n].row = row;
				grid_ev.gi[grid_ev.n].col = col;
				grid_ev.gi[grid_ev.n].dist = ddist;
				if (! ev->cloud)
					grid_ev.gi[grid_ev.n].ev_index = get_index(ddist, ev);
				if (DEBUG_DUMP)
					printlog("row %d col %d index %d\n", 
						row, col, grid_ev.gi[grid_ev.n].ev_index);
				grid_ev.n++;
			}
		}
	}
	print_progress(0, d->grid->rows);
	for (row = 0; row < d->grid->rows; row++) {
		for (col = 0; col < d->grid->cols; col++) {
			if ((dpa = d->grid->dpt[row][col]) != NULL) {
				for (i = 0; i < grid_ev.n; i++) {
					irow = row + grid_ev.gi[i].row;
					icol = col + grid_ev.gi[i].col;
					if (irow >= 0 && icol >= 0 && irow < d->grid->rows
							&& icol < d->grid->cols
							&& ((dpb = d->grid->dpt[irow][icol]) != NULL)) {
						ddist = grid_ev.gi[i].dist;
						head = dpa->attr;
						tail = dpb->attr;
						if (! ev->cloud) {
							index = grid_ev.gi[i].ev_index;
							if (gl_cressie)  /* sqrt abs diff */
								ev->gamma[index] += sqrt(fabs(head - tail));
							else {
								if (ev->evt == PRSEMIVARIOGRAM) {
									gam = 2.0 * (head - tail)/(head + tail);
									ev->gamma[index] += gam * gam;
								} else
									ev->gamma[index] += SQR(head - tail);
#ifdef ADJUST_VARIANCE
								if (d->colnvariance)
									ev->gamma[index] -= dpa->variance +
										 	dpb->variance;
#endif
							}
							ev->dist[index] += ddist;
							ev->pairs[index] = register_pairs(ev->pairs[index],
								ev->nh[index], dpa, dpb);
							ev->nh[index]++;
						} else { /* cloud: */
							if (gl_cressie)
								gamma = sqrt(fabs(head - tail));
							else if (ev->evt == PRSEMIVARIOGRAM) {
								gam = 2.0 * (head - tail)/(head + tail);
								gamma = gam * gam;
							} else {
								gamma = SQR(head - tail);
#ifdef ADJUST_VARIANCE
								if (d->colnvariance)
									gamma -= dpa->variance + dpb->variance;
#endif
							}
							ula = GET_INDEX(dpa);
							ulb = GET_INDEX(dpb);
							push_to_cloud(ev, gamma / 2.0, ddist, TO_NH(ula,ulb));
						} /* else !cloud */
					} /* if we have two non-NULL points */
				} /* for all possibly relevant pairs */
			} /* if this grid cell is non-NULL */
		} /* for all cols */
		print_progress(row + 1, d->grid->rows);
#ifdef USING_R
		R_CheckUserInterrupt();
#endif
	} /* for all rows */
	efree(grid_ev.gi);
	return ev;
}

/* covariograms: */
static SAMPLE_VGM *covariogram(DATA *d, SAMPLE_VGM *ev) {
	int i, j, index = 0;
	unsigned long uli, ulj;
	double gamma, ddist;

	ev->evt = COVARIOGRAM;
	ev = alloc_exp_variogram(d, NULL, ev);
	for (i = 0; i < d->n_sel; i++) {
		print_progress(i, d->n_sel);
#ifdef USING_R
		R_CheckUserInterrupt();
#endif
		for (j = 0; j <= (ev->map != NULL ? d->n_sel-1 : i); j++) {
			ddist = valid_distance(d->sel[i], d->sel[j], ev->cutoff, 1,
				d, d, (GRIDMAP *) ev->map);
			if (ddist >= 0.0) {
				if (! ev->cloud) {
					index = get_index(ddist, ev);
					ev->gamma[index] += d->sel[i]->attr * d->sel[j]->attr;
#ifdef ADJUST_VARIANCE
					if (d->colnvariance && i == j)
						ev->gamma[index] -= d->sel[i]->variance;
#endif
					ev->dist[index] += ddist;
					ev->pairs[index] = register_pairs(ev->pairs[index],
						ev->nh[index], d->sel[i], d->sel[j]);
					ev->nh[index]++;
				} else {
					if (! (ev->zero == ZERO_AVOID && ddist == 0.0)) {
						gamma = d->sel[i]->attr * d->sel[j]->attr;
#ifdef ADJUST_VARIANCE
						if (d->colnvariance && i == j)
							gamma -= d->sel[i]->variance;
#endif
						uli = i;
						ulj = j;
						push_to_cloud(ev, gamma, ddist, TO_NH(uli,ulj));
					}
				}
			}/* if ddist >= 0 */
		}  /* for j */
	} /* for i */
	print_progress(d->n_sel, d->n_sel);
	divide(ev);
	ev->recalc = 0;
	return ev;
} /* covariogram() */

static SAMPLE_VGM *cross_variogram(DATA *a, DATA *b, SAMPLE_VGM *ev) {
	int i, j, index = 0;
	unsigned long uli, ulj;
	double gamma, ddist;

	ev->evt = CROSSVARIOGRAM;
	ev = alloc_exp_variogram(a, b, ev);
	for (i = 0; i < a->n_sel; i++) {
		print_progress(i, a->n_sel);
#ifdef USING_R
		R_CheckUserInterrupt();
#endif
		for (j = 0; j < b->n_sel; j++) {
			ddist = valid_distance(a->sel[i], b->sel[j], ev->cutoff,
				gl_sym_ev || !ev->pseudo, a, b, (GRIDMAP *) ev->map); 
			if (ddist >= 0.0) {
				if (!ev->pseudo && i != j) {
					if (! ev->cloud) {
						index = get_index(ddist, ev);
						ev->gamma[index] += 
							(a->sel[i]->attr - a->sel[j]->attr) *
							(b->sel[i]->attr - b->sel[j]->attr);
						ev->dist[index] += ddist;
						ev->pairs[index] = register_pairs(ev->pairs[index],
							ev->nh[index], a->sel[i], a->sel[j]);
						ev->nh[index]++;
					} else if (!(ddist == 0.0 && ev->zero == ZERO_AVOID)) { 
						gamma = (a->sel[i]->attr - a->sel[j]->attr) *
							(b->sel[i]->attr - b->sel[j]->attr);
						uli = i; 
						ulj = j;
						push_to_cloud(ev, gamma / 2.0, ddist, TO_NH(uli,ulj));
					}
				} else if (ev->pseudo) {
					if (! ev->cloud) {
						index = get_index(ddist, ev);
						ev->gamma[index] += 
							SQR(a->sel[i]->attr - b->sel[j]->attr);
#ifdef ADJUST_VARIANCE
						if (a->colnvariance || b->colnvariance)
							ev->gamma[index] -= a->sel[i]->variance +
								 b->sel[j]->variance;
#endif
						ev->dist[index] += ddist;
						ev->pairs[index] = register_pairs(ev->pairs[index],
							ev->nh[index], a->sel[i], b->sel[j]);
						ev->nh[index]++;
					} else if (! (ev->zero == ZERO_AVOID && ddist == 0.0)) {
						gamma = SQR(a->sel[i]->attr - b->sel[j]->attr);
#ifdef ADJUST_VARIANCE
						if (a->colnvariance || b->colnvariance)
							gamma -= a->sel[i]->variance + b->sel[j]->variance;
#endif
						uli = i;
						ulj = j;
						push_to_cloud(ev, gamma / 2.0, ddist, TO_NH(uli,ulj));
					}
				}
			}/* if ddist >= 0 */
		}  /* for j */
	} /* for i */
	print_progress(a->n_sel, a->n_sel);
	divide(ev);
	ev->recalc = 0;
	return ev;
} /* cross_variogram */

static SAMPLE_VGM *cross_covariogram(DATA *a, DATA *b, SAMPLE_VGM *ev) {
	int i, j, index = 0;
	unsigned long uli, ulj;
	double gamma, ddist;

	ev->evt = CROSSCOVARIOGRAM;
	ev = alloc_exp_variogram(a, b, ev);
	for (i = 0; i < a->n_sel; i++) {      /* i -> a */
#ifdef USING_R
		R_CheckUserInterrupt();
#endif
		print_progress(i, a->n_sel);
		for (j = 0; j < b->n_sel; j++) {  /* j -> b */
			ddist = valid_distance(a->sel[i], b->sel[j], ev->cutoff,
				gl_sym_ev, a, b, (GRIDMAP *) ev->map); 
			if (ddist >= 0.0) {
				if (! ev->cloud) {
					index = get_index(ddist, ev);
					ev->gamma[index] += a->sel[i]->attr * b->sel[j]->attr;
					ev->dist[index] += ddist;
					ev->pairs[index] = register_pairs(ev->pairs[index],
						ev->nh[index], a->sel[i], b->sel[j]);
					ev->nh[index]++;
				} else if (! (ev->zero == ZERO_AVOID && ddist == 0.0)) {
					gamma = a->sel[i]->attr * b->sel[j]->attr;
					uli = i; 
					ulj = j;
					push_to_cloud(ev, gamma, ddist, TO_NH(uli,ulj));
				}
			}/* if ddist >= 0 */
		}  /* for j */
	} /* for i */
	print_progress(a->n_sel, a->n_sel);
	divide(ev);
	ev->recalc = 0;
	return ev;
} /* cross_covariogram() */

static double valid_distance(DPOINT *a, DPOINT *b, double max, 
		int symmetric, DATA *d1, DATA *d2, GRIDMAP *map) {
	double ddist, dX, dX2, inprod;
	DPOINT p;
	int /* mode = 0, */ i;
	unsigned int row, col;

	assert(a != NULL);
	assert(b != NULL);
	assert(d1 != NULL);
	assert(d2 != NULL);

	/* mode = d1->mode & d2->mode; */
/*
 * even if modes don't correspond, valid_direction() will 
 * calculate valid distances 
 */
	p.x = a->x - b->x;
	p.y = a->y - b->y;
	p.z = a->z - b->z;
	if (map && !gl_longlat) {
		/* transform here p to allow directional 2d cuts in a 3d world */
		if (map_xy2rowcol(map, p.x, p.y, &row, &col))
			return -1.0;
		else
			ddist = (1.0 * row) * map->cols + col + 0.5;
	} else {
		if (!gl_longlat && (p.x > max || p.y > max || p.z > max))
			return -1.0;
		/* Changed K.M. Fri Feb 27 15:56:57 1998 */
		/* if ddist < 0.0 then we don't need to check for dX! */
		if ((ddist = valid_direction(a, b, symmetric, d1)) > max || ddist < 0.0) 
			return -1.0;
	}
	dX = MIN(d1->dX, d2->dX);
	if (dX < DBL_MAX) {
		dX2 = dX * dX;
		/* allow only points for which 2-norm ||x_i-x_j|| < dX */
		if (d1->n_X != d2->n_X)
			ErrMsg(ER_IMPOSVAL, "valid_distance(): d1->n_X != d2->n_X");
		for (i = 0, inprod = 0.0; i < d1->n_X; i++) {
			inprod += SQR(a->X[i] - b->X[i]);
			/* printf("a->X[%d]: %g, b->X[%d]: %g", i, a->X[i], i, b->X[i]); */
		}
		if (inprod > dX2)
			ddist = -1.0;
		/* printf("dX2: %g, inprod: %g ddist: %g\n", dX2, inprod, ddist); */
	}
	/*
	if (d1->coln_id > 0 && d2->coln_id > 0 && strcmp())
		return -1.0;
	*/
	return ddist;
}

int is_directional(VARIOGRAM *v) {
	switch(v->ev->evt) {
		case CROSSCOVARIOGRAM:
			if (gl_sym_ev == 0) /* asymmetric cross(co)variances: */
				return (gl_tol_hor < 180.0 || gl_tol_ver < 180.0);
			else
				return (gl_tol_hor < 90.0 || gl_tol_ver < 90.0);
		case CROSSVARIOGRAM:
			if (v->ev->is_asym && gl_sym_ev == 0) /* asymm. cross(co)variances: */
				return (gl_tol_hor < 180.0 || gl_tol_ver < 180.0);
			else
				return (gl_tol_hor < 90.0 || gl_tol_ver < 90.0);
		default: /* symmetric (co)variances */
			return (gl_tol_hor < 90.0 || gl_tol_ver < 90.0);
	}
}

/*
 * this function should be changed--the mask map stack is misused as
 * to define the topology of variogram maps.
 *
 * use min/max coordinates for block diagonal as maximum cutoff
 * Returns: about 1/3 the max. dist between any two points in data. 
 */
void fill_cutoff_width(DATA *data /* pointer to DATA structure to derive
		the values from */,
		VARIOGRAM *v /* pointer to VARIOGRAM structure */)
{
	double d = 0.0;
	int i;
	GRIDMAP *m;
	DATA_GRIDMAP *dg;
	SAMPLE_VGM *ev;

	assert(data);
	assert(v);

	ev = v->ev;
	if ((get_n_masks() > 0 && get_method() != LSEM) || ev->S_grid != NULL) {
		m = new_map(READ_ONLY);
		if (ev->S_grid) {
			/* process S_grid to m */
			dg = (DATA_GRIDMAP *) ev->S_grid;
			m->x_ul = dg->x_ul;
			m->y_ul = dg->y_ul;
			m->cellsizex = dg->cellsizex;
			m->cellsizey = dg->cellsizey;
			m->rows = dg->rows;
			m->cols = dg->cols;
		} else {
#ifndef USING_R
			m->filename = get_mask_name(0);
			if ((m = map_read(m)) == NULL)
				ErrMsg(ER_READ, "cannot open map");
#endif
		}
		ev->iwidth = 1.0;
		ev->cutoff = m->rows * m->cols; 
			/* not a real cutoff, but rather the size of the container array */
		ev->map = m;
	} else if (gl_bounds != NULL) {
		i = 0;
		while (gl_bounds[i] >= 0.0) /* count length */
			i++;
		ev->cutoff = gl_bounds[i-1];
		ev->iwidth = ev->cutoff / i;
	} else {
		if (is_mv_double(&(ev->cutoff))) {
			if (gl_cutoff < 0.0) {
				d = data_block_diagonal(data);
				if (d == 0.0)
					ev->cutoff = 1.0; /* ha ha ha */
				else
					ev->cutoff = d * gl_fraction;
			} else
				ev->cutoff = gl_cutoff;
		}
		if (is_mv_double(&(ev->iwidth))) {
			if (gl_iwidth < 0.0)
				ev->iwidth = ev->cutoff / gl_n_intervals;
			else
				ev->iwidth = gl_iwidth;
		}
	}
}

#ifndef USING_R
void fprint_header_vgm(FILE *f, const DATA *a, const DATA *b,
			const SAMPLE_VGM *ev) {
	time_t tp;
	char *cp = NULL;
	/* char *pwd; */

	fprintf(f, "#gstat %s %s [%s]", GSTAT_OS, VERSION, command_line);
	/* if (pwd = getenv("PWD")) fprintf(f, "(in %s)", pwd); */
	fprintf(f, "\n");
	fprintf(f, "#sample %s%s\n",
		(ev->evt == CROSSVARIOGRAM && ev->pseudo ? "pseudo " : ""),
		vgm_type_str[ev->evt]);
	tp = time(&tp);
	fprintf(f, "#%s", asctime(localtime(&tp))); /* includes \n */
	cp = print_data_line(a, &cp);
	fprintf(f, "#data(%s): %s", name_identifier(a->id), cp);
	if (a != b) {
		cp = print_data_line(b, &cp);
		fprintf(f, " data(%s): %s", name_identifier(b->id), cp);
	}
	if (cp != NULL)
		efree(cp);
	fprintf(f, "\n");
	fprintf(f, "#[1] mean: %g variance: %g", a->mean, a->std * a->std);
	if (a != b)
		fprintf(f, " [2] mean: %g variance: %g", b->mean, b->std * a->std);
	fprintf(f, "\n");
	fprintf(f, "#cutoff: %g ", ev->cutoff);
	if (gl_bounds == NULL)
		fprintf(f,"%s %g\n","interval width:", ev->iwidth);
	else
		fprintf(f, "(fixed boundaries)\n");
	if (! ev->is_directional)
		fprintf(f, "#direction: total ");
	else
		fprintf(f,"#alpha <x,y>: %gd +/- %g; beta <alpha,z>: %gd +/- %g%s", 
			gl_alpha, gl_tol_hor, gl_beta, gl_tol_ver,
			gl_sym_ev ? " (symmetric)" : "");
	fprintf(f, "\n");
	if (ev->cloud) 
		fprintf(f, "# i  j distance %s cloud\n", vgm_type_str[ev->evt]);
	else
		fprintf(f, "#   from       to  n_pairs  av_dist %s\n",
			vgm_type_str[ev->evt]);
	return;
}
#endif

static SAMPLE_VGM *alloc_exp_variogram(DATA *a, DATA *b, SAMPLE_VGM *ev) {
	int i;
	double nd;

	assert(a != NULL);
	assert(ev != NULL);

	if (gl_zero_est != ZERO_DEFAULT && ev->zero != gl_zero_est)
		ev->zero = zero_int2enum(gl_zero_est);

	if (gl_gls_residuals) {
		if (a->calc_residuals)
			make_gls(a, 1);
		if (b != NULL && b->calc_residuals)
			make_gls(b, 1);
	} else {
		if (a->calc_residuals)
			make_residuals_lm(a);
		if (b != NULL && b->calc_residuals)
			make_residuals_lm(b);
	}

	if (ev->cloud) {
		ev->n_est = 0;
		return ev;
	}

	if (gl_bounds != NULL) {
		for (i = ev->n_est = 0; gl_bounds[i] >= 0.0; i++)
			ev->n_est++;
	} else {
/* check for overflow: */
		nd = floor(ev->cutoff / ev->iwidth) + 1;
		if (nd > INT_MAX) {
			pr_warning("choose a larger width or a smaller cutoff value");
			ErrMsg(ER_MEMORY, "(experimental variogram too large)");
		}
		ev->n_est = (int) nd;
	}
/*
 * zero est go to ev->gamma[ev->n_est - 1], ev->nh[ev->n_est - 1]
 */
	if (ev->zero)
		ev->n_est++;

	resize_ev(ev, ev->n_est);

/* initialize: */
	for (i = 0; i < ev->n_est; i++) {
		ev->gamma[i] = 0.0;
		ev->dist[i] = 0.0;
		ev->nh[i] = 0; 
		ev->pairs[i] = (DPOINT **) NULL;
	}
	return ev;
}

static void resize_ev(SAMPLE_VGM *ev, unsigned int size) {
	if (size > ev->n_max) {
		ev->n_max = size;
		ev->gamma = (double *) erealloc (ev->gamma, ev->n_max * sizeof(double));
		ev->dist = (double *) erealloc (ev->dist, ev->n_max * sizeof(double));
		ev->nh = (unsigned long *) erealloc (ev->nh, ev->n_max * sizeof(long));
		ev->pairs = (DPOINT ***) 
				erealloc(ev->pairs, ev->n_max * sizeof(DPOINT **));
	}
}

static void *register_pairs(void *pairs, unsigned long nh,
		DPOINT *a, DPOINT *b) {
/* 
 * while I'm here -- there may be a problem when ->list != ->sel on
 * the DATA used, but I don't know why. Probably will never be used.
 */
	/* resize pairs; add a and b to it */
	if (gl_register_pairs == 0)
		return NULL;
	if (nh % SEM_INCREMENT == 0)
		pairs = erealloc(pairs, 2 * (nh + SEM_INCREMENT + 1) * sizeof(DPOINT **));
	((DPOINT **) pairs)[2 * nh] = a;
	((DPOINT **) pairs)[2 * nh + 1] = b;
	return pairs;
}

static void push_to_cloud(SAMPLE_VGM *ev, double gamma, double dist,
	unsigned long index) {

	if (ev->n_est == ev->n_max) 
		resize_ev(ev, ev->n_max + SEM_INCREMENT);

	ev->gamma[ev->n_est] = gamma;
	ev->dist[ev->n_est]  = dist;
	ev->nh[ev->n_est]    = index;
	ev->pairs[ev->n_est] = NULL;
	ev->n_est++;
}

static int get_index(double dist, SAMPLE_VGM *ev) {
	double frac;
	int i = 0;

	if (dist == 0.0 && ev->zero != ZERO_INCLUDE)
		return ev->n_est - 1;

	if (gl_bounds != DEF_bounds) {
		for (i = 0; gl_bounds[i] >= 0.0; i++)
			if (dist <= gl_bounds[i])
				return i;
		assert(0);
	}

	if (ev->iwidth <= 0.0) {
		pr_warning("iwidth: %g", ev->iwidth);
		ErrMsg(ER_IMPOSVAL, "ev->iwidth <= 0.0");
	}

	frac = dist / ev->iwidth;
	if (dist > 0.0 && frac == floor(frac))
		return (int) (floor(frac)) - 1; 
	else
		return (int) floor(frac);
}

static void divide(SAMPLE_VGM *ev) {
	int i;

	if (ev->cloud)
		return; /* has been done in the first round */
	for (i = 0; i < ev->n_est; i++) {
		if (ev->nh[i]) {
			ev->dist[i] /= ev->nh[i];
			switch (ev->evt) {
				case SEMIVARIOGRAM:
					if (gl_cressie)
						ev->gamma[i] = 0.5 * pow(ev->gamma[i]/ev->nh[i], 4.0)
							/(0.457 + 0.494 / ev->nh[i]);
					else
						ev->gamma[i] /= (2.0 * ev->nh[i]);
					break;
				case CROSSVARIOGRAM:
					ev->gamma[i] /= (2.0 * ev->nh[i]);
					break;
				case COVARIOGRAM: /* BREAKTHROUGH */
				case CROSSCOVARIOGRAM:
					ev->gamma[i] /= (1.0 * ev->nh[i]);
					break;
				case PRSEMIVARIOGRAM:
					ev->gamma[i] /= (2.0 * ev->nh[i]);
					break;
				case NOTSPECIFIED: /* BREAKTHROUGH */
				default:
					assert(0);
					break;
			}
		}
	}
}

void fprint_sample_vgm(FILE *f, const SAMPLE_VGM *ev) {
#define EVFMT "%8g %8g %8lu %8g %8g\n"
	int i, n;
	double from, to;

	if (! ev->cloud) {
		/* start writing: distance 0 */
		if (ev->zero == ZERO_SPECIAL && ev->nh[ev->n_est-1])
#ifndef USING_R
			Rprintf(EVFMT, 0.0, 0.0, ev->nh[ev->n_est-1], 
				ev->dist[ev->n_est-1], ev->gamma[ev->n_est-1]);
#else 
			fprintf(f, EVFMT, 0.0, 0.0, ev->nh[ev->n_est-1], 
				ev->dist[ev->n_est-1], ev->gamma[ev->n_est-1]);
#endif
		/* continue writing: */
		if (ev->zero == ZERO_SPECIAL || ev->zero == ZERO_AVOID)
			n = ev->n_est - 1;
		else
			n = ev->n_est;
		for (i = 0; i < n; i++) {
			if (ev->nh[i] > 0) {
				if (gl_bounds == NULL) {
					from = i*ev->iwidth;
					to = (i+1)*ev->iwidth;
				} else {
					if (i == 0)
						from = 0.0;
					else
						from = gl_bounds[i-1];
					to = gl_bounds[i];
				}
				to = MIN(ev->cutoff, to);
#ifndef USING_R
				Rprintf(EVFMT, from, to, ev->nh[i],
					ev->dist[i], ev->gamma[i]);
#else
				fprintf(f, EVFMT, from, to, ev->nh[i],
					ev->dist[i], ev->gamma[i]);
#endif
				/*
				for (j = 0; j < ev->nh[i]; j++)
					fprintf(f, "[%d,%d] ",
						GET_INDEX(((DPOINT ***)ev->pairs)[i][2*j]),
						GET_INDEX(((DPOINT ***)ev->pairs)[i][2*j+1]));
				fprintf(f, "\n");
				*/
			}
		}
	} else {
		for (i = 0; i < ev->n_est; i++)
#ifndef USING_R
			Rprintf("%ld %ld %g %g\n", HIGH_NH(ev->nh[i]) + 1,
				LOW_NH(ev->nh[i]) + 1, ev->dist[i], ev->gamma[i]);
#else
			fprintf(f, "%ld %ld %g %g\n", HIGH_NH(ev->nh[i]) + 1,
				LOW_NH(ev->nh[i]) + 1, ev->dist[i], ev->gamma[i]);
#endif
	}
#ifndef USING_R
	fflush(f);
#endif
	return;
} /* fprint_sample_vgm */

static void ev2map(VARIOGRAM *v) {
	GRIDMAP *m1 = NULL, *m2 = NULL;
	unsigned int row, col, i;
	SAMPLE_VGM *ev;

	if (v->fname == NULL)
		return;
	ev = v->ev;
	m1 = map_dup(v->fname, ev->map);
	if (v->fname2 != NULL)
		m2 = map_dup(v->fname2, ev->map);
	for (row = i = 0; row < m1->rows; row++) {
		for (col = 0; col < m1->cols; col++) {
			if (ev->nh[i] > 0)
				map_put_cell(m1, row, col, ev->gamma[i]);
			if (m2 != NULL)
				map_put_cell(m2, row, col, 1.0 * ev->nh[i]);
			i++;
		}
	}
	m1->write(m1);
	if (m2 != NULL)
		m2->write(m2);
	return;
}

#ifndef USING_R
static SAMPLE_VGM *load_ev(SAMPLE_VGM *ev, const char *fname) {
	char *s = NULL, *tok;
	int i, size = 0, incr = 100;
	unsigned long l;
	FILE *f;

	f = efopen(fname, "r");
	if (ev == NULL)
		ev = init_ev();
	ev->evt = SEMIVARIOGRAM;
	for (i = 1; i <= 8; i++) {
		get_line(&s, &size, f);
		if (i == 6) {
			tok = strtok(s, " "); /* word */
			tok = strtok(NULL, " "); /* cutoff */
			if (read_double(tok, &(ev->cutoff))) {
				fclose(f); efree(s);
				pr_warning("file: %s, line: %d, token: %s", fname, i, tok);
				return NULL;
			}
			tok = strtok(NULL, " "); /* word */
			tok = strtok(NULL, " "); /* word */
			tok = strtok(NULL, " \n"); /* iwidth */
			if (tok != NULL) {
				if (read_double(tok, &(ev->iwidth))) {
					fclose(f); efree(s);
					pr_warning("file: %s, line: %d, token: %s", fname, i, tok);
					return NULL;
				}
			} /* else part: what to do with ev->iwidth? */
		}
	}
	while (get_line(&s, &size, f) != NULL) {
		ev->n_est++;
		if (ev->n_est >= ev->n_max) {
			ev->n_max += incr;
			ev->gamma = (double *) erealloc
				(ev->gamma, sizeof(double) * ev->n_max);
			ev->dist = (double *) erealloc
				(ev->dist, sizeof(double) * ev->n_max);
			ev->nh = (unsigned long *) erealloc
				(ev->nh, sizeof(long) * ev->n_max);
		}
		tok = strtok(s, " "); /* from */
		tok = strtok(NULL, " "); /* to */
		tok = strtok(NULL, " "); /* nh */
		if (read_ulong(tok, &l)) {
			fclose(f); efree(s);
			pr_warning("file: %s, line: %d, token: %s", fname, ev->n_est+8, tok);
			return NULL;
		}
		ev->nh[ev->n_est-1] = l;
		tok = strtok(NULL, " "); /* dist */
		if (read_double(tok, &(ev->dist[ev->n_est-1]))) {
			fclose(f); efree(s);
			pr_warning("file: %s, line: %d, token: %s", fname, ev->n_est+8, tok);
			return NULL;
		}
		tok = strtok(NULL, " \n"); /* semivariance or whatever */
		if (read_double(tok, &(ev->gamma[ev->n_est-1]))) {
			fclose(f); efree(s);
			pr_warning("file: %s, line: %d, token: %s", fname, ev->n_est+8, tok);
			return NULL;
		}
	}
	efree(s);
	efclose(f);
	return ev;
}
#endif
