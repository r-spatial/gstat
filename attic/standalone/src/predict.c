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
 Module for prediction or simulation:
 - loop over all prediction/simulation locations (regular or random path);
 - make selection at that location;
 - make an estimate (predict/simulate) at that location
   (if mask map has MV at the location, write MV)
 - write estimate to file
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "defs.h"
#include "mapio.h"
#include "userio.h"
#include "data.h"
#include "utils.h"
#include "debug.h"
#include "block.h"
#include "report.h"
#include "glvars.h"
#include "random.h"
#include "version.h"
#include "getest.h"
#include "msim.h"
#include "sim.h"
#include "select.h"
#include "predict.h"

typedef enum {
	AT_POINTS,
	AT_GRIDMAP
} PRED_AT;

static void init_predictions(PRED_AT w);
static DPOINT *next_location(DPOINT *loc, PRED_AT what, int random_path,
		unsigned int *row, unsigned int *col, DATA **data);
static void exit_predictions(PRED_AT w);
static void write_output(double *est, PRED_AT w, DPOINT *p, unsigned int row, unsigned int col);
static GRIDMAP *check_open(const char *name, int i);
static double *get_maskX(DATA **d, DPOINT *p, unsigned int row, unsigned int col);
static int get_map_location(DPOINT *loc, int random_path,
	unsigned int *row, unsigned int *col);
static DPOINT *get_point_location(int random_path);
static int get_random_cell(GRIDMAP *m, unsigned int *row, unsigned int *col);

static double *est = NULL;
static GRIDMAP **masks = NULL, **outmap = NULL;
static DATA *val_data = NULL;
static unsigned int n_done;
#ifdef WITH_SPIRAL
static DATA_GRIDMAP *mask_topology = NULL;
#endif

/* global variables: */
int strata_min;
unsigned int n_pred_locs = 0;
#define STRATUM(val) (floor(val - strata_min + 0.5))

void map_sign(GRIDMAP *m, const char *what);

#ifndef USING_R
void predict_all(DATA **data) {
	int i = 0, random_path = 0;
	DPOINT *here = NULL, *where = NULL;
	PRED_AT at_what;
	unsigned int row, col;

	n_done = 0;
	val_data = get_dataval();
	if (val_data->id > -1) {
		at_what = AT_POINTS;
		n_pred_locs = val_data->n_list;
		if (val_data->colns)
			strata_min = val_data->minstratum;
	} else if (get_n_masks() > 0) {
		at_what = AT_GRIDMAP;
		here = (DPOINT *) emalloc(sizeof(DPOINT));
		here->u.stratum = -2; /* only NON-MV cells */
		if (max_block_dimension(0) > 0.0)
			SET_BLOCK(here);
		else
			SET_POINT(here);
	} else /* what else ? */
		return;

	if (at_what == AT_GRIDMAP && get_n_outfile() == 0) {
		pr_warning("no output maps defined");
		return;
	}

	init_predictions(at_what);

	if (at_what == AT_GRIDMAP && !data[0]->dummy) { 
		if (data[0]->maxX < masks[0]->x_ul ||
				data[0]->minX > (masks[0]->x_ul + masks[0]->cols * masks[0]->cellsizex) ||
				data[0]->minY > masks[0]->y_ul ||
				data[0]->maxY < (masks[0]->y_ul - masks[0]->rows * masks[0]->cellsizey)) {
			pr_warning("ALL data are outside the map boundaries");
			printlog("data x[%g,%g], y[%g,%g]; map x[%g,%g], y[%g,%g]\n",
				data[0]->minX, data[0]->maxX, data[0]->minY, data[0]->maxY,
				masks[0]->x_ul,
				masks[0]->x_ul + masks[0]->cols * masks[0]->cellsizex,
				masks[0]->y_ul - masks[0]->rows * masks[0]->cellsizey,
				masks[0]->y_ul
				);
		}
		else if (map_xy2rowcol(masks[0], data[0]->minX, data[0]->minY, &row, &col) ||
				map_xy2rowcol(masks[0], data[0]->maxX, data[0]->minY, &row, &col) ||
				map_xy2rowcol(masks[0], data[0]->minX, data[0]->maxY, &row, &col) ||
				map_xy2rowcol(masks[0], data[0]->maxX, data[0]->maxY, &row, &col))
			pr_warning("at least some data are outside the map boundaries");
			/* this is not a sufficient test! */
	}

	if (gl_rp) /* Yes, by default */
		random_path = is_simulation(get_method());
	row = col = 0;
	while ((where = next_location(here, at_what, random_path,
			&row, &col, data)) != NULL) {
		for (i = 0; i < get_n_outfile(); i++)
			set_mv_double(&(est[i])); /* initialize estimates */
		if (where->u.stratum >= 0) {
			if (get_mode() != STRATIFY) {
				for (i = 0; i < get_n_vars(); i++)
					select_at(data[i], where);
			} else if (where->u.stratum < get_n_vars())
				select_at(data[where->u.stratum], where);
			get_est(data, get_method(), where, est);
		}
		/* printf("%g %g\n", est[0], est[1]); */
		write_output(est, at_what, where, row, col);
	}
	exit_predictions(at_what);
	if (here != NULL)
		efree(here);
	print_progress(100, 100);
}

static void init_predictions(PRED_AT w) {
	int i;
	DPOINT *bp;
	DATA **d = NULL;
#ifdef WITH_SPIRAL
	DATA_GRIDMAP *grid;
#endif

	est = (double *) emalloc(get_n_outfile() * sizeof(double));
	bp = get_block_p();
	d = get_gstat_data();
	switch (w) {
	  case AT_POINTS:
		if (o_filename == NULL) 
			ErrMsg(ER_VARNOTSET, "please specify output file");
		write_points(o_filename, val_data, NULL, NULL, get_n_outfile());
		if (bp->x == -1.0) { /* set default */
			bp->x = bp->y = 1.0;
			pr_warning("default block size set to: dx=1, dy=1");
		}
		break;
	  case AT_GRIDMAP:
		/* open mask files: */
		get_maskX(NULL, NULL, 0, 0); /* re-initializes static arrays */
		masks = (GRIDMAP **) emalloc(get_n_masks() * sizeof(GRIDMAP *));
		for (i = 0; i < get_n_masks(); i++)
			masks[i] = check_open(get_mask_name(i), i); /* read as float */
		if (n_pred_locs > 0)
			strata_min = floor(masks[0]->cellmin);
		outmap = (GRIDMAP **) emalloc(get_n_outfile() * sizeof(GRIDMAP *));
		printlog("initializing maps ");
		for (i = 0; i < get_n_outfile(); i++) {
			if (get_outfile_namei(i) != NULL) {
				printlog("."); /* creating maps ..... */
				if (get_method() == ISI)
					masks[0]->celltype = CT_UINT8;
				outmap[i] = map_dup(get_outfile_namei(i), masks[0]);
			} else
				outmap[i] = NULL;
		}
		printlog("\n");
		if (bp->x == -1.0) { /* set default to map cellsize */
			bp->x = masks[0]->cellsizex;
			bp->y = masks[0]->cellsizey;
			pr_warning("default block size set to dx=%g, dy=%g", bp->x, bp->y);
		}
		for (i = 0; i < get_n_vars(); i++) {
 			if (d[i]->dummy) {
 				d[i]->minX = masks[0]->x_ul + 0.5 * masks[0]->cellsizex;
 				d[i]->maxX = masks[0]->x_ul + masks[0]->cellsizex * 
						(masks[0]->cols - 0.5);
 				d[i]->maxY = masks[0]->y_ul - 0.5 * masks[0]->cellsizey;
 				d[i]->minY = masks[0]->y_ul - masks[0]->cellsizey * 
						(masks[0]->rows - 0.5);
 				d[i]->minZ = d[i]->maxZ = 0.0;
 			} 
 			if (d[i]->togrid)
 				datagrid_rebuild(d[i], 1);
		}
		break;
	}
	if (gl_nsim > 1)
		init_simulations(d);
	if (is_simulation(get_method()) && get_n_beta_set() != get_n_vars())
		setup_beta(d, get_n_vars(), gl_nsim);
} /* init_predictions() */

static void exit_predictions(PRED_AT what) {
	int i;

	if (gl_nsim > 1 && gl_lhs)
		lhs(get_gstat_data(), get_n_vars(), get_mode() == STRATIFY);

	switch (what) {
		case AT_POINTS:
			write_points(NULL, NULL, NULL, NULL, 0);
			if (gl_nsim > 1)
				save_simulations_to_ascii(o_filename);
			break;
		case AT_GRIDMAP:
			if (gl_nsim > 1) {
				if (DEBUG_DUMP)
					printlog("\nWriting results to files...");
				save_simulations_to_maps(masks[0]);
				if (DEBUG_DUMP)
					printlog("done");
			} else  {
				for (i = 0; i < get_n_outfile(); i++) {
					if (get_outfile_namei(i)) {
						map_sign(outmap[i], what_is_outfile(i));
						(outmap[i])->write(outmap[i]);
						map_free(outmap[i]);
					}
				}
			}
			for (i = 0; i < get_n_masks(); i++)
				map_free(masks[i]);
			efree(masks);
			efree(outmap);
			break;
	}
	print_orvc();
	efree(est);
} /* exit_predictions() */

static DPOINT *next_location(DPOINT *loc, PRED_AT what, int random_path,
		unsigned int *row, unsigned int *col, DATA **data) {

	double xc, yc;
	static unsigned int nr = 0;

	switch (what) {
		case AT_POINTS:
			if (DEBUG_TRACE) {
				nr++;
				printlog("\rbusy with loc: %3u", nr);
			}
			return get_point_location(random_path);
		case AT_GRIDMAP:
			if (get_map_location(loc, random_path, row, col)) {
				if (loc->u.stratum >= 0) { /* i.e., non-missing valued cell */
					if (DEBUG_TRACE)
						printlog("\rbusy with row: %3u col: %3u loc: %3u",
							*row + 1, *col + 1, nr + 1);
					map_rowcol2xy(masks[0], *row, *col, &xc, &yc);
					loc->x = xc;
					loc->y = yc;
					if (!is_mv_double(&gl_zmap))
						loc->z = gl_zmap;
					else
						loc->z = 0.0;
					loc->X = get_maskX(data, loc, *row, *col);
					nr++;
				}
				return loc;
			} 
			break;
	}
	return NULL;
} /* next_location() */ 

static void write_output(double *est, PRED_AT w, DPOINT *here,
		 unsigned int row, unsigned int col) {
	int i;

	switch (w) {
		case AT_POINTS:
			write_points(o_filename, val_data, here, est,
				get_mode() != STRATIFY ? get_n_outfile() : 2);
			break;
		case AT_GRIDMAP:
			for (i = 0; i < get_n_outfile(); i++)
				if (outmap[i] && !is_mv_double(&(est[i])))
					map_put_cell(outmap[i], row, col, est[i]);
			break;
	}
} /* write_output() */

static double *get_maskX(DATA **data, DPOINT *p,
		unsigned int row, unsigned int col) {
	static double *d = NULL;
	static int totX = 0, *posMask = NULL;
	int i, j, k, l;
	static DATA *bl = NULL;
	static GRIDMAP **local_masks = NULL;

	if (data == NULL) {
		if (d != NULL) {
			efree(d);
			d = NULL;
			efree(posMask);
			posMask = NULL;
			local_masks = NULL;
			totX = 0;
		}
		return NULL;
	}
	if (d == NULL) { /* first time calling */
		for (i = 0, totX = 0; i < get_n_vars(); i++) 
			totX += data[i]->n_X;
		posMask = (int *) emalloc(totX * sizeof(int));
		d = (double *) emalloc(totX * sizeof(double));
		for (i = 0, k = 0; i < get_n_vars(); i++) {
			for (j = 0; j < data[i]->n_X; j++) {
				if (data[i]->colX[j] > 0)
					posMask[k] = 1;
				if (data[i]->colX[j] == 0)
					posMask[k] = 0;
				if (data[i]->colX[j] < -1)
					posMask[k] = -1;
				k++;
			}
		}
		if (get_mode() == STRATIFY) {
			if (get_n_masks() > 1)
				local_masks = masks + 1; /* skip the first (= strata) map */
		} else
			local_masks = masks;
	}

 	bl = block_discr(bl, get_block_p(), p);
 	/* bl is a single point-list with p if IS_POINT(p) */
	for (i = 0, k = 0; i < get_n_vars(); i++) {
		for (j = 0; j < data[i]->n_X; j++) {
			if (data[i]->colX[j] < -1) {
				/* do eventual block averaging here: */
				for (l = 0, d[k] = 0.0; l < bl->n_list; l++)
					d[k] += bl->list[l]->u.weight *
							calc_polynomial(bl->list[l], data[i]->colX[j]);
			}
			k++;
		}
	}

	for (i = 0, j = 0; i < totX; i++) {
		switch(posMask[i]) {
			case -1:
				/* is done above */
				break;
			case 1:
				if (map_cell_is_mv(local_masks[j], row, col))
					ErrMsg(ER_IMPOSVAL, "missing value in one of the mask maps");
				d[i] = map_get_cell(local_masks[j], row, col);
				j++;
				break;
			case 0:
				d[i] = (double) 1.0;
				break;
			default:
				ErrMsg(ER_IMPOSVAL, "get_maskX()");
				break;
		}
	}
	return d;
}

int is_valid_strata_map(const char *name, int n_vars) {
	GRIDMAP *mp;
	int check_failed = 1;

	mp = new_map(READ_ONLY);
	mp->filename = name;
	if ((mp = map_read(mp)) == NULL)
		ErrMsg(ER_READ, name);
	/* 
	 * check min/max: enough strata? this check was:
	 * check_failed = (mp->cellmax - mp->cellmin < n_vars - 1);
	 */
	check_failed = (mp->cellmax - mp->cellmin < 1.0);
	if (! check_failed && mp->cellmax - mp->cellmin < n_vars - 1)
		pr_warning("fewer mask map categories than data variables present");
	map_free(mp);
	return (check_failed == 0);
}
#endif

unsigned int *get_n_sim_locs_table(unsigned int *size) {
	unsigned int i, j, *table;

	if (get_mode() == STRATIFY) {
		if (val_data->id > -1) {
			*size = val_data->maxstratum - val_data->minstratum + 2;
			table = (unsigned int *) emalloc(*size * sizeof(int));
			for (i = 0; i < *size; i++)
				table[i] = 0;
			for (i = 0; i < val_data->n_list; i++)
				table[val_data->list[i]->u.stratum]++;
		} else {
			*size = (int) STRATUM(masks[0]->cellmax) - 
					(int) STRATUM(masks[0]->cellmin) + 1;
			table = (unsigned int *) emalloc(*size * sizeof(int));
			for (i = 0; i < *size; i++)
				table[i] = 0;
			for (i = 0; i < masks[0]->rows; i++) {
				for (j = 0; j < masks[0]->cols; j++) {
					if (!map_cell_is_mv(masks[0], i, j))
						table[(int) STRATUM(map_get_cell(masks[0], i, j))]++;
				}
			}
		}
	} else {
		*size = (int) get_n_vars();
		table = (unsigned int *) emalloc(*size * sizeof(int));
		for (i = 0; i < *size; i++)
			table[i] = n_pred_locs;
	}
	return table;
}

#ifndef USING_R
static GRIDMAP *check_open(const char *name, int i) {
	GRIDMAP *mask;
	unsigned int r, c;

	mask = new_map(READ_ONLY);
	mask->filename = name;
	if ((mask = map_read(mask)) == NULL)
		ErrMsg(ER_READ, name);
	if (i == 0) {
		for (r = n_pred_locs = 0; r < mask->rows; r++)
			for (c = 0; c < mask->cols; c++)
				if (! map_cell_is_mv(mask, r, c))
					n_pred_locs++;
	} else if (! map_equal(mask, masks[0]))
		ErrMsg(ER_NULL, "mask map location is different from first mask map");
	return mask;
}

static DPOINT *get_point_location(int random_path) {
	static int current = 0;
	int i = 0, ri = 0; /* ri: random index */
	DPOINT *pt = NULL;

	if (current == val_data->n_list) {
		current = 0; /* reset for next run */
		return NULL;
	}
	if (current == 0 && random_path) { /* first time: randomize list order */
		for (i = 0; i < val_data->n_list; i++) {
			ri = floor(r_uniform() * (val_data->n_list));
			if (ri >= val_data->n_list) /* obsolete, but anyway... */
				ri = val_data->n_list - 1;
			/* now swap list pointers i and ri: */
			pt = val_data->list[i];
			val_data->list[i] = val_data->list[ri];
			val_data->list[ri] = pt;
		}
	} 
	if (DEBUG_TRACE)
		printf("[cell %d]\n", current); 
	pt = val_data->list[current];
	print_progress(current, val_data->n_list);
	SET_INDEX(pt, current);
	current++;
	return pt;
}

static int get_map_location(DPOINT *loc, int random_path,
	unsigned int *row, unsigned int *col) {
/* 
 * return a not-yet-visited map cell (in *row, *col), or
 * 0 if all cells are visited
 */
	double value;
	int at_end = 0;

	if (! random_path) {
		if (loc->u.stratum != -2) { /* don't increase when still at (0,0): */
			if (*col < masks[0]->cols - 1) /* goto next col in this row: */
					(*col)++;
			else if (*row < masks[0]->rows - 1) { /* start at next row */
				*col = 0;
				(*row)++;
			} else /* at end of map: */
				at_end = 1;
		}
	} else 
		at_end = (get_random_cell(masks[0], row, col) == 0);

	loc->u.stratum = -1; /* in case of MV in mask map: skip get_est() */
	if (!at_end && !map_cell_is_mv(masks[0], *row, *col)) {
		value = map_get_cell(masks[0], *row, *col);
		if (((int) STRATUM(value)) >= 0) 
			loc->u.stratum = (int) STRATUM(value);
		else {
			if (get_mode() == STRATIFY)
				pr_warning("negative stratum value %g set to zero", value);
			loc->u.stratum = 0;
		}
		est[0] = value;
		if (! DEBUG_TRACE)
			print_progress(n_done++, n_pred_locs);
	} 
	return (at_end == 0);
}
#endif

static int get_random_cell(GRIDMAP *m, unsigned int *row, unsigned int *col) {
	static char **u = NULL, *tmp = NULL;
	static unsigned int index = 0, i = 0, j = 0, p, q = 0;

/* 
 * SO, -- what's happening here?
 * We start with a coarse grid of the subset of grid of the smallest 2-power
 * that covers the whole map. E.g. at a 60x40 map, we first randomly locate
 * a 64x64 grid (having cells 64 map-cells apart), then we randomly locate
 * a 32x32 grid (...), than a 16x16, ... and finally a 1x1 grid (having
 * all cells not yet visited)
 */
	if (u == NULL) { /* first time -- allocate space */
		p = (unsigned int) ceil(log((double) MAX(m->rows, m->cols))/log(2.0));
		if (p > 32)
			ErrMsg(ER_RANGE, "map too large for GS/IS"); /* Whoow ;-) */
		u = (char **) emalloc(m->rows * sizeof(char *));
		tmp = (char *) emalloc(m->cols * m->rows * sizeof(char));
		for (i = 0; i < m->rows; i++)
			u[i] = &(tmp[i * m->cols]);
		for (i = 0; i < m->rows; i++) {
			for (j = 0; j < m->cols; j++)
				u[i][j] = 1; /* loop over all cells, also MV's */
		}
		i = 1;
		p = i << p; /* p = 2 ^ p */
		q = MIN(2, p);
		i = 0;
	}
	/* 
	 * find random path:
	 */
	while (q <= p) {
		/* find one random starting point: */
		if (i == 0)
			index = (unsigned int) floor(r_uniform() * (q * q - 1));
		while (i < q * q) {
			/* 
			 * get next point, use random path with every point visited once:
			 * therefore loop q * q times:
			 */
			index = (5 * index + 1) % (q * q); 
			*col = (index % q) * (p / q); /* colnr */
			*row = (index / q) * (p / q); /* rownr */
			if (*col < m->cols && *row < m->rows && u[*row][*col]) {
				i++;
				if (i == q * q) {
					q *= 2;
					i = 0;
				}
				u[*row][*col] = 0; /* skip this one the next time */
				return 1;
			} /* if *col .. */
			i++; /* loop i */
		} 
		q *= 2; /* next q */
		i = 0;  /* initialise i for the next q loop */
	} 
	/* when we're here, we're finished : */
	efree(u);
	u = NULL; 
	efree(tmp);
	tmp = NULL;
	return 0;
}

/* 
 * procedure map_sign() for putting history and description to maps (csf)
 * history is fixed for a session; description depends on the map contents
 */
void map_sign(GRIDMAP *m, const char *what) {
	char *pwd = NULL, *user = NULL, *timestring = NULL, *fname = NULL;
	time_t tm;

	if ((user = getenv("LOGNAME")) == NULL)
		user = ""; 
	if ((pwd = getenv("PWD")) == NULL)
		pwd = "";
	if ((fname = command_file_name) == NULL)
		fname = "";
	tm = time(&tm);
	if ((timestring = asctime(localtime(&tm))) == NULL)
		timestring = "";
	m->history = (char *) emalloc(1024 * sizeof(char));
#ifdef HAVE_SNPRINTF
	snprintf
#else 
	sprintf
#endif
		(m->history,
#ifdef HAVE_SNPRINTF
		1024,
#endif
		"%s%s\n%s %s %s\n%s%s\n%s %s\n%s %s\n%s %s (seed %lu)\n",
		*user ? 
		"creator:       " : "", user, 
		"program:      ", argv0, VERSION,
		*pwd ? 
		"path:          " : "", pwd,
		"command file: ", fname,
		"method used:  ", method_string(get_method()),
		"date:         ", timestring, get_seed());
	if (what == NULL) { /* paranoia */
		m->description = (char *) emalloc(sizeof(char));
		m->description[0] = '\0';
	} else {
		m->description = (char *) emalloc((strlen(what) + 1) * sizeof(char));
#ifdef HAVE_SNPRINTF
		snprintf(m->description, strlen(what) + 1, "%s\n", what);
#else
		sprintf(m->description, "%s\n", what);
#endif

	}
}

const void *get_mask0(void) {
	if (masks)
		return masks[0];
	return NULL;
}
