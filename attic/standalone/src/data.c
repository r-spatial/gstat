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
 *  data.c: basic i/o routines on DATA structure 
 */
#include <stdio.h>
#include <stdlib.h> /* qsort() */
#include <math.h>
#include <string.h>
#include <limits.h>
#include <float.h>

#include "defs.h"
#include "data.h"
#include "mapio.h"
#include "userio.h"
#include "utils.h"
#include "block.h"
#include "debug.h"
#include "read.h"
#include "glvars.h"
#include "polygon.h"
#include "random.h"
#include "defaults.h"
#include "matrix.h"
#include "lm.h" /* free_lm() */
#include "gls.h" /* free_glm() */
#include "nsearch.h"
#include "gcdist.h"

#ifdef HAVE_EXT_DBASE
#ifndef INCLUDED_EXT_DBASE
#include "ext_dbase.h"
#define INCLUDED_EXT_DBASE
#endif
#endif

const DATA_TYPE data_types[] = {
	{ DATA_UNKNOWN, "Unknown file type"},
	{ DATA_ASCII_TABLE, "Ascii table file"},
	{ DATA_EAS, "GeoEAS file"},
	{ DATA_IDRISI_VEC, "Idrisi ascii .vec"},
	{ DATA_IDRISI32_VEC, "Idrisi binary .vct"},
	{ DATA_IDRISI_BIN, "Idrisi binary .img"},
	{ DATA_IDRISI_ASCII, "Idrisi ascii .img"},
	{ DATA_IDRISI32_BIN, "Idrisi binary .rst"},
	{ DATA_IDRISI32_ASCII, "Idrisi ascii .rst"},
	{ DATA_GRIDASCII, "ArcInfo gridascii"},
	{ DATA_GRIDFLOAT, "ArcInfo gridfloat"},
	{ DATA_CSF, "PCRaster map"},
	{ DATA_T2, "T2 map"},
	{ DATA_ERMAPPER, "ER-Mapper file"},
	{ DATA_GNUPLOT, "Gnuplot binary"},
	{ DATA_GMT, "GMT netCDF format" },
	{ DATA_SURFER_DSAA, "Surfer DSAA ascii grid" },
	{ DATA_GSLIB, "GSLIB grid" },
	{ DATA_GRASS, "GRASS site list" },
	{ DATA_GRASS_GRID, "GRASS raster" },
	{ DATA_GDAL, "GDAL raster map" }
}; 

static int read_data_line(FILE *f, char *line, DATA *d, DPOINT *current,
	int *line_nr, int colmax);
static int is_one_integer(char *line, int *n);
static void init_dpoint(DATA *d, DPOINT *current);
static int double_is_mv(DATA *d, double *f);
static DATA *read_table(DATA *d);
static void mk_var_names(DATA *d);
static int read_eas_header(FILE *infile, DATA *d, int ncols);
static void calc_data_mean_std(DATA *d);
static void correct_strata(DATA *d);
static void field_error(char *fname, int line_nr, int fld, char *text);
static int average_duplicates(DATA *d);
#ifndef USING_R
static int read_data_from_map(DATA *d);
static int read_idrisi_points(DATA *d);
static int read_idrisi_point_data(DATA *d, const char *fname);
static int read_idrisi_point_header(DATA *d, const char *fname);
static int read_idrisi32_points(DATA *d);
static int read_idrisi32_point_data(DATA *d, const char *fname);
static int read_idrisi32_point_header(DATA *d, const char *fname);
#endif
static void transform_data(DATA *d);
static void grid_push_point(DATA *d, DPOINT *p, int adjust_to_gridcentrs);
static double point_norm_1D(const DPOINT *p);
static double point_norm_2D(const DPOINT *p);
static double point_norm_3D(const DPOINT *p);
static double pp_norm_1D(const DPOINT *a, const DPOINT *b);
static double pp_norm_2D(const DPOINT *a, const DPOINT *b);
static double pp_norm_3D(const DPOINT *a, const DPOINT *b);

/* great circle distances: */
static double point_norm_gc(const DPOINT *p);
static double pp_norm_gc2(const DPOINT *a, const DPOINT *b);
static double pb_norm_gc2(const DPOINT *where, BBOX bbox);

static void free_data_gridmap(DATA_GRIDMAP *t);
static void logprint_data_header(const DATA *d);

#ifdef HAVE_LIBGIS
#include "gis.h"
#include "site.h"
static DATA *read_grass_data(DATA * d);
#endif

static DPOINT min, max;
static int fix_minmax = 0;
const POLY_NM polynomial[N_POLY] =
    {{ POLY_X, "x",     1, X_BIT_SET},
     { POLY_Y, "y",     1, Y_BIT_SET},
     { POLY_Z, "z",     1, Z_BIT_SET},
     { POLY_X2, "x2",   2, X_BIT_SET},
     { POLY_Y2, "y2",   2, Y_BIT_SET},
     { POLY_Z2, "z2",   2, Z_BIT_SET},
     { POLY_XY, "xy",   2, Y_BIT_SET},
     { POLY_XZ, "xz",   2, Z_BIT_SET},
     { POLY_YZ, "yz",   2, Z_BIT_SET},
     { POLY_X3, "x3",   3, X_BIT_SET},
     { POLY_Y3, "y3",   3, Y_BIT_SET},
     { POLY_Z3, "z3",   3, Z_BIT_SET},
     { POLY_X2Y, "x2y", 3, Y_BIT_SET},
     { POLY_XY2, "xy2", 3, Y_BIT_SET},
     { POLY_X2Z, "x2z", 3, Z_BIT_SET},
     { POLY_XZ2, "xz2", 3, Z_BIT_SET},
     { POLY_Y2Z, "y2z", 3, Z_BIT_SET},
     { POLY_YZ2, "yz2", 3, Z_BIT_SET}};

/*
 * read_gstat_data():
 * Purpose       : read in data (x,y,z,value, ..) from a ASCII column file 
 * Created by    : Edzer J. Pebesma                                       
 * Date          : 6 april 1992                                          
 * Prerequisites : d,     x > 0, y > 0 , x != y, ..                     
 * Returns       : struct DATA (data.h)                                
 * Side effects  : emalloc's some memory                              
 *
 * function returns struct DATA (defined in data.h) 
 * reads from each line in file d->fname colums x, y, z, attr,
 * read_data allocates the neccesary memory for DATA.list;
 * if d->fname is "-", data is read from stdin (with a warning message);
 * in this case, it is assumed to be in GeoEAS format.
 */
#ifndef USING_R
DATA *read_gstat_data(DATA *d) {
	int check = 0;
	DATA *valdata;

	valdata = get_dataval();
	if (d == NULL)
		ErrMsg(ER_NULL, "ReadData()");
	d->minX = d->maxX = 0.0; 
	d->minY = d->maxY = 0.0;
	d->minZ = d->maxZ = 0.0;

	if (d->log)
		check++;
	if (!is_mv_double(&(d->Icutoff)))
		check++;
	if (d->nscore_table != NULL)
		check++;
	if (d->Category)
		check++;
	if (check > 1)
		ErrMsg(ER_IMPOSVAL, 
			"choose only one from log, indicator or nscore transform");

	if (d->id < 0) 
		ErrMsg(ER_IMPOSVAL, "read_gstat_data(): id < 0");
	if (d->fname == NULL || d->fname[0] == '\0' || almost_equals(d->fname, " "))
		d->dummy = 1;
	if (d->dummy) {
		d->n_list = d->n_sel = d->n_max = d->n_original = 0;
		/* determine mode: */
		if (get_n_masks() > 0)
			d->mode = X_BIT_SET | Y_BIT_SET | V_BIT_SET;
		else if (valdata->id != -1) {
			d->mode = V_BIT_SET;
			if (valdata->colnx)
				d->mode = d->mode | X_BIT_SET;
			if (valdata->colny)
				d->mode = d->mode | Y_BIT_SET;
			if (valdata->colnz)
				d->mode = d->mode | Z_BIT_SET;
		} else {
			if (d->colnx)
				d->mode = d->mode | X_BIT_SET;
			if (d->colny)
				d->mode = d->mode | Y_BIT_SET;
			if (d->colnz)
				d->mode = d->mode | Z_BIT_SET;
			if (d->mode == 0)
				ErrMsg(ER_IMPOSVAL, "cannot determine coordinate mode of data");
		}
		set_norm_fns(d);
		return d;
	} 
	/* CW trying is done by special fname format
	 * so it must be tried first
	 */
#ifdef HAVE_EXT_DBASE
	if (strncmp(d->fname,EXT_DBASE_FNAME_SIG,strlen(EXT_DBASE_FNAME_SIG))==0) {
		read_ext_dbase(d);
	} else
#endif
#ifdef HAVE_LIBGIS
	if (grass()) {
		DUMP(d->fname); 
		DUMP(": trying Grass site list/raster ... ");
		if (read_grass_data(d)) {
			DUMP("Yes; site list\n");
		} else if (read_data_from_map(d)) {
			DUMP("Yes; raster\n");
		} else {
			DUMP("no\n");
			ErrMsg(ER_READ, d->fname);
		}
	} else
#endif
	if (!read_data_from_map(d) && !read_idrisi_points(d) && 
				!read_idrisi32_points(d)) {
		/* last try assumes ascii/with geo-eas header */
		DUMP(d->fname); DUMP(": trying GeoEAS or ascii table... ");
		d = read_table(d);
		if (d->type.type == DATA_EAS) {
			DUMP("GeoEAS\n");
		} else {
			DUMP("ascii\n");
		}
	}
	set_norm_fns(d);
	d->n_original = d->n_list;
    if (d->type.type != DATA_EXT_DBASE) {
	 /* CW only access data to get some info
	  * not needed for current application of DATA_EXT_DBASE
	  * or read_ext_dbase should compute/return that info
	  */
		transform_data(d); 
		average_duplicates(d);
		calc_data_mean_std(d);
		correct_strata(d);
    }
	return d;
}
#endif

void set_norm_fns(DATA *d) {

	if (d->mode & Z_BIT_SET) {
		d->point_norm = point_norm_3D;
		d->pp_norm2 = pp_norm_3D;
		d->pb_norm2 = pb_norm_3D;
	} else if (d->mode & Y_BIT_SET) {
		if (gl_longlat) {
			d->point_norm = point_norm_gc;
			d->pp_norm2 = pp_norm_gc2;
			d->pb_norm2 = pb_norm_gc2;
			/* if (gl_split != DEF_split)
				pr_warning("longlat data cannot do quadtree, setting split to %d", INT_MAX); */
			gl_split = INT_MAX;
		} else {
			d->point_norm = point_norm_2D;
			d->pp_norm2 = pp_norm_2D;
			d->pb_norm2 = pb_norm_2D;
		}
	} else {
		d->point_norm = point_norm_1D;
		d->pp_norm2 = pp_norm_1D;
		d->pb_norm2 = pb_norm_1D;
	}
}

void init_data_minmax(void) {
	fix_minmax = 0;
	set_mv_double(&(min.x));
	set_mv_double(&(min.y));
	set_mv_double(&(min.z));
	set_mv_double(&(max.x));
	set_mv_double(&(max.y));
	set_mv_double(&(max.z));
}

void setup_data_minmax(DATA *d) {

	if (fix_minmax)
		ErrMsg(ER_NULL, "min and max should be fixed");

	if (d->id == 0) {
		min.x = d->minX; min.y = d->minY; min.z = d->minZ;
		max.x = d->maxX; max.y = d->maxY; max.z = d->maxZ;
	} else {
		min.x = MIN(min.x, d->minX);
		min.y = MIN(min.y, d->minY);
		min.z = MIN(min.z, d->minZ);
		max.x = MAX(max.x, d->maxX);
		max.y = MAX(max.y, d->maxY);
		max.z = MAX(max.z, d->maxZ);
	}
}

#ifndef USING_R
static DATA *read_table(DATA *d) {
	int colmax = 0, i, line_size = 0, line2_size = 0, line_nr = 1, ncols;
	static DPOINT current;
	static int sizeof_currentX = 0;
	char *line = NULL, *line2 = NULL;
	FILE *infile = NULL;

	current.u.stratum = 0;
	if (sizeof_currentX == 0)
		current.X = NULL;
	colmax = MAX(d->colnx, MAX(d->colny, MAX(d->colnz, d->colnvalue)));
	colmax = MAX(colmax, MAX(d->colnvariance, d->colns));
	colmax = MAX(colmax, d->coln_id);
    
	if (d->colnx > 0)
		d->mode = d->mode | X_BIT_SET;
	if (d->colny > 0)
		d->mode = d->mode | Y_BIT_SET;
	if (d->colnz > 0)
		d->mode = d->mode | Z_BIT_SET;
	if (d->colns > 0) 
		d->mode = d->mode | S_BIT_SET;
	if (d->colnvalue > 0) 
		d->mode = d->mode | V_BIT_SET;
	setup_polynomial_X(d);
	/* UK, LM: */
	for (i = 0; i < d->n_X; i++)
		colmax = MAX(colmax, d->colX[i]);
	if (colmax < 1)
		ErrMsg(ER_VARNOTSET, "please specify data columns");

	infile = efopen(d->fname, "r");
	if (get_line(&line, &line_size, infile) == NULL) {
		pr_warning("empty data file: %s", d->fname);
		ErrMsg(ER_READ, d->fname);
	}
	line_nr++;
	d->type = data_types[DATA_ASCII_TABLE]; /* try to disprove: */
	if (get_line(&line2, &line2_size, infile) != NULL) { /* check on header */
		line_nr++;
		if (!is_one_integer(line, &ncols) && is_one_integer(line2, &ncols)) {
			d->type = data_types[DATA_EAS];
			line_nr += read_eas_header(infile, d, ncols); /* checks errors */
		} 
	}
	mk_var_names(d);

	d->n_list = d->n_max = 0; 
	if (d->id == ID_OF_VALDATA && max_block_dimension(0) > 0.0) {
		SET_BLOCK(&current);
		if (IS_POINT(&current)) {
			pr_warning("bitfield:[%u]\n", current.bitfield);
			ErrMsg(ER_IMPOSVAL, "IS_BLOCK");
		}
	} else
		SET_POINT(&current);
	/* X's */
	if (d->n_X > sizeof_currentX) { /* and therefore > 0: */
		if (sizeof_currentX == 0)
			current.X = (double *) emalloc(d->n_X * sizeof(double));
		else
			current.X = (double *) erealloc(current.X, d->n_X * sizeof(double));
		sizeof_currentX = d->n_X;
	}
	if (d->type.type == DATA_ASCII_TABLE) { /* re-read first two lines: */
		line_nr -= 2;
		if (! read_data_line(infile, line, d, &current, &line_nr, colmax))
			push_point(d, &current);
		if (line2_size && !read_data_line(infile, line2, d, &current, &line_nr, colmax)) {
			push_point(d, &current);
			efree(line2);
		}
	}
	while (get_line(&line, &line_size, infile) != NULL) {
		if (! read_data_line(infile, line, d, &current, &line_nr, colmax))
			push_point(d, &current);
	} /* while get_line() */
	if (infile != (FILE *) stdin)
		efclose(infile);
	efree(line);
	return d;
}
#endif

static void mk_var_names(DATA *d) {
	char tmp1[100], tmp2[100]; 
	/* make variable names, if not read from EAS header: */
	if (d->variable == NULL || d->variable[0] == '\0')
		sprintf(tmp1, "col[%d]", d->colnvalue);
	else {
		if (strlen(d->variable) > 50)
			d->variable[50] = '\0';
		sprintf(tmp1, "%s", d->variable);
	}
	if (d->log)
		sprintf(tmp2, "log(%s)", tmp1);
	else if (!is_mv_double(&(d->Icutoff)))
		sprintf(tmp2, "I(%s,%g)", tmp1, d->Icutoff);
	else if (d->Category)
		sprintf(tmp2, "I(%s='%s')", tmp1, d->Category);
	else if (d->nscore_table != NULL)
		sprintf(tmp2, "%s (nscores)", tmp1);
	else /* no transform, copy: */
		sprintf(tmp2, "%s", tmp1);
	d->variable = string_dup(tmp2);
	if (d->x_coord == NULL && d->colnx) {
		sprintf(tmp1, "x_%d", d->colnx);
		d->x_coord = string_dup(tmp1);
	}
	if (d->y_coord == NULL && d->colny) {
		sprintf(tmp1, "y_%d", d->colny);
		d->y_coord = string_dup(tmp1);
	}
	if (d->z_coord == NULL && d->colnz)  {
		sprintf(tmp1, "z_%d", d->colnz);
		d->z_coord = string_dup(tmp1);
	}
	if (d->V_coord == NULL && d->colnvariance)  {
		sprintf(tmp1, "var_%d", d->colnvariance);
		d->V_coord = string_dup(tmp1);
	}
	if (d->s_coord == NULL && d->colns) {
		sprintf(tmp1, "stra_%d", d->colns);
		d->s_coord = string_dup(tmp1);
	}
    if (d->id_name == NULL && d->coln_id) {
        sprintf(tmp1, "ID_%d", d->colns);
        d->id_name = string_dup(tmp1);
    }
}

static void init_dpoint(DATA *d, DPOINT *current) {
	int i;

	current->x = current->y = current->z = 0.0;
	current->attr = current->variance = 0.0;
    
	for (i = 0; i < d->n_X; i++) {
		if (d->colX[i] == 0) 
			current->X[i] = 1.0;
		else
			current->X[i] = 0.0;
	}
}

#ifndef USING_R
static int read_data_line(FILE *f, char *line, DATA *d, DPOINT *current,
		int *line_nr, int colmax) {
/*
 * reads data from line;
 * if this one is a missing value, return 1, or else 0.
 */
	int this_one_is_mv, field, line_size = 0, i;
	char *token = NULL;
	char *category = NULL;

	init_dpoint(d, current);
	if (line[0] == '#') {
		this_one_is_mv = 1; /* skip this line */
		pr_warning("skipping line %d in data file %s", *line_nr, d->fname);
	} else
		this_one_is_mv = 0;
	(*line_nr)++;
	for (field = 1; field <= colmax && !this_one_is_mv; field++) {
		token = strtok(field == 1 ? line : NULL, DELIMITERS);
		if (token == NULL && field == 1) { /* empty line: at the bottom? */
			while (get_line(&line, &line_size, f) != NULL) {
				if (strtok(line, DELIMITERS) != NULL) { /* no, there's more: */
					message("%s:%d: too few records on line\n", d->fname, *line_nr);
					ErrMsg(ER_READ, (const char *) d->fname);
				}
			}
			this_one_is_mv = 1; /* skip this line */
			/* if we do get out of this while(), break from the for loop: */ 
			break;
		} else if (token == NULL) {
			message("%s:%d: too few records on line\n", d->fname, *line_nr);
			ErrMsg(ER_READ, (const char *) d->fname);
		}
		if (field == d->colnx) { 
			if (read_double(token, &(current->x)))
				field_error(d->fname, *line_nr, field, token);
			if (! this_one_is_mv)
				this_one_is_mv = double_is_mv(d, &(current->x));
		}
		if (field == d->colny) {
			if (read_double(token, &(current->y)))
				field_error(d->fname, *line_nr, field, token);
			if (! this_one_is_mv)
				this_one_is_mv = double_is_mv(d, &(current->y));
		}
		if (field == d->colnz) {
			if (read_double(token, &(current->z)))
				field_error(d->fname, *line_nr, field, token);
			if (! this_one_is_mv)
				this_one_is_mv = double_is_mv(d, &(current->z));
		}
		if (field == d->colnvariance) {
			if (read_double(token, &(current->variance)))
				field_error(d->fname, *line_nr, field, token);
			if (! this_one_is_mv)
				this_one_is_mv = double_is_mv(d, &(current->variance));
			if (! this_one_is_mv && current->variance <= 0.0)
				ErrMsg(ER_IMPOSVAL, "only positive variances allowed");
		}
		if (field == d->colnvalue) { 
			if (d->Category)
				category = string_dup(token);
			else {
				if (read_double(token, &(current->attr)))
					field_error(d->fname, *line_nr, field, token);
				if (! this_one_is_mv) 
					this_one_is_mv = double_is_mv(d, &(current->attr));
			}
		}
		if (field == d->colns) {
			if (read_int(token, &(current->u.stratum)))
				field_error(d->fname, *line_nr, field, token);
		}
        if (field == d->coln_id) {
            d->point_ids=erealloc(d->point_ids, (d->n_list+1)*sizeof(char**));
            d->point_ids[d->n_list] = string_dup(token);
        }

		/* uk, lm: */
		if (d->n_X > 0) {
			for (i = 0; i < d->n_X; i++) 
				if (field == d->colX[i]) { /* field is ge 1 */
					if (read_double(token, &(current->X[i])))
						field_error(d->fname, *line_nr, field, token);
					if (! this_one_is_mv) 
						this_one_is_mv = double_is_mv(d, &(current->X[i]));
				}
		}
	} /* for field ... */
	if (! this_one_is_mv) { /* nothing misses: */
		if (d->n_list == 0) {
			d->minX = d->maxX = current->x; 
			d->minY = d->maxY = current->y; 
			d->minZ = d->maxZ = current->z; 
			d->minvariance = d->maxvariance = current->variance;
			d->minstratum = d->maxstratum = current->u.stratum;
		} else {
			d->minX = MIN(d->minX, current->x); 
			d->maxX = MAX(d->maxX, current->x);
			d->minY = MIN(d->minY, current->y); 
			d->maxY = MAX(d->maxY, current->y); 
			d->minZ = MIN(d->minZ, current->z); 
			d->maxZ = MAX(d->maxZ, current->z);
			d->minvariance = MIN(d->minvariance, current->variance);
			d->maxvariance = MAX(d->maxvariance, current->variance);
			d->minstratum = MIN(d->minstratum, current->u.stratum);
			d->maxstratum = MAX(d->maxstratum, current->u.stratum);
		}
		if (category != NULL) {
			current->attr = strcmp(d->Category, category) ? 0 : 1;
			efree(category);
			category = NULL;
		}
	} 
	return this_one_is_mv;
} 
#endif

static int double_is_mv(DATA *d, double *f) {
	if (is_mv_double(&(d->mv)))
		return (is_mv_double(f));
	else 
		return ((is_mv_double(f)) || (*f == d->mv));
}

DATA *get_area_centre(DATA *area, DATA *d) {
	int i, j;
	DPOINT p;

	d->n_list = d->n_max = 0; 
	d->variable = area->variable;
	d->x_coord = area->x_coord;
	d->y_coord = area->y_coord;
	d->z_coord = area->z_coord;
	d->type = data_types[area->type.type];
	d->fname = "";

	p.x = p.y = p.z = 0.0;
	p.u.stratum = 0;
	d->n_X = area->n_X;
	if (area->n_X > 0) {
		p.X = (double *) emalloc(area->n_X * sizeof(double));
		d->colX = (int *) emalloc(area->n_X * sizeof(int));
		for (j = 0; j < area->n_X; j++) {
			p.X[j] = 0.0;
			d->colX[j] = area->colX[j];
		}
	} else {
		p.X = NULL;
		d->colX = NULL;
	}
	for (i = 0; i < area->n_list; i++) {
		p.x += area->list[i]->x;
		p.y += area->list[i]->y;
		p.z += area->list[i]->z;
		for (j = 0; j < area->n_X; j++)
			p.X[j] += area->list[i]->X[j]; 
	}
	p.x /= area->n_list;
	p.y /= area->n_list;
	p.z /= area->n_list;
	for (j = 0; j < area->n_X; j++)
		p.X[j] /= area->n_list;
	p.attr = 0.0;

	printlog("prediction centre at x=%g, y=%g, z=%g",p.x,p.y,p.z);
	if (d->n_X) {
		printlog(" where x0 averages [");
		for (j = 0; j < area->n_X; j++)
			printlog("%g%s", p.X[j], j<area->n_X-1?",":"");
		printlog("]\n");
	} else
		printlog("\n");
	push_point(d, &p);
	d->minX = d->maxX = p.x;
	d->minY = d->maxY = p.y;
	d->minZ = d->maxZ = p.z;
	d->mode = area->mode;
	d->n_X = area->n_X;
	calc_data_mean_std(d);
	return d;
}

void centre_area(DATA *area) {
	int i;
	DPOINT p;

	p.x = p.y = p.z = 0.0;
	for (i = 0; i < area->n_list; i++) {
		p.x += area->list[i]->x;
		p.y += area->list[i]->y;
		p.z += area->list[i]->z;
	}
	p.x /= area->n_list;
	p.y /= area->n_list;
	p.z /= area->n_list;
	for (i = 0; i < area->n_list; i++) {
		area->list[i]->x -= p.x;
		area->list[i]->y -= p.y;
		area->list[i]->z -= p.z;
	}
	area->minX -= p.x;
	area->maxX -= p.x;
	area->minY -= p.y;
	area->maxY -= p.y;
	area->minZ -= p.z;
	area->maxZ -= p.z;
}

void report_data(const DATA *d) {
	int i, j;
	char tmp[11];
/*
 * OUTPUT information on data read: 
 */

 	tmp[10] = '\0';

 	if (d->dummy) {
 		printlog("[dummy variable]\n");
 		return;
 	}
/*
gstat 2.0a (December 1997) for Linux
Copyright (C) 1992, 1997 Edzer J. Pebesma
data(a):                  xx     (GeoEAS file)
attribute:         zinc, ppm     [x:] xcoord, m  : [    178605,    181390]
n:                       155     [y:] ycoord, m  : [    329714,    333611]
sample mean:         469.716     [z:] zinc, ppm  : [       113,      1839]
sample std.:         367.074     [V:] zinc, ppm  : [       113,      1839]
data():                  xx     (GeoEAS file)
attribute:            col. 0     [x:] xcoord, m  : [    178605,    181390]
n:                       155     [y:] ycoord, m  : [    329714,    333611]
                                 [z:] zinc, ppm  : [       113,      1839]
                                 [V:] zinc, ppm  : [       113,      1839]
                                 [s:] zinc, ppm  [       113,      1839]
*/
	i = (d->id == ID_OF_VALDATA ? 0 : strlen(name_identifier(d->id)));
	j = 21 - strlen(d->fname);
	for ( ; i < j; i++)
		printlog(" ");
	printlog("%s     (%s)\n", d->fname, d->type.name);
	printlog("attribute: %18s     ", NULS(d->variable));

	if (d->mode & X_BIT_SET)
		printlog("[x:] %-10s : [%10g,%10g]", 
			strncpy(tmp, NULS(d->x_coord), 10), d->minX, d->maxX);

	printlog("\n");

	printlog("n: ");
	if (d->n_averaged > 0)
		printlog(" [avgd. %4d] %12d     ", d->n_averaged, d->n_list); 
	else
		printlog("%26d     ", d->n_list);

	if (d->mode & Y_BIT_SET)
		printlog("[y:] %-10s : [%10g,%10g]",
			strncpy(tmp, NULS(d->y_coord), 10), d->minY, d->maxY);
	printlog("\n");

	if (d->mode & V_BIT_SET) {
		printlog("sample mean: %16g     ", d->mean);
		if (d->mode & Z_BIT_SET) {
			printlog("[z:] %-10s : [%10g,%10g]\n",
				strncpy(tmp, NULS(d->z_coord), 10), d->minZ, d->maxZ);
		} 
		printlog("sample std.: %16g     ", d->std);
		if (!d->colnvariance || !(d->mode & Z_BIT_SET))
			printlog("\n");
	}

	if ((d->mode & Z_BIT_SET) && !(d->mode & V_BIT_SET)) {
		printlog("%-33s", " ");
		printlog("[z:] %-10s : [%10g,%10g]\n",
			strncpy(tmp, NULS(d->z_coord), 10), d->minZ, d->maxZ);
	}

	if (d->colnvariance) {
		if (!(d->mode & V_BIT_SET))
			printlog("%-33s", " ");
			printlog("[V:] %-10s : [%10g,%10g]\n",
			strncpy(tmp, NULS(d->V_coord), 10), d->minvariance, d->maxvariance);
	}

	if (d->colns)
		printlog("%-33s[s:] %-10s : [%10d,%10d]\n", " ",
		strncpy(tmp, NULS(d->s_coord), 10), d->minstratum, d->maxstratum);

	if (d->n_X > 1 || (d->n_X == 1 && d->colX[0] != 0)) {
		printlog("base functions:"); 
		for (i = 0; i < d->n_X; i++) {
			printlog("%s", i ? ", " : " ");
			if (d->colX[i] <= 0) {
				if (d->colX[i] == 0)
					printlog("intercept");
				else
					printlog("%s", POLY_NAME(d->colX[i]));
			} else
				printlog("column %d", d->colX[i]);
		}
		printlog("\n"); 
	}
	/* fflush(dest); */
	return;
}

static void calc_data_mean_std(DATA *d) {
/* 
 * Calculates fields mean and std of d with mean and standard dev. (/(n-1))
 */
	int i;

	if (d->standard == 2) { /* we did this already.. */
		for (i = 0; i < d->n_list; i++)
			d->list[i]->attr *= d->std;
	}

	d->mean = 0.0; 
	d->std = 0.0;
	if (d->n_list <= 0) {
		pr_warning("calc_data_mean_std: n_list <= 0: %d", d->n_list);
		return;
	}
	for (i = 0; i < d->n_list; i++) 
		d->mean += d->list[i]->attr;
	d->mean /= d->n_list;
	if (d->n_list == 1)
		return;
	for (i = 0; i < d->n_list; i++) 
		d->std += SQR(d->list[i]->attr - d->mean);
	d->std = sqrt((d->std)/(d->n_list - 1));

	if (d->standard > 0) {
		for (i = 0; i < d->n_list; i++)
			d->list[i]->attr /= d->std;
		d->standard = 2;
	}
	return;
}

static void correct_strata(DATA *d) {
	int i;

	if (d->colns == 0)
		return;
	for (i = 0; i < d->n_list; i++)
		d->list[i]->u.stratum -= d->minstratum;
}

static void transform_data(DATA *d) {
	Double_index *values;
	double nscore, q;
	FILE *f = NULL;
	int i, j, tie_length, tie_total = 0;
	DPOINT *p;

	if (d->log) {
		for (i = 0; i < d->n_list; i++) {
			p = d->list[i];
			if (p->attr <= 0.0) 
				ErrMsg(ER_IMPOSVAL, "log of non-positive value");
			p->attr = log(p->attr);
		}
	} else if (! is_mv_double(&(d->Icutoff))) {
		for (i = 0; i < d->n_list; i++) {
			p = d->list[i];
			if (p->attr <= d->Icutoff)
				p->attr = 1.0;
			else
				p->attr = 0.0;
		}
	} else if (d->nscore_table) {
#ifndef USING_R
		if (strlen(d->nscore_table) > 0) {
			f = efopen(d->nscore_table, "w");
			fprintf(f, "normal scores for %s\n", d->variable);
			fprintf(f, "2\nobserved value\nnormal score\n");
		}
#endif
		values = emalloc(d->n_list * sizeof(Double_index));
		for (i = 0; i < d->n_list; i++) {
			values[i].d = d->list[i]->attr; /* */
			values[i].index = i; /* unique identifiers */
		}
		qsort((void *) values, (size_t) d->n_list, sizeof(Double_index),
			(int CDECL (*)(const void *, const void *)) double_index_cmp);
		for (i = 0; i < d->n_list; i += tie_length) { /* ignore ties: */
			/* assume they're all ties, with min run length 1: */
			tie_length = 1;
			while (i + tie_length < d->n_list && 
					values[i].d == values[i+tie_length].d)
				tie_length++;
			q = (i + 0.5 * tie_length) / d->n_list;
			nscore = q_normal(q);
			for (j = 0; j < tie_length; j++) {
#ifndef USING_R
				if (f)
					fprintf(f, "%g %g\n", d->list[values[i+j].index]->attr, 
							nscore);
#endif
				/* transform: */
				d->list[values[i+j].index]->attr = nscore;
			}
			tie_total += (tie_length - 1);
		}
		efree(values);
		if (f)
			efclose(f);
		if (tie_total)
			pr_warning("tied normal scores are assigned to tied data (%.2g%%)", 
					100.0 * tie_total / d->n_list);
	}
}

void setup_polynomial_X(DATA *d) {

	int i, j, degree;

	degree = d->polynomial_degree;
	if (degree < 0 || degree > 3)
		ErrMsg(ER_SYNTAX, "polynomial degree n, `d=n', should be in [0..3]");
	for (i = 1; i <= degree; i++)
		for (j = 0; j < N_POLY; j++)
			if (polynomial[j].degree == i && (d->mode & polynomial[j].mode))
				data_add_X(d, polynomial[j].poly_nr);
}

void data_add_X(DATA *d, int col) {
	int i;

	for (i = 0; d->id != ID_OF_VALDATA && i < d->n_X; i++)
		if (d->colX[i] == col)
			ErrMsg(ER_IMPOSVAL, "X-variable: column appears twice");
	d->n_X++;
	d->colX = (int *) erealloc(d->colX, d->n_X * sizeof(int));
	d->colX[d->n_X - 1] = col;
}

void calc_polynomials(DATA *d) {
	int i, j, do_block;

#define CHECK_BITX if(!(d->mode & X_BIT_SET)) ErrMsg(ER_VARNOTSET,"x coordinate not set")
#define CHECK_BITY if(!(d->mode & Y_BIT_SET)) ErrMsg(ER_VARNOTSET,"y coordinate not set")
#define CHECK_BITZ if(!(d->mode & Z_BIT_SET)) ErrMsg(ER_VARNOTSET,"z coordinate not set")

	for (j = 0; j < d->n_X; j++) {
		if (d->colX[j] < -1) {
			switch(d->colX[j]) {
 				case POLY_X: case POLY_X2: case POLY_X3: CHECK_BITX; break;
 				case POLY_Y: case POLY_Y2: case POLY_Y3: CHECK_BITY; break;
 				case POLY_Z: case POLY_Z2: case POLY_Z3: CHECK_BITZ; break;
 				case POLY_XY: CHECK_BITX; CHECK_BITY; break; 
 				case POLY_XZ: CHECK_BITX; CHECK_BITZ; break;
 				case POLY_YZ: CHECK_BITY; CHECK_BITZ; break;
 				case POLY_X2Y: CHECK_BITX; CHECK_BITY; break;
 				case POLY_XY2: CHECK_BITX; CHECK_BITY; break;
 				case POLY_X2Z: CHECK_BITX; CHECK_BITZ; break;
 				case POLY_XZ2: CHECK_BITX; CHECK_BITZ; break;
 				case POLY_Y2Z: CHECK_BITY; CHECK_BITZ; break;
 				case POLY_YZ2: CHECK_BITY; CHECK_BITZ; break;
 				default: ErrMsg(ER_IMPOSVAL, "unknown polynomial number"); break;
 			}
 		}
 	}
	for (j = do_block = 0; !do_block && j < d->n_X; j++)
		do_block = (d->colX[j] < -1);
	for (i = 0; do_block && i < d->n_list; i++)
		/* bl is a single point-list if IS_POINT(d->list[i]) */
		calc_polynomial_point(d, d->list[i]);
}

void calc_polynomial_point(DATA *d, DPOINT *pt) {
	static DATA *bl = NULL;
	int j, k;

	bl = block_discr(bl, get_block_p(), pt);
	for (j = 0; j < d->n_X; j++) {
		if (d->colX[j] < -1)
			/* do eventual block averaging here: */
			for (k = 0, pt->X[j] = 0.0; k < bl->n_list; k++)
				pt->X[j] += bl->list[k]->u.weight *
						calc_polynomial(bl->list[k], d->colX[j]);
	}
}

double calc_polynomial(DPOINT *p, int colX) {
/*
 * fills polynomial field (x, y, z, x2, y2, z2, xy, xz, yz) 
 * with standardized values
 *
 * Counting on the following behaviour:
 * first, all data + valdata are passed through setup_data_minmax(), in 
 * order to get the right values into min and max;
 * then, the routines calc_polynomial* are called.
 * changing min or max inbetween would result in rubbish.
 */
	double x, y, z;

	if (fix_minmax == 0)
		fix_minmax = 1; /* stop touching it */

	x = ((min.x==max.x) ? p->x : (p->x - min.x)/(max.x - min.x));
	y = ((min.y==max.y) ? p->y : (p->y - min.y)/(max.y - min.y));
	z = ((min.z==max.z) ? p->z : (p->z - min.z)/(max.z - min.z));
	switch(colX) {
 		case POLY_X:   return (x); 
 		case POLY_X2:  return (x * x);
 		case POLY_X3:  return (x * x * x);
 		case POLY_Y:   return (y);
 		case POLY_Y2:  return (y * y);
 		case POLY_Y3:  return (y * y * y);
 		case POLY_Z:   return (z);
 		case POLY_Z2:  return (z * z);
 		case POLY_Z3:  return (z * z * z);
 		case POLY_XY:  return (x * y);
 		case POLY_XZ:  return (x * z);
 		case POLY_YZ:  return (y * z);
 		case POLY_X2Y: return (x * x * y);
 		case POLY_XY2: return (x * y * y);
 		case POLY_X2Z: return (x * x * z);
 		case POLY_XZ2: return (x * z * z);
 		case POLY_Y2Z: return (y * y * z);
 		case POLY_YZ2: return (y * z * z);
 		default:      
 			ErrMsg(ER_IMPOSVAL, "unknown polynomial number");
 			break;
 	}
 	return 1.0; /* will never happen */
}

static int average_duplicates(DATA *d) { 
/*
 * average duplicate records (equal coordinates) 
 * and correct number of records;
 * new variance = (sum variances)/(n * n)
 */
	int n_tot = 0, n_dupl, i, j;
	double dzero2;

	if (! d->average)
		return 0;
	dzero2 = gl_zero * gl_zero;
	for (i = 0; i < d->n_list; i++) {
		n_dupl = 0;
		j = i + 1;
		while (j < d->n_list) { /* can't use a for() loop here */
			if (d->pp_norm2(d->list[i], d->list[j]) < dzero2) {
				d->list[i]->attr += d->list[j]->attr;
				d->list[i]->variance += d->list[j]->variance;
				pop_point(d, j); /* decrements d->n_list */
				n_dupl++;
			} else 
				j++;
		} /* while j */
		if (n_dupl > 0) {
			n_tot += n_dupl;
			d->list[i]->attr /= (n_dupl + 1.0);
			d->list[i]->variance /= SQR(n_dupl + 1.0);
		}
	}
	return d->n_averaged = n_tot;
}

#ifndef USING_R
static int read_eas_header(FILE *infile, DATA *d, int ncols) {
	char *line = NULL, *cp;
	int i, size = 0; 

	for (i = 1; i <= ncols; i++) {
		if (get_line(&line, &size, infile) == NULL) {
			pr_warning("GeoEAS file %s: error in header line %d",
				d->fname, 2 + i);
			ErrMsg(ER_READ, d->fname);;
		}
		if ((cp = strtok(line, "\r\n")) != NULL) /* EJP, 27/01/01 */
			; /* do nothing: just want the side effect of strtok() on line */
		if (i == d->colnvalue)
			d->variable = string_dup(line);
		if (i == d->colnx)
			d->x_coord = string_dup(line);
		if (i == d->colny)
			d->y_coord = string_dup(line);
		if (i == d->colnz)
			d->z_coord = string_dup(line);
		if (i == d->colnvariance)
			d->V_coord = string_dup(line);
		if (i == d->colns)
			d->s_coord = string_dup(line);
        if (i == d->coln_id)
            d->id_name = string_dup(line);
	}
	if (line != NULL)
		efree(line);
	return ncols;
}
#endif

static void field_error(char *fname, int line_nr, int fld, char *text) {
	message("gstat:%s:%d: read error on field %d\n", fname, line_nr, fld);
	ErrMsg(ER_RDFLT, text);
}

void free_data(DATA *d) {
	int i;

	assert(d);
	if (DEBUG_FORCE) /* let atexit(qtree_print) do it's job... */
		return;
	if (d->P_base) { /* segmented: */
		efree(d->P_base);
		if (d->n_X && d->X_base)
			efree(d->X_base);
	} else { /* non-segmented */
		if (d->list) /* CW at all MV on output both P_base and d_list are 0 */
			for (i = d->n_list - 1; i >= 0; i--)
				pop_point(d, i);
	}

	if (d->sel != NULL && d->sel != d->list)
		efree(d->sel);
	if (d->list) 
		efree(d->list);
	if (d->colX) 
		efree(d->colX);

	if (d->qtree_root != NULL)
		qtree_free(d->qtree_root);

	if (d->lm) 
		free_lm(d->lm);
	if (d->glm) 
		free_glm(d->glm);

	if (d->grid)
		free_data_gridmap(d->grid);

    if (d->point_ids)
        for (i = d->n_list - 1; i >= 0; i--)
            efree(d->point_ids[i]);
#ifdef HAVE_EXT_DBASE
	if (d->ext_dbase)
	  unlink_ext_dbase(d);
#endif
    
	efree(d);
	return;
}

static int is_one_integer(char *line, int *n) {
	char *dup, *cp;
	int ret_value = 1;

	*n = 0;
	dup = string_dup(line);
	if ((cp = strtok(dup, " \r\t\n")) == NULL) /* # tokens is 0 */
		ret_value = 0;
	if (ret_value && (strtok(NULL, " \r\t\n") != NULL)) /* # tokens > 1 */
		ret_value = 0;
	if (ret_value)
		ret_value = (read_int(cp, n) == 0); /* error */
	efree(dup);
	return ret_value;
}

DATA *init_one_data(DATA *data) {

	if (data == NULL)
		data = (DATA *) emalloc(sizeof(DATA));

	data->colnvalue = 0;
	data->colnx = 0;
	data->colny = 0;
	data->colnz = 0;
	data->colns = 0;
	data->coln_id = 0;
	data->colnvariance = 0;
	data->n_list = -1;
	data->n_max = -1;
	data->nsim_at_data = 0;
	data->init_max = 0;
	data->n_sel = -1;
	data->n_sel_max = 0;
	data->id = -1;
	data->log = 0;
	data->standard = 0;
	data->what_is_u = U_UNKNOWN;
	data->centre = 0;
	data->region = 0;
	data->mode = 0;
	data->dummy = 0;
	data->force = 0;
	data->vdist = 0;
	data->square = 0;
	data->average = 0;
	data->every = 1;
	data->offset = 0;
	data->skip = 0;
	data->prob = 1.0;
	data->lambda = 1.0;
	data->calc_residuals = 1;
	data->is_residual = 0;
	data->polynomial_degree = 0;
	data->n_averaged = 0;
	data->fname = NULL;
	data->type = data_types[DATA_UNKNOWN];
	data->variable = NULL;
	data->x_coord = NULL;
	data->y_coord = NULL;
	data->z_coord = NULL;
	data->s_coord = NULL;
	data->V_coord = NULL;
	data->point_ids = NULL;
	data->id_name = NULL;
	data->mean = 0.0;
	data->std = 0.0;
	data->sel_rad = DBL_MAX;
	data->dX = DBL_MAX;
	data->sel_max = INT_MAX;
	data->sel_min = 0;
	data->oct_max = 0; /* default: don't use octant search */
	data->oct_filled = 1; /* in case of no octants */
	data->list = NULL;
	data->sel = NULL;
	data->P_base = NULL;
	data->X_base = NULL;
	data->lm = NULL;
	data->glm = NULL;
	data->n_merge = 0;
	data->mtbl = NULL;
	data->n_X = 0; 
	data->colX = NULL;
	data_add_X(data, 0); 
	/* add intercept -->> UK, defaulting to ordinary kriging */
	data->qtree_root = NULL;
	data->grid = NULL;
	data->togrid = 0;
	data->point_norm = NULL;
	data->pp_norm2 = NULL;
	data->pb_norm2 = NULL;
	data->beta = NULL;
	data->Category = NULL;
	data->var_fn_str = NULL;
	data->nscore_table = NULL;
	data->variance_fn = NULL;

	set_mv_double(&(data->Icutoff));
	set_mv_double(&(data->mv));

#ifdef HAVE_EXT_DBASE
	data->ext_dbase = NULL;
#endif

	return data;
}

void print_data(const DATA *d, int list) {
	int i;

	printlog("\ndata id: %d\n", d->id);
	if (! is_mv_double(&(d->Icutoff))) 
		printlog("ind. cutoff: %g\n", d->Icutoff);
	if (d->Category) 
		printlog("category: %s\n", d->Category);
	if (! is_mv_double(&(d->mv))) 
		printlog("missing value: %g\n", d->mv);
	if (d->beta) {
		printlog("beta: [");
		for (i = 0; i < d->beta->size; i++)
			printlog(" %g", d->beta->val[i]);
		printlog("]\n");
	}
	printlog("sel_radius %g sel_max %d sel_min %d\n",
		d->sel_rad, d->sel_max, d->sel_min);
	if (d->n_X > 0) {
		for (i = 0; i < d->n_X; i++) {
			printlog("X[%d]: ", i);
			if (d->colX[i] == 0)
				printlog("intercept ");
			if (d->colX[i] < 0) 
				printlog("%s ", POLY_NAME(d->colX[i]));
			if (d->colX[i] > 0) 
				printlog("%d ", d->colX[i]);
		}
		printlog("\n");
	}
	printlog("n_list %d n_max %d n_sel %d\n", d->n_list, d->n_max, d->n_sel); 
	if (list) {
		printlog("current list:\n");
		logprint_data_header(d);
		if (d->n_list) {
#ifdef HAVE_EXT_DBASE
         if (d->type.type == DATA_EXT_DBASE)
			printlog("<extdbase>\n");
		 else
#endif

			for (i = 0; i < d->n_list; i++)
				logprint_point(d->list[i], d);
		} else
			printlog("<empty>\n");
	} else {
		printlog("current selection:\n");
		logprint_data_header(d);
		if (d->n_sel) {
			for (i = 0; i < d->n_sel; i++)
				logprint_point(d->sel[i], d);
		} else 
			printlog("<empty>\n");
	}
}

static void logprint_data_header(const DATA *d) {
	printlog("\nidx x:%s;", save_string(d->x_coord));
	printlog("y:%s;", save_string(d->y_coord));
	printlog("z:%s;", save_string(d->z_coord));
	printlog("v:%s;\n", save_string(d->variable));
}

void logprint_point(const DPOINT *p, const DATA *d) {
/*
 * print contents of p (non zero: use d) to log device
 */
	int j;

	printlog("%3d ", GET_INDEX(p));
	if (d->mode & X_BIT_SET) printlog("x: %4g ", p->x);
	if (d->mode & Y_BIT_SET) printlog("y: %4g ", p->y);
	if (d->mode & Z_BIT_SET) printlog("z: %4g ", p->z);
	if (d->mode & V_BIT_SET) printlog("v: %4g ", p->attr);
	switch (d->what_is_u) {
		case U_UNKNOWN: break;
		case U_ISDIST: printlog("dist: %4g ", sqrt(p->u.dist2)); break;
		case U_ISWEIGHT: printlog("weight: %4g ", p->u.weight); break;
		case U_ISSTRATUM: printlog("stratum: %d ", p->u.stratum); break;
	}
	for (j = 0; j < d->n_X; j++)
		printlog("X[%d]: %6g ", j, p->X[j]);
    if (d->point_ids)
        printlog("ID: %s ", save_string(d->point_ids[GET_INDEX(p)]));
	printlog("\n");
}

void push_point(DATA *d, const DPOINT *p) {
 	int i;
/*
 * add one point p to the data structure d
 * [counts on the fact that erealloc(NULL,size) calls malloc(size)]
 */

 	if (d->prob < 1.0) {
		if (r_uniform() > d->prob)
			return;
	} else if (d->every > 1) {
		/* EJP: WAS  if ((d->n_list + d->offset) % d->every != 0) */
		if ((d->n_list + d->skip + 1 - d->offset) % d->every != 0) {
			d->skip++;
			return;
		}
	}

	if (d->n_list < 0) {
		message("push_point: n_list < 0: %d (%s)\n", d->n_list, d->fname);
		ErrMsg(ER_NULL, "push_point(): n_list < 0");
	}
	if (d->n_max < 0) {
		message("push_point: n_max < 0: %d (%s)\n", d->n_max, d->fname);
		ErrMsg(ER_NULL, "push_point(): n_max < 0");
	}
	/*
	 * use rather large blocks of memory for points:
	 */
	if (d->n_list == d->n_max) { /* increase memory: */
		/* resize d->n_max: */
		if (d->list == NULL) {
			if (d->init_max > 0)
				d->n_max = d->init_max;
			else
				d->n_max = MAX_DATA;
		} else {
			d->n_max += MAX_DATA; /* or else: d->n_max *= 2; */
			if (d->init_max > 0 && DEBUG_DUMP)
				pr_warning("exceeding nmax, now %d", d->n_max);
		}

		/* resize blocked memory bases P_base and X_base, and list: */
		d->P_base = (DPOINT *) erealloc(d->P_base, d->n_max * sizeof(DPOINT));
		if (d->n_X > 0) {
			if (intercept_only(d)) { /* create a single instance of the X row: */
				if (d->X_base == NULL) { /* first time */
					d->X_base = (double *) emalloc(sizeof(double));
					d->X_base[0] = 1.0;
				}
			} else /* each point needs it's own X row: */
				d->X_base = (double *) erealloc(d->X_base, d->n_max * d->n_X * sizeof(double));
		}
		d->list = (DPOINT **) erealloc(d->list, d->n_max * sizeof(DPOINT *));

		/*
		 * realloc'ing may have moved P_base or X_base, so reset all pointers:
		 */
		for (i = 0; i < d->n_list; i++) {
			d->list[i] = &(d->P_base[i]);
			if (d->n_X) {
				if (intercept_only(d))
					/* d->P_base[i].X = d->X_base; */
					d->list[i]->X = d->X_base;
				else
					/* d->P_base[i].X = &(d->X_base[d->n_X * i]); */
					d->list[i]->X = &(d->X_base[d->n_X * i]);
			} else
				/* d->P_base[i].X = NULL; */
				d->list[i]->X = NULL;
		}
		for (i = d->n_list; i < d->n_max; i++)
			d->list[i] = NULL; /* for savety */
		/* rebuild qtree_root: this is avoided by setting nmax */
		qtree_rebuild(d);
		datagrid_rebuild(d, 0);
	}
	/*
	 * copy information on current point into P_base and X_base arrays:
	 */

#ifdef SLOW
	d->P_base[d->n_list] = *p;
#else
	memcpy(&(d->P_base[d->n_list]), p, sizeof(DPOINT));
#endif

	if (d->n_X > 0 && !intercept_only(d)) { 
#define SLOW 1
#ifdef SLOW
	/* slow... copy X row */
		for (i = 0; i < d->n_X; i++)
			d->X_base[d->n_X * d->n_list + i] = p->X[i];
#else
		memcpy(&(d->X_base[d->n_X * d->n_list]), p->X, d->n_X * sizeof(double));
#endif
	} 

	/*
	 * adjust list and X pointer to this copy:
	 */
	d->list[d->n_list] = &(d->P_base[d->n_list]);

	if (intercept_only(d))
		d->list[d->n_list]->X = d->X_base;
	else
		d->list[d->n_list]->X = &(d->X_base[d->n_X * d->n_list]);

	SET_INDEX(d->list[d->n_list], d->n_list);

	qtree_push_point(d, d->list[d->n_list]);
	grid_push_point(d, d->list[d->n_list], 0);
	/* 
	 * this will be ignored during read_gstat_data(), the tree structure will
	 * be filled only during the first call to qtree_quick_select().
	 * Later on, it does have effect if simulated points are pushed.
	 */
	d->n_list++;
	return;
}

void pop_point(DATA *d, int list_nr)
/* 
 * removes DPOINT list_nr, and makes it point to the last DPOINT
 * also changes d->n_list
 */
{
	if (list_nr >= d->n_list) {
		message("pop_point: list_nr >= n_list: %d %d\n", list_nr, d->n_list);
		ErrMsg(ER_NULL, "pop_point():");
	}
	qtree_pop_point(d->list[list_nr], d);
	if (d->P_base == NULL) {
/* 
 * free this one:
 */
		if (d->n_X > 0 && !(intercept_only(d)))
			efree(d->list[list_nr]->X);
		efree(d->list[list_nr]); 
	}
/* 
 * change the last pointer to this:
 */
	if (list_nr != d->n_list - 1) /* we didn't pop the last: */
		d->list[list_nr] = d->list[d->n_list - 1]; 
	d->list[d->n_list - 1] = NULL;
	d->n_list--;
}

static double point_norm_1D(const DPOINT *p) {
/* calculate norm of vector (p->x) */
	return fabs(p->x);
}

static double point_norm_2D(const DPOINT *p) {
/* calculate norm of vector (p->x, p->y, p->z) */
	return sqrt(p->x * p->x +  p->y * p->y);
}

static double point_norm_3D(const DPOINT *p) {
/* calculate norm of vector (p->x, p->y, p->z) */
	return sqrt(p->x * p->x +  p->y * p->y + p->z * p->z);
}

static double pp_norm_1D(const DPOINT *a, const DPOINT *b) {
/* calculate 2-norm of vector (p->x) */
	double dx;
	dx = a->x - b->x;
	return dx * dx;
}

static double pp_norm_2D(const DPOINT *a, const DPOINT *b) {
/* calculate 2-norm of vector (p->x, p->y) */
	double dx, dy;
	dx = a->x - b->x; dy = a->y - b->y;
	return dx * dx + dy * dy;
}

static double pp_norm_3D(const DPOINT *a, const DPOINT *b) {
/* calculate 2-norm of vector (p->x, p->y, p->z) */
	double dx, dy, dz;
	dx = a->x - b->x; dy = a->y - b->y; dz = a->z - b->z;
	return dx * dx + dy * dy + dz * dz;
}

static double point_norm_gc(const DPOINT *p) {
/* calculate norm of vector (p->x, p->y, p->z) */
	ErrMsg(ER_IMPOSVAL, "long/lat: this function should never be called?");
	return gstat_gcdist(p->x, 0.0, p->y, 0.0);
}

double pp_norm_gc(const DPOINT *a, const DPOINT *b) {
	return gstat_gcdist(a->x, b->x, a->y, b->y); /* dist */
}

static double pp_norm_gc2(const DPOINT *a, const DPOINT *b) {
	return pow(gstat_gcdist(a->x, b->x, a->y, b->y), 2.0); /* squared dist */
}

static double pb_norm_gc2(const DPOINT *where, BBOX bbox) {
	/* ErrMsg(ER_IMPOSVAL, "great circle distances cannot be combined with quadtree"); */
	return 0.0; /* always inside, no quadtree */
}


int coordinates_are_equal(const DATA *a, const DATA *b) {
	int i, equal = 1 /* try to disprove equality */;

	if (a->n_list != b->n_list)
		return 0;
	i = 0;
	while (equal && i < a->n_list) {
		equal = ((a->list[i]->x == b->list[i]->x) &&
			(a->list[i]->y == b->list[i]->y) &&
			(a->list[i]->z == b->list[i]->z));
		i++;
	}
	return equal;
}

#ifndef USING_R
static int read_data_from_map(DATA *d) {
	DPOINT current;
	unsigned int i, j, first_cell = 1;
	GRIDMAP *m = NULL;

	current.z = 0.0;
	current.bitfield = 0;
    
	if (d->id == ID_OF_VALDATA && max_block_dimension(0) > 0.0)
		SET_BLOCK(&current);
	else
		SET_POINT(&current);
	current.X = (double *) emalloc(sizeof(double));
	current.X[0] = 1.0;
	m = new_map(READ_ONLY);
	m->filename = d->fname;
	if (map_read(m) == NULL) {
		map_free(m);
		return 0;
	}
	d->grid = gsetup_gridmap(m->x_ul, m->y_ul, m->cellsizex, m->cellsizey, 
			m->rows, m->cols);
	switch(m->type) {
		case MT_UNKNOWN: break;
		case MT_CSF: d->type = data_types[DATA_CSF]; break;
		case MT_ARCGRID: d->type = (m->is_binary ? data_types[DATA_GRIDFLOAT] :
			data_types[DATA_GRIDASCII]); break;
		case MT_IDRISI: 
			d->type = (m->is_binary ? data_types[DATA_IDRISI_BIN] :
				data_types[DATA_IDRISI_ASCII]); break;
		case MT_IDRISI32: 
			d->type = (m->is_binary ? data_types[DATA_IDRISI32_BIN] :
				data_types[DATA_IDRISI32_ASCII]); break;
		case MT_GNUPLOT: d->type = data_types[DATA_GNUPLOT]; break;
		case MT_T2: d->type = data_types[DATA_T2]; break;
		case MT_ERMAPPER: d->type = data_types[DATA_ERMAPPER]; break;
		case MT_GMT: d->type = data_types[DATA_GMT]; break;
		case MT_SURFER: d->type = data_types[DATA_SURFER_DSAA]; break;
		case MT_GSLIB: d->type = data_types[DATA_GSLIB]; break;
		case MT_GRASS: d->type = data_types[DATA_GRASS_GRID]; break;
		case MT_GDAL: d->type = data_types[DATA_GDAL]; break;
		default:
			ErrMsg(ER_IMPOSVAL, "m->type value not evaluated");
	}
	d->n_list = d->n_max = 0;
	d->mode = X_BIT_SET | Y_BIT_SET | V_BIT_SET;
	setup_polynomial_X(d);
	for (i = 0; i < m->rows; i++) {
		for (j = 0; j < m->cols; j++) {
			if (! map_cell_is_mv(m, i, j)) {
				map_rowcol2xy(m, i, j, &(current.x), &(current.y));
				current.attr = map_get_cell(m, i, j);
				push_point(d, &current);
				/* does also d->grid->dpt[i][j] = d->list[d->n_list-1]; */
				if (first_cell) {
					first_cell = 0;
					d->maxX = d->minX = current.x;
					d->maxY = d->minY = current.y;
				} else {
					d->maxX = MAX(d->maxX, current.x);
					d->maxY = MAX(d->maxY, current.y);
					d->minX = MIN(d->minX, current.x);
					d->minY = MIN(d->minY, current.y);
				}
			}
		}
	}
	d->colnx = d->colny = d->colnvalue = d->colnz = 0;
	d->x_coord = string_dup("x (map)");
	d->y_coord = string_dup("y (map)");
	if (m->description)
		d->variable = string_dup(m->description);
	else
		d->variable = string_dup(d->fname);
	map_free(m);
	return 1;
}

static int read_idrisi_points(DATA *d) {
/*
 * read idrisi point file
 */
 	DUMP(d->fname); DUMP(": trying idrisi point format... ");
	if (read_idrisi_point_header(d, string_cat(d->fname, ".dvc")) ||
			read_idrisi_point_data(d, string_cat(d->fname, ".vec"))) {
 		DUMP("no\n");
		return 0;
	}
 	DUMP("yes\n");
 	d->type = data_types[DATA_IDRISI_VEC];
	return 1;
}

static int read_idrisi_point_header(DATA *d, const char *fname) {
	FILE *f;
	char *tok1 = NULL, *tok2 = NULL, *line_buf = NULL;
	int ok = 1, line_size = 0;

	if ((f = fopen(fname, "r")) == NULL)
		return 1;
	while (ok && get_line(&line_buf, &line_size, f) != NULL) {
		tok1 = line_buf; /* first word */
		if (strlen(line_buf) <= 13) {
			pr_warning("line `%s'", line_buf);
			ErrMsg(ER_READ, fname);
		}
		tok2 = line_buf + 14; /* second word */
		tok2 = strtok(tok2, "\n\r");
		line_buf[13] = '\0'; /* cut words */
		if (string_casecmp(tok1, "file title  :") == 0)
			d->variable = string_dup(tok2);
		else if (string_casecmp(tok1, "id type     :") == 0) {
			if (!strstr(tok2, "real")) {
				pr_warning("%f: `id type' should be `real'", fname);
				ok = 0;
			}
		} else if (string_casecmp(tok1, "file type   :") == 0) {
			if (!strstr(tok2, "ascii")) {
				pr_warning("%f: `file type' should be `ascii'", fname);
				ok = 0;
			} 
		} else if (string_casecmp(tok1, "object type :") == 0) {
			if (!strstr(tok2, "point")) {
				pr_warning("%f: `object type' should be `point'", fname);
				ok = 0;
			} 
		} else if (string_casecmp(tok1, "min. X      :") == 0) {
			if (read_double(tok2, &(d->minX)))
				ok = 0;
		} else if (string_casecmp(tok1, "max. X      :") == 0) {
			if (read_double(tok2, &(d->maxX)))
				ok = 0;
		} else if (string_casecmp(tok1, "min. Y      :") == 0) {
			if (read_double(tok2, &(d->minY)))
				ok = 0;
		} else if (string_casecmp(tok1, "max. Y      :") == 0) {
			if (read_double(tok2, &(d->maxY)))
				ok = 0;
		}/* else ignore */
	}
	efclose(f);
	if (!ok) {
		pr_warning("error idrisii point header in %s (error: `%s %s')",
			fname, tok1, tok2);
		return 1;
	} 
	return 0;
}

static int read_idrisi_point_data(DATA *d, const char *fname) {
	FILE *f;
	DPOINT current;
	int n_points = 1, line_size = 0, line_is_xy = 0, ok = 1, nr = 0;
	char *tok1, *tok2, *line = NULL;

	if ((f = fopen(fname, "r")) == NULL)
		return 1;
    current.z = 0.0;
	if (d->id == ID_OF_VALDATA && max_block_dimension(0) > 0.0)
		SET_BLOCK(&current);
	else
		SET_POINT(&current);
	current.X = (double *) emalloc(sizeof(double));
	current.X[0] = 1.0;
	d->n_list = d->n_max = 0;
	d->n_X = 1;  /* just in case... */
	if (d->n_X != 1)
		ErrMsg(ER_IMPOSVAL, "cannot read X columns from idrisi point files");
	d->colX[0] = 0;
	while (get_line(&line, &line_size, f) && n_points == 1 && ok) {
		nr++;
		tok1 = strtok(line, " ");
		tok2 = strtok(NULL, " \n");
		if (line_is_xy) {
			if (read_double(tok1, &(current.x)) || 
					read_double(tok2, &(current.y))) {
				ok = 0;
				pr_warning("%s: read error on line %u: %s %s\n", 
					fname, nr, tok1, tok2);
			}
			push_point(d, &current);
			line_is_xy = 0;
		} else {
			if (read_double(tok1, &(current.attr)) || read_int(tok2, &n_points)
					|| n_points > 1) {
				ok = 0;
				pr_warning("%s: read error on line %u: %s %s\n", fname, nr, tok1, tok2);
			}
			line_is_xy = 1;
		}
	}
	if (! ok)
		ErrMsg(ER_READ, fname);
	d->colnx = d->colny = d->colnvalue = d->colnz = 0;
	d->x_coord = string_dup("x (idr.pt)");
	d->y_coord = string_dup("y (idr.pt)");
	d->mode = X_BIT_SET | Y_BIT_SET | V_BIT_SET;
	efree(line);
	return 0;
}


static int read_idrisi32_points(DATA *d) {
/*
 * read idrisi point file
 */
 	DUMP(d->fname); DUMP(": trying idrisi32 point format... ");
	if (read_idrisi32_point_header(d, string_cat(d->fname, ".vdc")) ||
			read_idrisi32_point_data(d, string_cat(d->fname, ".vct"))) {
 		DUMP("no\n");
		return 0;
	}
 	DUMP("yes\n");
 	d->type = data_types[DATA_IDRISI32_VEC];
	return 1;
}

static int read_idrisi32_point_header(DATA *d, const char *fname) {
	FILE *f;
	char *tok1 = NULL, *tok2 = NULL, *line_buf = NULL;
	unsigned int ok = 1; /* line_size = 0; KS changed to signed */
    int line_size = 0;

    if ((f = fopen(fname, "r")) == NULL)
		return 1;
	while (ok && get_line(&line_buf, &line_size, f) != NULL) {
		tok1 = line_buf; /* first word */
		if (strlen(line_buf) <= 13) {
            pr_warning("line `%s'", line_buf);
			ErrMsg(ER_READ, fname);
		}
		tok2 = line_buf + 14; /* second word */
		tok2 = strtok(tok2, "\n\r");
        if (!tok2) tok2 = "unknown";
		line_buf[13] = '\0'; /* cut words */
		if (string_casecmp(tok1, "file title  :") == 0)
			d->variable = string_dup(tok2);
/*KS*//*	else if (string_casecmp(tok1, "id type     :") == 0) {
			if (!strstr(tok2, "real")) {
				pr_warning("%f: `id type' should be `real'", fname);
				ok = 0;
			}
		} else if (string_casecmp(tok1, "file type   :") == 0) {
			if (!strstr(tok2, "ascii")) {
				pr_warning("%f: `file type' should be `ascii'", fname);
				ok = 0;
			}*/ /*KS replaced */

/*KS*/   else if (string_casecmp(tok1, "id type     :") == 0) {
			if (!strstr(tok2, "integer")) {
                d->datatype = 1;
			}
		} else if (string_casecmp(tok1, "file type   :") == 0) {
			if (!strstr(tok2, "ascii")) {
				d->filetype = 1; /*binary*/
			} /*KS replacement for above 1/18/99*/

		} else if (string_casecmp(tok1, "object type :") == 0) {
			if (!strstr(tok2, "point")) {
                pr_warning("%f: `object type' should be `point'", fname);
				ok = 0;
			}
		} else if (string_casecmp(tok1, "min. X      :") == 0) {
			if (read_double(tok2, &(d->minX)))
				ok = 0;
        } else if (string_casecmp(tok1, "min.X       :") == 0) {
          	if (read_double(tok2, &(d->minX)))
				ok = 0;
		} else if (string_casecmp(tok1, "max. X      :") == 0) {
			if (read_double(tok2, &(d->maxX)))
				ok = 0;
        } else if (string_casecmp(tok1, "max.X       :") == 0) {
          	if (read_double(tok2, &(d->maxX)))
				ok = 0;
		} else if (string_casecmp(tok1, "min. Y      :") == 0) {
			if (read_double(tok2, &(d->minY)))
				ok = 0;
        } else if (string_casecmp(tok1, "min.Y       :") == 0) {
          	if (read_double(tok2, &(d->minY)))
				ok = 0;
		} else if (string_casecmp(tok1, "max. Y      :") == 0) {
			if (read_double(tok2, &(d->maxY)))
				ok = 0;
        } else if (string_casecmp(tok1, "max.Y       :") == 0) {
          	if (read_double(tok2, &(d->maxY)))
				ok = 0;
		}/* else ignore */
	}
	efclose(f);
	if (!ok) {
        pr_warning("error idrisii point header in %s (error: `%s %s')",
			fname, tok1, tok2);
		return 1;
	}
	return 0;
}

static int read_idrisi32_point_data(DATA *d, const char *fname) {
	FILE *f;
	DPOINT current;
	unsigned int /*n_points = 1, line_size = 0,*/ line_is_xy = 0, ok = 1, nr = 0;
    int n_points = 1, line_size =0; /* KS changed from unsigned */
	char *tok1, *tok2, *line = NULL;
    int total, check; double *vecptr; 
	int i; 
	long *n_points2; char *tmp;
		/*KS added*/ /*1/18/99*/

/*KS*/
	if (d->filetype == 0) {    /*KS added*/ /*1/18/99*/
		if (NULL == (f = fopen(fname, "r")))
			return 1;
/*KS*/
	} else {
		/* KS: int fb; if ((fb = open(fname, O_RDONLY|O_BINARY)) == -1) */
		if ((f = fopen(fname, "rb")) == NULL)
			return 1;  /*KS added*/
	}
    current.z = 0.0;
    if (d->id == ID_OF_VALDATA && max_block_dimension(0) > 0.0)
		SET_BLOCK(&current);
	else
		SET_POINT(&current);
	current.X = (double *) emalloc(sizeof(double));
	current.X[0] = 1.0;
	d->n_list = d->n_max = 0;
	d->n_X = 1;  /* just in case... */
	d->colX[0] = 0;

    /*for Idrisi for Windows Version 2.0 ASCII vector files*/
    /*KS*/ 
	if (d->filetype == 0) {    /*KS added*//*1/18/99*/
	while (get_line(&line, &line_size, f) && n_points == 1 && ok) {
		nr++;
		tok1 = strtok(line, " ");
		tok2 = strtok(NULL, " \n");
		if (line_is_xy) {
			if (read_double(tok1, &(current.x)) || read_double(tok2, &(current.y))) {
				ok = 0;
                pr_warning("%s: read error on line %u: %s %s\n", 
					fname, nr, tok1, tok2);
			}
			/* transform_logI(d, &current); */
			push_point(d, &current);
			line_is_xy = 0;
		} else {
			if (read_double(tok1, &(current.attr)) || read_int(tok2, &n_points)
					|| n_points > 1) {
				ok = 0;
                pr_warning("%s: read error on line %u: %s %s\n", 
						fname, nr, tok1, tok2);
			}
			line_is_xy = 1;
		}
	}
	if (! ok)
		ErrMsg(ER_READ, fname);
		efree(line);
    }

   /*for Idrisi for Windows Version 3.0 vector files*/
/*KS*/ 
/*EJP: moved all open()-like fb calls to fopen()-like f calls */
	if (d->filetype == 1) {
    		/* KS: total = filelength(fb); */
			total = file_size(fname);
            total = total - 261;
            if (total < 65520)
    			vecptr = (double *) (emalloc((unsigned int) total));
			else 
				vecptr = (double *) (emalloc(65520));
            tmp = (char *) (emalloc(261));
            n_points2 = (long *) (emalloc(4));
			/* KS:
    		read(fb,tmp,1);
            read(fb,n_points2,4);
            read(fb, tmp, 256);
			*/
			if (fread(tmp, 1, 1, f) != 1 ||
				fread(n_points2, 4, 1, f) != 1 ||
				fread(tmp, 256, 1, f) != 1)
					ErrMsg(ER_READ, fname);

            efree(n_points2); 
			efree(tmp);
			/* EJP: what's all this about? */

            if (total < 65520) {
            	/* KS: read(fb,vecptr,total); */
				if (fread(vecptr, (unsigned int) total, 1, f) != 1)
					ErrMsg(ER_READ, fname);
            	total = (int)(total/8);
            	i = 0;
    			while (i < total) {
              		current.attr = vecptr[i];
              		i++;
              		current.x = vecptr[i];
              		i++;
              		current.y = vecptr[i];
              		i++;
              		/* transform_logI(d, &current); */
			  		push_point(d, &current);
    			}
            }
            else {
            	total = (int)(total/3);
            	/* KS: while (eof(fb) == 0) { */
            	while (feof(f) == 0) {
                	/* KS: check = read(fb,vecptr,65520); */
					check = fread(vecptr,1,65520,f);
                    i = 0;
                    if (check == 65520) {
    					while (i < 8190) {
              				current.attr = vecptr[i];
              				i++;
              				current.x = vecptr[i];
                            i++;
              				current.y = vecptr[i];
              				i++;
              				/* transform_logI(d, &current); */
			  				push_point(d, &current);
    					}
                   } else if (check > 0) {
                   		check = (int)(check/8);
                        while (i < check) {
              				current.attr = vecptr[i];
              				i++;
              				current.x = vecptr[i];
                            i++;
              				current.y = vecptr[i];
              				i++;
			  				push_point(d, &current);
    					}
                   }
               	} /* while */
			}/* else */
         efree(vecptr);
         /* KS: close(fb); */
		 fclose(f);
	}  /*KS added section for binary .vct files*//*1/18/99*/
	d->colnx = d->colny = d->colnvalue = d->colnz = 0;
	d->x_coord = string_dup("x (idr.pt)");
	d->y_coord = string_dup("y (idr.pt)");
	d->mode = X_BIT_SET | Y_BIT_SET | V_BIT_SET;
	return 0;
}

char *print_data_line(const DATA *d, char **to) {
/* 
 * allocates memory for *to, prints the datadescription on *to and returns it
 */
#define ADD_STR(fmt, value) { sprintf(s, fmt, value); strcat(*to, s); }
	char s[100]; /* max length for a single item */
	int est_length = 0, i, start;
	est_length = 10 + 10 * (20 + d->n_X);
	if (d->fname)
		est_length += strlen(d->fname);

	*to = (char *) erealloc(*to, (unsigned int) est_length);
	*to[0] = '\0';
	if (d->fname) {
		strcat(*to, "'");
		strcat(*to, d->fname);
		strcat(*to, "'");
	} else if (d->dummy)
		ADD_STR("%s", "dummy");
	if (d->colnx)
		ADD_STR(", x=%d", d->colnx);
	if (d->colny)
		ADD_STR(", y=%d", d->colny);
	if (d->colnz)
		ADD_STR(", z=%d", d->colnz);
	if (d->colnvalue)
		ADD_STR(", v=%d", d->colnvalue);
    if (d->coln_id)
        ADD_STR(", ID=%d", d->coln_id);
	if (d->colns)
		ADD_STR(", s=%d", d->colns);
	if (d->colnvariance)
		ADD_STR(", V=%d", d->colnvariance);
	if (d->sel_min != 0)
		ADD_STR(", min=%d", d->sel_min); 
	if (d->sel_max != INT_MAX)
		ADD_STR(", max=%d", d->sel_max); 
	if (d->sel_rad < DBL_MAX) 
		ADD_STR(", radius=%g", d->sel_rad); 
	if (d->every != 1)
		ADD_STR(", every=%d", d->every); 
	if (d->offset != 0)
		ADD_STR(", offset=%d", d->offset); 
	if (d->prob != 1.0)
		ADD_STR(", prob=%g", d->prob); 
	if (d->beta) {
		if (intercept_only(d)) {
			ADD_STR(", sk_mean=%g", d->beta->val[0]); 
		} else {
			strcat(*to, ", beta = [");
			for (i = 0; i < d->beta->size; i++) {
				ADD_STR("%g", d->beta->val[i]);
				if (i < d->beta->size - 1)
					strcat(*to, ", ");
			}
			strcat(*to, "]");
		}
	}
	if (d->log) 
		strcat(*to, ", log");
	if (! is_mv_double(&(d->Icutoff)))
		ADD_STR(", I=%g", d->Icutoff); 
	if (d->Category)
		ADD_STR(", Category='%s'", d->Category); 
	if (! is_mv_double(&(d->mv)))
		ADD_STR(", mv=%g", d->mv);
	start = 1;
	if (d->polynomial_degree > 0)
		ADD_STR(", d=%d", d->polynomial_degree);
	if (d->n_X == 0) {
		sprintf(s, ", X=-1");
		strcat(*to, s);
	}
	if (d->n_X > 0 && d->colX[0] != 0) { /* not the default intercept */
		start = 0;
		sprintf(s, ", X=-1");
		strcat(*to, s);
		if (d->colX[0] < 0) { /* -1 is not a polynomial! */
			if (POLY_DEGREE(d->colX[0]) <= d->polynomial_degree)
				sprintf(s, "%s", ""); /* ignore */
			else
				sprintf(s, "&%s", POLY_NAME(d->colX[0]));
		} else
			sprintf(s, "&%d", d->colX[0]);
		strcat(*to, s);
	} 
	if (d->n_X > 1) {
		for (i = 1; i < d->n_X; i++) {
			if (d->colX[i] < 0) {
				/* was this term included by the current pol. degree? */
				if (POLY_DEGREE(d->colX[i]) <= d->polynomial_degree) { /* yes */
					sprintf(s, "%s", ""); /* ignore */
					start++;
				} else
					sprintf(s, "%s%s", i == start ? ", X=" : "&", POLY_NAME(d->colX[i]));
			} else
				sprintf(s, "%s%d", i == start ? ", X=" : "&", d->colX[i]);
			strcat(*to, s);
		}
	}
	if (d->dX < DBL_MAX)
		ADD_STR(", dX=%g", d->dX);
	if (d->square)
		strcat(*to, ", square"); 
	if (d->vdist) 
		strcat(*to, ", vdist");
	if (d->average) 
		strcat(*to, ", average");
	if (! d->calc_residuals) 
		strcat(*to, ", noresiduals");
	strcat(*to, ";"); 
	return *to;
} /* print_data_line */
#endif

int push_to_merge_table(DATA *d, int to_var, int col_this_X, int col_other_X) {
	int i;
	DATA **data;

	data = get_gstat_data();
	if (to_var >= d->id) { /* should not occur by construction */
		pr_warning("use push_to_merge_table only backwards (%d >= %d)",
				to_var, d->id);
		return 1;
	}
	if (col_this_X >= d->n_X || col_other_X >= data[to_var]->n_X) {
		pr_warning("merge error: variable out of range");
		return 1;
	}
	if (d->beta || data[to_var]->beta) {
		pr_warning("cannot merge to or from fixed (known) parameters");
		return 1;
	}
	for (i = 0; i < d->n_merge; i++) {
		if (col_this_X == d->mtbl[i].col_this_X) {
			pr_warning("merge error: cannot merge column twice");
			return 1;
		}
	}
	d->n_merge++;
	d->mtbl = (MERGE_TABLE *) erealloc(d->mtbl,
		d->n_merge * sizeof (MERGE_TABLE));
	d->mtbl[d->n_merge - 1].to_var = to_var;
	d->mtbl[d->n_merge - 1].col_this_X = col_this_X;
	d->mtbl[d->n_merge - 1].col_other_X = col_other_X;
	return 0;
}

DATA_GRIDMAP *gsetup_gridmap(double x_ul, double y_ul, double cellsizex, 
			double cellsizey, unsigned int rows, unsigned int cols) {

	DATA_GRIDMAP *t;
	int i, j;

	t = (DATA_GRIDMAP *) emalloc(sizeof(DATA_GRIDMAP));
	t->x_ul = x_ul;
	t->y_ul = y_ul;
	t->cellsizex = cellsizex;
	t->cellsizey = cellsizey;
	t->rows = rows;
	t->cols = cols;
	t->dpt = (DPOINT ***) emalloc(t->rows * sizeof(DPOINT **));
	t->grid_base = (DPOINT **) emalloc(t->rows * t->cols * sizeof(DPOINT *));
	for (i = 0; i < t->rows; i++)
		t->dpt[i] = &(t->grid_base[i * t->cols]);
	for (i = 0; i < t->rows; i++)
		for (j = 0; j < t->cols; j++)
			t->dpt[i][j] = NULL;
	return t;
}

static void free_data_gridmap(DATA_GRIDMAP *t) {
	efree(t->grid_base);
	efree(t->dpt);
}

static void grid_push_point(DATA *d, DPOINT *p, int adjust_to_gridcentres) {
	int row, col;

	if (d->grid) {
		row = floor((d->grid->y_ul - p->y)/d->grid->cellsizey);
		col = floor((p->x - d->grid->x_ul)/d->grid->cellsizex);
		row = MAX(0, row);
		row = MIN(row, d->grid->rows - 1);
		col = MAX(0, col);
		col = MIN(col, d->grid->cols - 1);
		d->grid->dpt[row][col] = p;
		if (adjust_to_gridcentres) {
			p->x = d->grid->x_ul + (col + 0.5) * d->grid->cellsizex;
			p->y = d->grid->y_ul - (row + 0.5) * d->grid->cellsizey;
		}
	}
	return;
}

void datagrid_rebuild(DATA *d, int adjust_to_gridcentres) {
	int i;

	if (d->grid)
		for (i = 0; i < d->n_list; i++)
			grid_push_point(d, d->list[i], adjust_to_gridcentres);
	return;
}

double data_block_diagonal(DATA *data) {
	DPOINT a, b;

	a.x = data->maxX;
	b.x = data->minX;
	if (data->mode & Y_BIT_SET) {
		a.y = data->maxY;
		b.y = data->minY;
	} else {
		a.y = 0.0;
		b.y = 0.0;
	}
	if (data->mode & Z_BIT_SET) {
		a.z = data->maxZ;
		b.z = data->minZ;
	} else {
		a.z = 0.0;
		b.z = 0.0;
	}
	return sqrt(data->pp_norm2(&a, &b));
}

D_VECTOR *push_d_vector(double d, D_VECTOR *v) {
	if (v == NULL) {
		v = (D_VECTOR *) emalloc(sizeof(D_VECTOR));
		v->size = v->max_size = 0;
		v->val = NULL;
	}
	v->size++;
	if (v->size > v->max_size) { /* (re)allocate v->val */
		if (v->val == NULL)
			v->val = (double *) emalloc(v->size * sizeof(double));
		else
			v->val = (double *) erealloc(v->val, v->size * sizeof(double));
		v->max_size = v->size;
	}
	v->val[v->size - 1] = d;
	return v;
}

void free_d_vector(D_VECTOR *v) {
	if (v != NULL) {
		if (v->size > 0)
			efree(v->val);
		efree(v);
	}
}

int intercept_only(const DATA *d) {
	assert(d != NULL);
	return (d->n_X == 1 && d->colX[0] == 0);
}

double v_mu(double mu) {
	return mu;
}

double v_mu2(double mu) {
	return mu * mu;
}

double v_mu3(double mu) {
	return mu * mu * mu;
}

double v_bin(double mu) {
	return (mu * (1.0 - mu));
}

double v_identity(double mu) {
	return 1.0;
}

#ifndef USING_R
#ifdef HAVE_LIBGIS
static DATA *read_grass_data(DATA * d) {
/* From: "main.c" for GRASS Program "v.bubble". 
  1 feb 20000 : Job Spijker spijker@geo.uu.nl
*/
	char *vect_mapset;
	char *site_mapset;
	int dims = 0, cat = 0, strs = 0, dbls = 0, nsites, ret, ignored, outside;
	FILE *fd;
	DPOINT point;
	struct Cell_head region;
	Site_head info;
	Site *site;

	point.u.stratum = 0;
/* Make sure that the current projection is UTM or defined-99 or  */
/* unreferenced XY projection.                       */
	if (!gl_longlat) { /* user apparently knows what she/he does! */
		if ((G_projection() != 0) && (G_projection() != 1)
				&& (G_projection() != 99))
			G_fatal_error(
			"%s:  Projection must be either Unreferenced XY (value 0) or \
	UTM (value 1).  Change the value \"proj\" in file \"WIND\" to either \
	0 or 1 and then re-execute \"%s\".",
				G_program_name(), G_program_name());
	}

	/* Obtain the mapset name for the chosen site_lists file d->fname */
	/*
	 * EJP, 01/18/05; grass 6.0beta1
	if ((site_mapset = G_find_file("site_lists", d->fname, "")) == NULL)
	*/
	if ((site_mapset = G_find_sites(d->fname, "")) == NULL)
		/*
		G_fatal_error("%s: Site_list file: <%s> does not exist.",
				G_program_name(), d->fname);
		*/
		return NULL;

	/* Get "mapset". */
	vect_mapset = G_mapset();

	/* Open site_lists file. */
	if ((fd = G_fopen_sites_old(d->fname, site_mapset))== NULL)
		/*
		G_fatal_error("%s: Unable to open site_lists file <%s> in mapset <%s>",
				G_program_name(), d->fname, site_mapset);
		*/
		return NULL;

	if (G_site_describe(fd, &dims, &cat, &strs, &dbls) != 0)
		G_fatal_error("%s: G_site_describe() failed to guess format!",
				G_program_name());

	/* Provide warning message if projection is Unreferenced XY. */
	if (G_projection() == 0)
		pr_warning("%s: projection is 0 (is this a problem?)", 
			G_program_name());

	G_site_get_head(fd, &info);
	G_site_put_head(stderr, &info);

	fprintf(stderr, "GRASS site list %s: %d cat, %d dim, %d str, %d dbl.\n", 
			d->fname, cat, dims, strs, dbls);
	site = G_site_new_struct(cat, dims, strs, dbls);
	/* cat is CELL_TYPE, FCELL_TYPE, DCELL_TYPE, or -1 (no category) */

	if (site->dbl_alloc > 0 && d->colnvalue > site->dbl_alloc)
		pr_warning("gstat/grass doubles in site list are ignored, check data v field!");

	/* Get the region/window */
	G_get_window(&region);
	
	if (d->colnvalue < 1) /* Big Revolution: gstat takes a sensible default!! */
		d->colnvalue = 1;

	d->n_list = d->n_max = 0; 
	d->minX = region.west;
	d->maxX = region.east;
	d->minY = region.south;
	d->maxY = region.north;

	d->mode = X_BIT_SET | Y_BIT_SET | V_BIT_SET;
	if (d->colnz > dims)
		pr_warning(
			"gstat/grass: ignoring z coordinate (%d) larger than dims (%d)", 
			d->colnz, dims);
	else if (d->colnz == 1 || d->colnz == 2)
		pr_warning(
			"gstat/grass: ignoring z coordinate (should not be 1 or 2)");
	else if (dims > 2) { /* deal with z? */
		if (d->colnz <= 0)
			pr_warning("gstat/grass: specify z coordinate if you want to use it\n");
		else {
			printlog("using coordinate column %d as z-coordinate", d->colnz);
			/* d->mode = d->mode & Z_BIT_SET; */
			d->mode = d->mode | Z_BIT_SET;
		}
	}
	d->colnx = 1; 
	d->colny = 2;
	mk_var_names(d);

	/* do something with bsite here */
	point.X = (double *) emalloc(sizeof(double));
	point.X[0] = 1;

	nsites = ret = ignored = outside = 0;
	while ((ret = G_site_get(fd, site)) != -1) {
		nsites++;
		/* printf("[%s]\n", site->str_att[0]); */
		switch (ret) {
			case -2: /* fatal error or insufficient data */
				pr_warning(
					"fatal error/insufficient data (-2) on site %d (ignored)", 
					nsites);
				ignored++;
				break;
			case 1: /* format mismatch (extra data) */
				pr_warning(
					"format mismatch (extra data) (-1) on site %d (ignored)", 
					nsites);
				ignored++;
				break;
			case 0: /* success */
				if (d->region == 0 || G_site_in_region(site, &region)) {
					point.x = site->east;
					point.y = site->north;
					if (d->colnz) {
						point.z = site->dim[d->colnz - 3];
						printf("dim:[%g]\n", point.z);
						if (d->n_list == 0) {
							d->minZ = point.z;
							d->maxZ = point.z;
						} else {
							d->minZ = MIN(d->minZ, point.z);
							d->maxZ = MAX(d->maxZ, point.z);
						}
					}
					/* now try to get something sensible here... */
					if (d->colnvalue <= site->dbl_alloc) {
						point.attr = site->dbl_att[d->colnvalue - 1];
					} else if (d->colnvalue <= site->str_alloc) {
						if (d->Category)
							point.attr = strcmp(d->Category, 
								site->str_att[d->colnvalue - 1]) ? 0 : 1;
						else if (read_double(site->str_att[d->colnvalue - 1],
								&(point.attr)))
							ErrMsg(ER_RDFLT, site->str_att[d->colnvalue - 1]);
					} else { 
						if (d->colnvalue > 1)
							ErrMsg(ER_IMPOSVAL, 
								"data v field not found in sites list");
						switch (site->cattype) {
							case CELL_TYPE:
								point.attr = site->ccat;
								break;
							case FCELL_TYPE:
								point.attr = site->fcat;
								break;
							case DCELL_TYPE:
								point.attr = site->dcat;
								break;
							case -1:
							default:
								if (d->id != ID_OF_VALDATA)
									ErrMsg(ER_IMPOSVAL, 
							"no attribute information available from sites file");
								break;
						}
					}
					push_point(d, &point);
				} else 
					outside++;
				break;
			default:
				ErrMsg(ER_IMPOSVAL, "G_site_get(): unknown return value");
				break;
		} /* switch ret */
	} /* while */

	if (ignored)
		printlog("gstat/grass: %d sites ignored (bad format).\n", ignored);
	if (outside)
		printlog("gstat/grass: ignored %d sites outside region.\n", outside);
	printlog("gstat/grass: %d sites read successfully.\n", d->n_list);
	d->type = data_types[DATA_GRASS];
	G_site_free_struct(site);
	/* EJP, 01/18/05, grass 6.0beta1
	fclose(fd);
	*/
	G_sites_close(fd);
	return d;
}
#endif
#endif
