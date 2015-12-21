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

/*! \file mapio.c 
	\brief generic grid map I/O library functions 

mapio.c offers a generic grid (raster) map i/o interface to gstat.
It offers functions for reading and writing files, and access to cell
values (get/put). Raster map format is auto-detected.

Current implementation supports the following formats:
	arcinfo gridascii
	arcinfo gridfloat
	idrisi ascii
	idrisi binary real
	PCRaster/csf (all formats in, REAL4 out)
	ER-Mapper (all formats; 4-byte real output)
	GMT (netcdf)
	T2 (mike-she grid maps) (when #ifdef HAVE_T2_GRIDFORMAT)
	Surfer (DSAA/ascii)

TODO:
	? implement optional arcinfo/idrisi map read/writing on a per cel basis
	- integrate map topology information in local searches, esp. for cs.
	[[[ for simulation, do use a storage on a grid basis, not in a
	DATA structure: efficient because of grid ordering -> easy searching]]]

	GSLIB grid definition (is it useful? read as a DATA structure with a tag?):
	xmn ymn zmn # coordinates centre of first block
	nx ny nz    # number of blocks;
	xsiz ysiz zsiz # block size in x,y and z direction
	value1      # values following, x cycles fastest, then y, then z
	value2      # loc = (iz-1)nx.ny + (iy-1)nx + ix
	...	 # iz = 1 + int(loc/(nx.ny))
		    # iy = 1 + int((loc-(iz-1)nx.ny)/nx)
		    # ix = loc - (iz-1)nx.ny-(iy-1)nx
*/

/* 
 * if (gl_rowwise != 0), no complete map is not kept in memory.
 * take care that the complete map should be blanked (RputAllMV())
 * before writing takes place, when applying this method to other
 * formats (e.g. GDAL), as map_set_row() is only applied to WRITE_ONLY
 * maps for rows that contain non-missing valued cells!
 * SECOND: WRITE_ONLY maps should be open in read-write (r+) state!!
 * */

#include <stdio.h>
#include <math.h>				/* floor() */
#include <float.h>				/* FLT_MAX */
#include <string.h>				/* strtok() */
#include <ctype.h>				/* isalpha() */
#include <stdlib.h>				/* exit */

#include "defs.h"

#ifdef HAVE_LIBCSF
#include "csf.h"
#endif

#ifdef HAVE_LIBGIS
# include "gis.h"
#endif

/*! if defined, a stand-alone mapio library can be compiled */
#ifndef MAPIO_LIB
# include "glvars.h"
#else							/* standalone version: */
# include "defaults.h"
int debug_level = 1, gl_secure = 0, gl_rowwise = 1;
double gl_zero = DEF_zero;
char *gl_mv_string = "NA";
#endif

#include "utils.h"
#include "debug.h"
#include "read.h"
#include "userio.h"
#include "mapio.h"

/* functions, specific to a package */
static GRIDMAP *read_arcgrid(GRIDMAP * m);
static GRIDMAP *write_arcgrid(GRIDMAP * m);
static GRIDMAP *read_idrisi_image(GRIDMAP * m);
static GRIDMAP *write_idrisi(GRIDMAP * m);
static GRIDMAP *read_idrisi32_image(GRIDMAP * m);
static GRIDMAP *write_idrisi32(GRIDMAP * m);
static GRIDMAP *write_gnuplot_binary(GRIDMAP * m);
static GRIDMAP *read_ermapper(GRIDMAP * m);
static GRIDMAP *write_ermapper(GRIDMAP * m);
static GRIDMAP *read_surfer(GRIDMAP * m);
static GRIDMAP *write_surfer(GRIDMAP * m);
#ifdef WITH_GSLIB
static GRIDMAP *read_gslib(GRIDMAP * m);
#endif
static GRIDMAP *write_gslib(GRIDMAP * m);
static GRIDMAP *write_error(GRIDMAP * m);

#ifdef HAVE_T2_GRIDFORMAT
static GRIDMAP *read_T2(GRIDMAP * m);
static GRIDMAP *write_T2(GRIDMAP * m);
#endif

#ifdef HAVE_LIBCSF
static GRIDMAP *read_csf(GRIDMAP * m);
static GRIDMAP *write_csf(GRIDMAP * m);
static GRIDMAP *dup_csf(GRIDMAP * m, GRIDMAP * dup);
void CsfReadRow(GRIDMAP *m, float *buf, unsigned int row);
void CsfWriteRow(GRIDMAP *m, float *buf, unsigned int row);
#endif

#ifdef HAVE_LIBNETCDF
static GRIDMAP *read_gmt(GRIDMAP * m);
static GRIDMAP *write_gmt(GRIDMAP * m);
static int cdf_read_grd_info(GRIDMAP * m, double *x_max, double *y_min,
							 double *z_scale_factor, double *z_add_offset,
							 int *node_offset);
#endif

#ifdef HAVE_LIBGDAL
#include "gdal.h"
#include "cpl_error.h"
#include "cpl_string.h"
static GRIDMAP *read_gdal(GRIDMAP * m);
static GRIDMAP *write_gdal(GRIDMAP * m);
#endif

#ifdef HAVE_LIBGIS
static GRIDMAP *read_grass(GRIDMAP * m);
static GRIDMAP *write_grass(GRIDMAP * m);
#endif

/* generic functions: */
static int read_arcgrid_header(GRIDMAP * m, FILE * f);
static int write_arcgrid_header(GRIDMAP * m, FILE * f);
static int read_idrisi_header(GRIDMAP * m, const char *fname);
static int write_idrisi_header(GRIDMAP * m, const char *fname);
static int read_idrisi32_header(GRIDMAP * m, const char *fname);
static int write_idrisi32_header(GRIDMAP * m, const char *fname);
static int read_ermapper_header(GRIDMAP * m, FILE * f);
static void write_ermapper_header(GRIDMAP * m, FILE * f);
static int read_surfer_header(GRIDMAP * m, FILE * f);
static int read_ascii_grid(GRIDMAP * m, FILE * f,
						   int first_line_in_buffer);
static void write_ascii_grid(GRIDMAP * m, FILE * f, int as_rows);
static int read_multiformat_bgrid(GRIDMAP * m, const char *fname,
								  int swap);
static void write_binary_grid(GRIDMAP * m, const char *fname, int swap);
static void swap_floats(unsigned char *b, unsigned int n);
static void swap_multiformat(unsigned char *b, unsigned int m,
							 unsigned int n);
static unsigned int sizeof_ct(CellType ct);
static void alloc_mv_grid(GRIDMAP * m);
static void map_set_row(GRIDMAP *m, unsigned int row);

static float default_misval_idrisi = 0.0;

#define SWAP_N(a,n) swap_floats((unsigned char *)a,n)
#define SWAP_M_N(a,m,n) swap_multiformat((unsigned char *)a,m,n)


#define CHECK_ROWS     1
#define CHECK_COLS     2
#define CHECK_CELLSIZE 4
#define CHECK_X_UL     8
#define CHECK_Y_UL    16
#define CHECK_SUM     31		/* sum of all checks */

#define BINARY_NATIVE      1
#define BINARY_NON_NATIVE  2
#define DEFAULT_MISVAL -9999.0
#define SURFER_MISVAL 1.70141E+38

static char *line_buf = NULL;
static int line_size = 0;

static char *er_datum = NULL;
static char *er_projection = NULL;

/*
 * create a new GRIDMAP structure
 * allocates memory and initializes all fields for a GRIDMAP structure
 * returns: pointer to GRIDMAP structure
 */
GRIDMAP *new_map(MAP_READ_STATUS status)
{
	GRIDMAP *map;

	map = (GRIDMAP *) emalloc(sizeof(GRIDMAP));
	map->status = status;
	map->type = MT_UNKNOWN;
	map->history = NULL;
	map->description = NULL;
	map->filename = NULL;
	map->rows = 0;
	map->cols = 0;
	map->base_size = 0;
	map->grid = NULL;
	map->base = NULL;
	map->first_time_row = NULL;
	map->is_binary = 0;
	map->celltype = CT_UNKNOWN;
	map->misval = DEFAULT_MISVAL;	/* only for arcgrid */
	map->cellmin = map->cellmax = FLT_MAX;
	map->CSF_MAP = NULL;
#ifdef HAVE_LIBGDAL
	map->GeoTransform = (double *) emalloc(6 * sizeof(double));
#endif
	map->write = write_error;
	map->read_row = map->write_row = NULL;
	map->current_row = 0;
	return map;
}

static void map_set_row(GRIDMAP *m, unsigned int new_row)
{
	int start = 0, i;

	if (m->status == READ_ONLY && m->read_row == NULL)
		return;
	if (m->status == WRITE_ONLY && m->write_row == NULL)
		return;

	if (m->grid == NULL) {
		m->grid = (float **) emalloc(m->rows * sizeof(float *));
		m->first_time_row = (unsigned int *) emalloc(m->rows * 
					sizeof(unsigned int));
		m->current_row = 0;
		m->grid[0] = (float *) emalloc(m->cols * sizeof(float));
		memset(m->grid[0], 0xFF, m->cols * sizeof(float));
		for (i = 1; i < m->rows; i++) {
			m->grid[i] = NULL;
			m->first_time_row[i] = 1;
		}
		start = 1;
	}

	assert(m->grid[m->current_row] != NULL);
	assert(new_row >= 0 && new_row < m->rows);

	switch (m->status) {
		case READ_ONLY:
			if (m->current_row != new_row || start)
				m->read_row(m, m->grid[m->current_row], new_row);
			break;
		case WRITE_ONLY:
			if (new_row != m->current_row) { /* flush old row to disk: */
				m->write_row(m, m->grid[m->current_row], m->current_row);
				if (m->first_time_row[new_row]) {
					memset(m->grid[m->current_row], 0xFF, 
							m->cols * sizeof(float));
					m->first_time_row[new_row] = 0;
				} else { /* been here before: re-read the written content: */
					m->read_row(m, m->grid[m->current_row], new_row);
					/* adjust this pointer to [new_row] a few lines down: */
				}
			}
			break;
		default:
			ErrMsg(ER_IMPOSVAL, "unknown switch");
			break;
	}
	/* now adjust buffer pointer; */
	if (new_row != m->current_row) { 
		m->grid[new_row] = m->grid[m->current_row];
		m->grid[m->current_row] = NULL;
		m->current_row = new_row;
	}
	assert(m->grid[m->current_row] != NULL);
}

/*
 * read in a grid map and fill all necessary fields;
 * at least m->filename should be filled before returns: 
 * a pointer to a GRIDMAP read, NULL in case the map is not
 * one of the formats supported.
 */
#ifndef USING_R
GRIDMAP *map_read(GRIDMAP * m)
{

	assert(m);
	assert(m->filename);
	assert(m->status == READ_ONLY);

#ifdef HAVE_LIBGIS /* try grass: */
	if (read_grass(m))
		return m;
	DUMP("no\n")
#endif

#ifdef HAVE_LIBGDAL /* try gdal: */
	if (read_gdal(m))
		return m;
	DUMP("no\n")
#endif

#ifdef HAVE_LIBCSF /* try csf: */
		if (read_csf(m))
		return m;
	DUMP("no\n")
#endif

#ifdef HAVE_LIBNETCDF /* try GMT grid: */
		if (read_gmt(m))
		return m;
	DUMP("no\n")
#endif

	/* try ERMAPPER dataset: */
	if (read_ermapper(m))
		return m;
	DUMP("no\n")

	/* try IDRISI: */
	if (read_idrisi_image(m))
		return m;
	DUMP("no\n")

	/* try IDRISI32: */
	if (read_idrisi32_image(m))
		return m;
	DUMP("no\n")

	/* try ARCGRID */
	if (read_arcgrid(m))
		return m;
	DUMP("no\n");

	/* try SURFER DSAA */
	if (read_surfer(m))
		return m;
	DUMP("no\n");

#ifdef TRY_GSLIB /* try GSLIB grid */
	if (read_gslib(m))
		return m;
	DUMP("no\n");
#endif

#ifdef HAVE_T2_GRIDFORMAT
	/* try T2 */
	if (read_T2(m))
		return m;
	DUMP("no\n");
#endif

	/* no success: */
	return NULL;
}
#endif

static GRIDMAP *write_error(GRIDMAP * m)
{
	pr_warning("%s: writing this map format is not supported", m->filename);
	assert(0);
	return NULL;
}

static void alloc_mv_grid(GRIDMAP * m)
{
	unsigned int i;

	m->base_size = m->rows * m->cols;
	m->base = (float *) emalloc(m->base_size * sizeof(float));
	memset(m->base, 0xFF, m->base_size * sizeof(float));

	m->grid = (float **) emalloc(m->rows * sizeof(float *));
	for (i = 0; i < m->rows; i++) /* set up handles: */
		m->grid[i] = &(m->base[i * m->cols]);
}

/*
 * give x,y coordinate of cell center for cell [row, col]
 * libcsf has it's own function; other formats assume increasing x
 * for increasing cols and decreasing y for increasing rows
 * returns: non-zero if row or col are outside map limits
 */
int map_rowcol2xy(GRIDMAP * m,	/* pointer to gridmap */
				  unsigned int row,	/* current row number */
				  unsigned int col,	/* current column number */
				  double *x,	/* return value: pointer to x-coordinate */
				  double *y /* return value: pointer to y-coordinate */ )
{
	assert(m);
	assert(x);
	assert(y);

	if (row >= m->rows || col >= m->cols)
		return 1;
#if defined(HAVE_LIBCSF)
	if (m->CSF_MAP)
		RgetCoords((MAP *) m->CSF_MAP, 1, row, col, x, y);
	else
#endif
	{
		*x = m->x_ul + (col + 0.5) * m->cellsizex;
		*y = m->y_ul - (row + 0.5) * m->cellsizey;
	}
	return 0;
}

/*
 * converts x and y coordinate to (row,col) pair.
 *
 * see comment for map_rowcol2xy()
 *
 * returns: non-zero if x or y are outside map limits
 */
int map_xy2rowcol(GRIDMAP * m /* pointer to map */ ,
				  double x,		/* x-coordinate */
				  double y,		/* y-coordinate */
				  unsigned int *row,	/* output value: pointer to row number */
				  unsigned int *col
				  /* output value: pointer to column number */ )
{
	assert(m);
	assert(row);
	assert(col);

#if defined(HAVE_LIBCSF)
	if (m->CSF_MAP) {			/* handle possible non-zero map angle, CSF 2+: */
		if (RgetRowCol
			((MAP *) m->CSF_MAP, x, y, (size_t *) row,
			 (size_t *) col) != 1) return 1;
	} else
#endif
	{
		if (x < m->x_ul || x > m->x_ul + m->cols * m->cellsizex ||
			y > m->y_ul || y < m->y_ul - m->rows * m->cellsizey)
			return 1;
		*row = (unsigned int) floor((m->y_ul - y) / m->cellsizey);
		*col = (unsigned int) floor((x - m->x_ul) / m->cellsizex);
		if (*row == m->rows) /* on the bottom edge */
			*row = *row - 1;
		if (*col == m->cols) /* on the right edge */
			*col = *col - 1;
	}
	return 0;
}

/*
 * check whether a cell contains a missing value
 * returns: non-zero if cell is missing value
 */
int map_cell_is_mv(GRIDMAP * m /* pointer to map */ ,
				   unsigned int row /* row number */ ,
				   unsigned int col /* column number */ )
{
	assert(m);
	assert(row < m->rows && col < m->cols);
	map_set_row(m, row);

	return is_mv_float(&(m->grid[row][col]));
}

/*
 * return cell value of a given map cell
 * this function assumes that the given cell is not a missing value cell
 * returns: cell value
 */
float map_get_cell(GRIDMAP * m /* pointer to GRIDMAP */ ,
				   unsigned int row,	/* row number */
				   unsigned int col /* column number */ )
{
	assert(m);
	assert(row < m->rows && col < m->cols);
	map_set_row(m, row);
	assert(!is_mv_float(&(m->grid[row][col])));

	return m->grid[row][col];
}

/* 
 * write cell value to map
 * keeps track of map minimum and maximum values
 * returns: 0
 */
int map_put_cell(GRIDMAP * m,	/* pointer to GRIDMAP structure */
		 unsigned int row,	/* row number */
		 unsigned int col,	/* column number */
		 float value /* cell value to write */ )
{
	assert(m);
	assert(row < m->rows && col < m->cols);

	map_set_row(m, row);

	if (m->grid == NULL)
		alloc_mv_grid(m);
	/* assert(m->grid); */
	assert(m->grid[row]);
	m->grid[row][col] = value;
	if (m->cellmin == FLT_MAX)
		m->cellmin = m->cellmax = value;
	else {
		m->cellmin = MIN(m->cellmin, value);
		m->cellmax = MAX(m->cellmax, value);
	}
	return 0;
}

/*
 * duplicates map 
 * makes a ``shallow'' copy, no map contents are copied
 *
 * returns: the duplicate map
 */
GRIDMAP *map_dup(const char *fname, GRIDMAP *m)
{
	GRIDMAP *dup = NULL;
#ifdef HAVE_LIBGDAL
	char **papszOptions = NULL, *drv = NULL;
#endif

	assert(m);
	assert(fname);

	dup = new_map(WRITE_ONLY);
	*dup = *m;					/* copy all fields; incl status */
#ifdef HAVE_LIBGDAL
	/*
	if (strcmp(GDALGetDriverShortName(m->hDriver), "PCRaster") == 0) {
		papszOptions = CSLSetNameValue( papszOptions, "PCRASTER_VALUESCALE", "VS_SCALAR" );
		papszOptions = CSLSetNameValue( papszOptions, "PCRASTER_CELLREPRESENTATION", "CR_REAL4" );
	}
	*/
	CPLPushErrorHandler(CPLQuietErrorHandler);
	dup->hDataset = GDALCreateCopy(m->hDriver, fname, m->hDataset, FALSE, papszOptions, NULL, NULL);
	if (dup->hDataset == NULL) {
		CPLPopErrorHandler();
		CPLErrorReset();
		pr_warning("GDAL cannot CreateCopy this format; using asciigrid for %s", fname);
		dup->is_binary = 0;
		dup->write = write_arcgrid;
	} else 
		CPLPopErrorHandler();
	/* if (GDALGetRasterAccess(dup->hDataset) != GA_Update)
		printf("NO ACCESS 0\n");
	*/

#endif
	dup->cellmin = dup->cellmax = FLT_MAX;	/* re-initialize */
	dup->CSF_MAP = NULL;		/* copying would only cause trouble! */
	dup->filename = string_dup(fname);
	dup->status = WRITE_ONLY;
	dup->current_row = 0;
	if (dup->write_row == NULL) /* do the whole block in memory */
		alloc_mv_grid(dup);
	else
		dup->grid = NULL;
#ifdef HAVE_LIBCSF
	if (dup->type == MT_CSF && (dup = dup_csf(m, dup)) == NULL)
		ErrMsg(ER_WRITE, fname);
#endif
	return dup;
}

/* 
 * free map structure and contents
 */
void map_free(GRIDMAP * m /* pointer to GRIDMAP structure */ )
{

	assert(m);

	if (m->grid)
		efree(m->grid);
	if (m->base)
		efree(m->base);

#ifdef HAVE_LIBCSF
	if (m->type == MT_CSF && m->CSF_MAP) {
		if (Mclose((MAP *) m->CSF_MAP)) {
			Mperror(m->filename);
			ErrMsg(ER_WRITE, m->filename);
		}
		m->CSF_MAP = NULL;
	}
#endif
	efree(m);
}

#ifndef USING_R
/*
 * check whether two maps have equal topologies 
 * compares cellsizes, location and dimension 
 * returns: 1 if maps are equal, else 0
 */
int map_equal(GRIDMAP * a, GRIDMAP * b)
{
	assert(a);
	assert(b);

#ifdef HAVE_LIBCSF
	if (a->CSF_MAP && b->CSF_MAP)
		return Rcompare((MAP *) a->CSF_MAP, (MAP *) b->CSF_MAP);
	/* Rcompare() checks also for angle, projection, ... */
#endif
	if (a->cellsizex != b->cellsizex ||
		a->cellsizey != b->cellsizey ||
		a->x_ul != b->x_ul ||
		a->y_ul != b->y_ul || a->rows != b->rows || a->cols != b->cols)
		return 0;
	return 1;
}

/* 
 * read map as arcgrid map
 * first checks for grid map format resulting from exporting a grid
 * map in ArcGrid with the gridfloat command (filename and filename.hdr
 * should both exist in this case); then checks for mapformat resulting
 * from exporting with the gridascii command (ascii grid with several
 * header lines).
 * returns: pointer to map structure on success, else NULL 
 */
static GRIDMAP *read_arcgrid(GRIDMAP * m)
{
	FILE *f;
	int flt = 0;
	const char *flt_name, *hdr_name;

	assert(m);

	DUMP(m->filename);
	DUMP(": trying ArcInfo format... ");

	hdr_name = string_cat(m->filename, ".hdr");
	flt_name = string_cat(m->filename, ".flt");
	flt = file_exists(string_cat(m->filename, ".flt"));

	if (file_exists(hdr_name) && (flt || file_exists(m->filename))) {
/* floatgrid/gridfloat: file.hdr and <file.flt or file>*/
		m->is_binary = 1;
		f = efopen(string_cat(m->filename, ".hdr"), "r");
		if (read_arcgrid_header(m, f)) {
			efclose(f);
			return NULL;
		}
		m->celltype = CT_IEEE4;
		if (read_multiformat_bgrid(m, flt ? flt_name : m->filename,
								   m->is_binary == BINARY_NON_NATIVE)) {
			m->celltype = CT_UNKNOWN;
			return NULL;
		}
		DUMP("floatgrid/gridfloat\n");
	} else if (file_exists(m->filename)) {	/* gridascii: one single file */
		m->is_binary = 0;
		f = efopen(m->filename, "r");
		if (read_arcgrid_header(m, f) || read_ascii_grid(m, f, 1)) {
			efclose(f);
			return NULL;
		}
		DUMP("asciigrid/gridascii\n");
	} else						/* m->filename and m->filename.hdr both don't exist: */
		return NULL;
	m->type = MT_ARCGRID;
	m->write = write_arcgrid;
	efclose(f);
	return m;
}

static GRIDMAP *write_arcgrid(GRIDMAP * m)
{
	FILE *f;

	if (!m->is_binary) {
		f = efopen(m->filename, "w");
		write_arcgrid_header(m, f);
		write_ascii_grid(m, f, 1);
		efclose(f);
	} else {
		f = efopen(string_cat(m->filename, ".hdr"), "w");
		write_arcgrid_header(m, f);
		efclose(f);
		write_binary_grid(m, string_cat(m->filename, ".flt"), 0);
	}
	return m;
}

static int read_arcgrid_header(GRIDMAP * m, FILE * f)
{
	char *tok1 = NULL, *tok2 = NULL;
	unsigned int check = 0, ok = 1, centerx = 0, centery = 0;
	int i;

	while (ok && get_line(&line_buf, &line_size, f) != NULL &&
		   isalpha(*line_buf)) {
		tok1 = strtok(line_buf, " \t\n\r");	/* first word */
		tok2 = strtok(NULL, " \t\n\r");	/* second word */
		if (tok1 == NULL || tok2 == NULL) {
			ok = 0;
		} else if (string_casecmp(tok1, "ncols") == 0) {
			if (read_int(tok2, &i))
				ok = 0;
			m->cols = i;
			check = check | CHECK_COLS;
		} else if (string_casecmp(tok1, "nrows") == 0) {
			if (read_int(tok2, &i))
				ok = 0;
			m->rows = i;
			check = check | CHECK_ROWS;
		} else if (string_casecmp(tok1, "xllcenter") == 0) {
			if (read_double(tok2, &(m->x_ul)))
				ok = 0;
			centerx = 1;
			check = check | CHECK_X_UL;
		} else if (string_casecmp(tok1, "yllcenter") == 0) {
			if (read_double(tok2, &(m->y_ul)))
				ok = 0;
			centery = 1;
			check = check | CHECK_Y_UL;
		} else if (string_casecmp(tok1, "xllcorner") == 0) {
			if (read_double(tok2, &(m->x_ul)))
				ok = 0;
			check = check | CHECK_X_UL;
		} else if (string_casecmp(tok1, "yllcorner") == 0) {
			if (read_double(tok2, &(m->y_ul)))
				ok = 0;
			check = check | CHECK_Y_UL;
		} else if (string_casecmp(tok1, "cellsize") == 0) {
			if (read_double(tok2, &(m->cellsizex)))
				ok = 0;
			check = check | CHECK_CELLSIZE;
		} else if (string_casecmp(tok1, "nodata_value") == 0) {
			if (read_float(tok2, &(m->misval)))
				ok = 0;
		} else if (string_casecmp(tok1, "byteorder") == 0) {
			if (string_casecmp(tok2, "msbfirst") == 0) {
				if (cpu_is_little_endian())
					m->is_binary = BINARY_NON_NATIVE;
				else
					m->is_binary = BINARY_NATIVE;
			} else if (string_casecmp(tok2, "lsbfirst") == 0) {
				if (cpu_is_little_endian())
					m->is_binary = BINARY_NATIVE;
				else
					m->is_binary = BINARY_NON_NATIVE;
			} else {
				pr_warning("byteorder `%s' not recognized", tok2);
				ErrMsg(ER_IMPOSVAL,
					   "only lsbfirst and msbfirst are supported");
			}
		} else {
			if (check)
				pr_warning("unrecognized word in arcgrid file `%s': %s",
						   m->filename, tok1);
			return 1;
		}
	}
	if (check != CHECK_SUM) {
		if (check)
			pr_warning("uncomplete arcgrid header in %s", m->filename);
		return 1;
	}
	m->cellsizey = m->cellsizex;
	if (centery)
		m->y_ul += (m->rows - 0.5) * m->cellsizey;
	else
		m->y_ul += m->rows * m->cellsizey;
	if (centerx)
		m->x_ul -= 0.5 * m->cellsizex;
	return 0;
}

static int write_arcgrid_header(GRIDMAP * m, FILE * f)
{
	fprintf(f, "NCOLS     %12d\n", m->cols);
	fprintf(f, "NROWS     %12d\n", m->rows);
	fprintf(f, "XLLCORNER %12.12g\n", m->x_ul);
	fprintf(f, "YLLCORNER %12.12g\n", m->y_ul - m->cellsizey * m->rows);
	fprintf(f, "CELLSIZE  %12.12g\n", SQUARECELLSIZE(m));
	/* if (m->misval != DEFAULT_MISVAL) */
	/* fprintf(f, "NODATA_VALUE %9.9g\n", m->misval); */
	fprintf(f, "NODATA_VALUE %g\n", m->misval);
	if (m->is_binary)
		fprintf(f, "BYTEORDER %s\n", cpu_is_little_endian()?
				"LSBFIRST" : "MSBFIRST");
	return 0;
}

static GRIDMAP *read_idrisi_image(GRIDMAP * m)
{
	FILE *f = NULL;

	DUMP(m->filename);
	DUMP(": trying Idrisi image format... ");
	if (read_idrisi_header(m, string_cat(m->filename, ".doc")))
		return NULL;
	if (m->is_binary) {
		if (read_multiformat_bgrid(m, string_cat(m->filename, ".img"),
								   !cpu_is_little_endian()))
			return NULL;
	} else {
		if ((f = fopen(string_cat(m->filename, ".img"), "r")) == NULL)
			return NULL;
		if (read_ascii_grid(m, f, 0)) {
			efclose(f);
			return NULL;
		}
		efclose(f);
	}
	m->type = MT_IDRISI;
	m->write = write_idrisi;
	DUMP("yes\n");
	return m;
}

static GRIDMAP *write_idrisi(GRIDMAP * m)
{
	FILE *f = NULL;

	write_idrisi_header(m, string_cat(m->filename, ".doc"));
	if (m->is_binary)
		write_binary_grid(m, string_cat(m->filename, ".img"),
						  !cpu_is_little_endian());
	else {
		f = efopen(string_cat(m->filename, ".img"), "w");
		write_ascii_grid(m, f, 0);
		efclose(f);
	}
	return m;
}

static int read_idrisi_header(GRIDMAP * m, const char *fname)
{
	FILE *f;
	char *tok1 = NULL, *tok2 = NULL;
	unsigned int check = 0, ok = 1;
	enum { IS_UNKNOWN, IS_REAL, IS_BYTE } is_what = IS_UNKNOWN;
	int i;

	if ((f = fopen(fname, "r")) == NULL)
		return 1;
	while (ok && get_line(&line_buf, &line_size, f) != NULL) {
		tok1 = line_buf;		/* first word */
		if (strlen(line_buf) <= 13) {
			pr_warning("line `%s'", line_buf);
			ErrMsg(ER_READ, fname);
		}
		tok2 = line_buf + 14;	/* second word */
		tok2 = strtok(tok2, "\n\r");
		line_buf[13] = '\0';	/* cut words */
		if (DEBUG_DUMP)
			printlog("[%s][%s]\n", tok1, tok2);
		if (string_casecmp(tok1, "data type   :") == 0) {
			if (strstr(tok2, "real"))
				is_what = IS_REAL;
			else if (strstr(tok2, "byte"))
				is_what = IS_BYTE;
			else if (!strstr(tok2, "integer"))
				ok = 0;
		} else if (string_casecmp(tok1, "file title  :") == 0) {
			if (tok2 != NULL)
				m->description = string_dup(tok2);
		} else if (string_casecmp(tok1, "file type   :") == 0) {
			if (strstr(tok2, "ascii"))
				m->is_binary = 0;
			else if (strstr(tok2, "binary")) {
				if (is_what == IS_UNKNOWN)
					ErrMsg(ER_IMPOSVAL,
						   "can only read binary files of data type `real' or `byte'");
				m->is_binary = (cpu_is_little_endian()? BINARY_NATIVE :
								BINARY_NON_NATIVE);
				if (is_what == IS_BYTE)
					m->celltype = CT_UINT8;
				else
					m->celltype = CT_IEEE4;
			} else
				ok = 0;
		} else if (string_casecmp(tok1, "columns     :") == 0) {
			if (read_int(tok2, &i))
				ok = 0;
			m->cols = i;
			check = check | CHECK_COLS;
		} else if (string_casecmp(tok1, "rows        :") == 0) {
			if (read_int(tok2, &i))
				ok = 0;
			m->rows = i;
			check = check | CHECK_ROWS;
		} else if (string_casecmp(tok1, "min. X      :") == 0) {
			if (read_double(tok2, &(m->x_ul)))
				ok = 0;
			check = check | CHECK_X_UL;
		} else if (string_casecmp(tok1, "max. Y      :") == 0) {
			if (read_double(tok2, &(m->y_ul)))
				ok = 0;
			check = check | CHECK_Y_UL;
		} else if (string_casecmp(tok1, "resolution  :") == 0) {
			if (read_double(tok2, &(m->cellsizex)))
				ok = 0;
			check = check | CHECK_CELLSIZE;
		} else if (string_casecmp(tok1, "min. value  :") == 0) {
			if (read_float(tok2, &(m->cellmin)))
				ok = 0;
		} else if (string_casecmp(tok1, "max. value  :") == 0) {
			if (read_float(tok2, &(m->cellmax)))
				ok = 0;
		} else if (string_casecmp(tok1, "flag value  :") == 0) {
			if (read_float(tok2, &(m->misval)))	/* no numeric value: */
				m->misval = DEFAULT_MISVAL;	/* re-initialize */
		}						/* else ignore */
	}
	fclose(f);					/* DO NOT use efclose */
	m->cellsizey = m->cellsizex;
	if (!ok || check != CHECK_SUM) {
		if (check)
			pr_warning
				("uncomplete idrisii header in %s (error: `%s %s'), %d",
				 m->filename, tok1, tok2, check);
		return 1;
	}
	/* check: */
	return 0;
}

static int write_idrisi_header(GRIDMAP * m, const char *fname) {
	FILE *f;
	char *cp;
	int start = 1;

	f = efopen(fname, "w");
	fprintf(f, "file title  : %s",
			m->description ? m->description : "none\n");
	if (m->celltype == CT_UINT8 && m->is_binary)
		fprintf(f, "data type   : byte\n");
	else
		fprintf(f, "data type   : real\n");
	if (m->is_binary)
		fprintf(f, "file type   : binary\n");
	else
		fprintf(f, "file type   : ascii\n");
	fprintf(f, "columns     : %u\n", m->cols);
	fprintf(f, "rows        : %u\n", m->rows);
	fprintf(f, "ref. system : plane\n");
	fprintf(f, "ref. units  : m\n");
	fprintf(f, "unit dist.  : 1\n");
	fprintf(f, "min. X      : %.7f\n", m->x_ul);
	fprintf(f, "max. X      : %.7f\n", m->x_ul + m->cellsizex * m->cols);
	fprintf(f, "min. Y      : %.7f\n", m->y_ul - m->cellsizey * m->rows);
	fprintf(f, "max. Y      : %.7f\n", m->y_ul);
	fprintf(f, "pos'n error : unknown\n");
	fprintf(f, "resolution  : %.7f\n", m->cellsizex);
	fprintf(f, "min. value  : %.7f\n", m->cellmin);
	fprintf(f, "max. value  : %.7f\n", m->cellmax);
	fprintf(f, "value units : no\n");	/* no classes, change into ??? */
	fprintf(f, "value error : unknown\n");	/* is a map */
	if (m->celltype == CT_UINT8 && m->is_binary)
		fprintf(f, "flag value  : 255\n");
	else
		fprintf(f, "flag value  : %.7f\n", m->misval);
	fprintf(f, "flag def'n  : missing data\n");
	fprintf(f, "legend cats : 0\n");
	if (m->history) {
		while ((cp = strtok(start ? m->history : NULL, "\n")) != NULL) {
			start = 0;
			fprintf(f, "comment     : %s\n", cp);
		}
	}
	efclose(f);
	return 0;
}

static GRIDMAP *read_idrisi32_image(GRIDMAP *m) { 
		/* KS is Kirstin Schneider -- in 1999 she did the Idrisi32 conversion
		 * and released her code under GPL. Her changes to my original code
		 * are marked with KS.  Tue Jul  9 14:48:30 CEST 2002,
		 * EJP did the merge and added ``32'' to all routines
		 * and constants; and tried to restore support for non-Intel platforms */
		/*KS changed Idrisi for Windows Version 3.0 file formats*/
		/*.doc -> .rdc     .img -> .rst  */
	FILE *f = NULL;                               

	DUMP("trying idrisi32 image format... ");
	if (read_idrisi32_header(m, string_cat(m->filename, ".rdc")))
		return NULL;
	if (m->is_binary) {
		/* KS: if (read_multiformat_bgrid2(m, string_cat(m->filename, ".rst"), 0)) */
		if (read_multiformat_bgrid(m, string_cat(m->filename, ".rst"), !cpu_is_little_endian()))
			return NULL;
	} else {
		if (NULL == (f = fopen(string_cat(m->filename, ".rst"), "r"))) /*KS switched NULL order*/
			return NULL;
		if (read_ascii_grid(m, f, 0)) {
			efclose(f);
			return NULL;
		}
		efclose(f);
	}
	m->type = MT_IDRISI32; /* EJP -> added 32 */
    m->write = write_idrisi32; /* EJP -> added 32 */
	DUMP("yes\n");
	return m;
}

static GRIDMAP *write_idrisi32(GRIDMAP *m) {    
		/*KS changed Idrisi for Windows Version 3.0 file formats*/
		/*.doc -> .rdc     .img -> .rst  */
	FILE *f = NULL;                             

	if (write_idrisi32_header(m, string_cat(m->filename, ".rdc")))
		return NULL;
	if (m->is_binary) {
		write_binary_grid(m, string_cat(m->filename, ".rst"),
				!cpu_is_little_endian());
			return NULL;
	} else {
		f = efopen(string_cat(m->filename, ".rst"), "w");
		write_ascii_grid(m, f, 0);
        efclose(f);
		return NULL;
	}
	return m;
}

static int read_idrisi32_header(GRIDMAP *m, const char *fname) {
	FILE *f;
	char *tok1 = NULL, *tok2 = NULL;
	unsigned int check = 0, ok = 1 /* , is_real = 0 */ ;
    double unit_dist; /* KS added */
	double ymin, xmax; /* EJP added, removed from GRIDMAP structure, 
						 where KS put them in */
    enum { IS_UNKNOWN, IS_REAL, IS_BYTE } is_what = IS_UNKNOWN;

	if (NULL == (f = fopen(fname, "r")))  /*KS switched NULL*/
		return 1;
	while (ok && NULL != get_line(&line_buf, &line_size, f)) {/*KS switched NULL*/
        tok1 = NULL; /* KS added because NT did not like*/
		tok1 = line_buf; /* first word */
		if (strlen(line_buf) <= 13) {
        	if (string_casecmp(tok1, "\n") != 0) { /*KS added line*/
               pr_warning("line `%s'", line_buf);
			   ErrMsg(ER_READ, fname);
            }
		}
        tok2 = NULL; /*KS added because NT did not like*/
		tok2 = line_buf + 14; /* second word */
		tok2 = strtok(tok2, "\n\r");
		line_buf[13] = '\0'; /* cut words */
		/*
		if (DEBUG_DUMP)
			printf("[%s][%s]\n", tok1, tok2);
		*/
		if (string_casecmp(tok1, "data type   :") == 0) {
			if (!strstr(tok2, "integer") && !strstr(tok2, "real"))
				ok = 0;
			if (strstr(tok2, "real")) {
				/* is_real = 1;  */
				is_what = IS_REAL;
             }
		} else if (string_casecmp(tok1, "file title  :") == 0) {
		/*	m->description = string_dup(tok2); */
		} else if (string_casecmp(tok1, "file type   :") == 0) {
            if (strstr(tok2, "ascii"))
				m->is_binary = 0;
			else if (strstr(tok2, "binary")) {
				switch(is_what) {
					case IS_UNKNOWN:
						ErrMsg(ER_IMPOSVAL,
						"can only read binary files of data type `real' or `byte'");
						break;
					case IS_BYTE:
						m->celltype = CT_UINT8;
						break;
					case IS_REAL:
						m->celltype = CT_IEEE4;
						break;
				}
				/* KS did: m->is_binary = 1; EJP: changed back to: */
				m->is_binary = (cpu_is_little_endian() ? BINARY_NATIVE :
					BINARY_NON_NATIVE);
			} else
				ok = 0;
		} else if (string_casecmp(tok1, "columns     :") == 0) {
			if (read_uint(tok2, &(m->cols))) /* EJP: -> read_uint */
				ok = 0;
			check = check | CHECK_COLS;
		} else if (string_casecmp(tok1, "rows        :") == 0) {
			if (read_uint(tok2, &(m->rows))) /* EJP -> read_uint */
				ok = 0;
			check = check | CHECK_ROWS;
/*KS*/  } else if (string_casecmp(tok1, "unit dist.  :") == 0) {
            if (read_double(tok2, &(unit_dist)))
				ok = 0;  /*KS added unit distance*/
		} else if (string_casecmp(tok1, "min. X      :") == 0) {
			if (read_double(tok2, &(m->x_ul)))
				ok = 0;
			check = check | CHECK_X_UL;
/*KS*/ 	} else if (string_casecmp(tok1, "min. Y      :") == 0) {
			if (read_double(tok2, &(ymin)))
				ok = 0;
        } else if (string_casecmp(tok1, "max. Y      :") == 0) {
			if (read_double(tok2, &(m->y_ul)))
				ok = 0;
			check = check | CHECK_Y_UL;
/*KS*/  } else if (string_casecmp(tok1, "max. X      :") == 0) {
			if (read_double(tok2, &(xmax)))
				ok = 0;  /*KS added wxmax*/
		} else if (string_casecmp(tok1, "resolution  :") == 0) {
/*KS*/    /*		if (read_double(tok2, &(m->cellsize)))
            	ok = 0;   KS removed*/
/*KS*/      m->cellsizex = ((xmax - m->x_ul)/m->cols) * unit_dist;
            m->cellsizey = ((m->y_ul - ymin)/m->rows) * unit_dist; /*KS added*/
            check = check | CHECK_CELLSIZE;
		} else if (string_casecmp(tok1, "min. value  :") == 0) {
			if (read_float(tok2, &(m->cellmin)))
				ok = 0;
		} else if (string_casecmp(tok1, "max. value  :") == 0) {
			if (read_float(tok2, &(m->cellmax)))
				ok = 0;
		} else if (string_casecmp(tok1, "flag value  :") == 0) {
			if (read_float(tok2, &(m->misval))) /* no numeric value: */
			m->misval = default_misval_idrisi; /* re-initialize */
		}/* else ignore */
	}
	fclose(f); /* EJP: don't use efclose */
    /*efree(tok2); efree(tok1);*/
	if (!ok || check != CHECK_SUM) {
		if (check) 
            pr_warning("uncomplete idrisi header in %s (error: `%s %s')",
					m->filename, tok1, tok2);
		return 1;
	}
	/* check: */
    return 0;
}

static int write_idrisi32_header(GRIDMAP *m, const char *fname) {
	FILE *f;
	char *cp;
	int start = 1;

	f = efopen(fname, "w");
    fprintf(f, "file format : IDRISI Raster A.1\n"); /*KS added for new Idrisi rdc format*/
	fprintf(f, "file title  : %s", m->description ? m->description : "none\n");
	fprintf(f, "data type   : real\n");
	if (m->is_binary)
	fprintf(f, "file type   : binary\n");
	else
	fprintf(f, "file type   : ascii\n");
	fprintf(f, "columns     : %u\n", m->cols);
	fprintf(f, "rows        : %u\n", m->rows);
	fprintf(f, "ref. system : plane\n");
	fprintf(f, "ref. units  : m\n");
	fprintf(f, "unit dist.  : 1\n");
	fprintf(f, "min. X      : %.7f\n", m->x_ul);
	fprintf(f, "max. X      : %.7f\n", m->x_ul + m->cellsizex * m->cols);
	fprintf(f, "min. Y      : %.7f\n", m->y_ul - m->cellsizey * m->rows);
	fprintf(f, "max. Y      : %.7f\n", m->y_ul);
	fprintf(f, "pos'n error : unknown\n");
	fprintf(f, "resolution  : %.7f\n", m->cellsizex);
    if (m->misval<m->cellmin) fprintf(f, "min. value  : %.7f\n", m->misval);
		else fprintf(f, "min. value  : %.7f\n", m->cellmin);
	fprintf(f, "max. value  : %.7f\n", m->cellmax);
    fprintf(f, "display min : %.7f\n", m->cellmin); /*KS added for new Idrisi rdc format*/
	fprintf(f, "display max : %.7f\n", m->cellmax); /*KS added for new Idrisi rdc format*/
	fprintf(f, "value units : no\n"); /* no classes, change into ??? */
	fprintf(f, "value error : unknown\n"); /* is a map */
	fprintf(f, "flag value  : %.7f\n", m->misval);
	fprintf(f, "flag def'n  : missing data\n");
	fprintf(f, "legend cats : 0\n");
	while ((cp = strtok(start ? m->history : NULL, "\n")) != NULL) {
		start = 0;
		fprintf(f, "lineage     : %s\n", cp);
	}
	efclose(f);
	return 0;
}

static int read_ascii_grid(GRIDMAP * m, FILE * f, int first_line_in_buffer)
{
	unsigned int cells = 0, ok = 1, row, col, at_start;
	char *tok;

	cells = 0;
	alloc_mv_grid(m);
	if (!first_line_in_buffer)
		get_line(&line_buf, &line_size, f);
	if (line_size <= 0)
		return 1;
	do {
		at_start = 1;			/* we're at a new line */
		while (ok && (tok = strtok(at_start ? line_buf : NULL, " \t\n\r"))
			   != NULL) {
			at_start = 0;
			row = cells / m->cols;
			col = cells % m->cols;
			if (cells >= m->rows * m->cols ||
				read_float(tok, &(m->grid[row][col])))
				ok = 0;
			if (ok && m->grid[row][col] == m->misval)
				set_mv_float(&(m->grid[row][col]));
			else {
				if (m->cellmin == FLT_MAX)
					m->cellmin = m->cellmax = m->grid[row][col];
				else {
					m->cellmin = MIN(m->cellmin, m->grid[row][col]);
					m->cellmax = MAX(m->cellmax, m->grid[row][col]);
				}
			}
			cells++;
		}
	} while (ok && get_line(&line_buf, &line_size, f) != NULL);
	if (!ok || cells != m->rows * m->cols) {
		if (cells > 0)
			pr_warning("`%s' has too %s values (read: %u, expected: %u)\n",
					   m->filename,
					   cells < m->rows * m->cols ? "few" : "many", 
					   cells, m->rows * m->cols);
		return 1;
	}
	return 0;
}

static void write_ascii_grid(GRIDMAP * m, FILE * f, int as_rows)
{
	unsigned int i, j;
	char *fmt;

	fmt = as_rows ? "%s%g" : "%s%g\n";
	for (i = 0; i < m->rows; i++) {
		for (j = 0; j < m->cols; j++)
			fprintf(f, fmt, (j && as_rows) ? " " : "",
					is_mv_float(&(m->grid[i][j])) ? m->misval : m->grid[i][j]);
		if (as_rows)
			fprintf(f, "\n");
	}
}

static void write_binary_grid(GRIDMAP * m, const char *fname, int swap)
{
	int i, j;
	FILE *f;
	unsigned char b;			/* single byte */

	assert(fname);

	f = efopen(fname, "wb");
	for (i = 0; i < m->rows; i++) {
		if (m->celltype == CT_UINT8) {
			for (j = 0; j < m->cols; j++) {	/* write out this row */
				if (map_cell_is_mv(m, i, j))
					b = 0xFF;	/* 255, or what else? */
				else
					b = m->grid[i][j];
				if (fwrite(&b, 1, 1, f) != 1)
					ErrMsg(ER_WRITE, fname);
			}
		} else {				/* default -- write IEEE4 (float) grid */
			for (j = 0; j < m->cols; j++)
				if (map_cell_is_mv(m, i, j))
					m->grid[i][j] = m->misval;	/* replace 0xFFFF's w. misval */
			if (swap)
				SWAP_N(m->grid[i], m->cols);
		}
	}
	if (m->celltype != CT_UINT8) {
		if (fwrite(m->base, sizeof(float), m->base_size, f) !=
			m->base_size) ErrMsg(ER_WRITE, fname);
	}
	efclose(f);
	if (DEBUG_DUMP)
		printf("binary map %s: %d cells\n", fname, m->base_size);
	return;
}

/* 
 * write map in gnuplot binary format
 *
 * This was written to convert maps to the gnuplot binary format
 * in order to display them as a semi 3-D surface (splot) plot.
 * I didn't get it to work because (by the time I tested) gnuplot's
 * splot couldn't handle the missing values well (or the Inf's they
 * are converted to)
 *
 * I never wrote a function for reading gnuplot binary grids 
 * since nobody seems to use them for maps. No surprise--they don't
 * contain any topology information.
 *
 * returns: NULL on error, else pointer to gridmap 
 */
static GRIDMAP *write_gnuplot_binary(GRIDMAP * m)
{
	FILE *bfile = NULL;
	float value;
	double x, y;
	int j, k;
#ifndef HAVE_LIBCSF
	typedef int INT4;
#endif
	const INT4 Inf = 0x7f800000;	/* bits for the float Infinity value */

	if ((bfile = fopen(m->filename, "wb")) == NULL)
		return NULL;
	value = (float) m->cols;
	fwrite(&value, 4, 1, bfile);	/* <ncols> */
	for (j = 0; j < m->cols; j++) {
		map_rowcol2xy(m, 0, j, &x, &y);
		value = (float) x;
		fwrite(&value, 4, 1, bfile);	/* <x0> <x1>,...,<xncols-1> */
	}
	for (j = 0; j < m->rows; j++) {	/* for each row j */
		map_rowcol2xy(m, j, 0, &x, &y);
		value = (float) y;
		fwrite(&value, 4, 1, bfile);	/* coordinate of row */
		for (k = 0; k < m->cols; k++) {	/* followed by row values */
			if (map_cell_is_mv(m, j, k) /* || m->grid[j][k] == m->misval */
				)
				/* set_mv_float(&value); */
				memcpy(&value, &Inf, sizeof(float));
			else
				value = map_get_cell(m, j, k);
			fwrite(&value, 4, 1, bfile);
		}
	}
	fclose(bfile);
	return m;
}

/* swap endian-ness of n floats held in b */
static void swap_floats(unsigned char *b, unsigned int n)
{
	unsigned char tmp;
	int i;

	/* Order: 0123 => 3210 */
	for (i = 0; i < n; i++) {
		tmp = b[0];
		b[0] = b[3];
		b[3] = tmp;				/* swap 0 and 3 */
		tmp = b[1];
		b[1] = b[2];
		b[2] = tmp;				/* swap 1 and 2 */
		b += 4;
	}
}

#ifdef HAVE_LIBCSF
static GRIDMAP *read_csf(GRIDMAP * m)
{
	MAP *m_tmp = NULL;			/* CSF map pointer */

	DUMP(m->filename);
	DUMP(": trying csf map format... ");
	if ((m_tmp = Mopen(m->filename, M_READ)) != NULL) {
		if (MgetProjection(m_tmp) != PT_YDECT2B) {
			pr_warning("for file %s:", m->filename);
			ErrMsg(ER_IMPOSVAL, 
			"projection type should be: y increases from bottom to top");
		}
		m->cellsizex = m->cellsizey = RgetCellSize(m_tmp);
		m->x_ul = RgetXUL(m_tmp);
		m->y_ul = RgetYUL(m_tmp);
		if (RuseAs(m_tmp, CR_REAL4) != 0)
			ErrMsg(ER_IMPOSVAL, "cannot read map as REAL4");
		else
			m->celltype = CT_IEEE4;
		if (MgetMapDataType(m_tmp) != T_RASTER)
			ErrMsg(ER_IMPOSVAL, "cannot read non-raster maps");
		m->rows = RgetNrRows(m_tmp);
		m->cols = RgetNrCols(m_tmp);
		RgetMinVal(m_tmp, &(m->cellmin));
		RgetMaxVal(m_tmp, &(m->cellmax));
		m->type = MT_CSF;

		if (gl_rowwise) {
			m->read_row = CsfReadRow;
			m->write_row = CsfWriteRow;
		} else {
			alloc_mv_grid(m);
			if (RgetSomeCells(m_tmp, 0, m->base_size, m->base) != m->base_size)
				ErrMsg(ER_READ, m->filename);
		}

		m->CSF_MAP = (void *) m_tmp;
		m->write = write_csf;
		DUMP("yes\n");
		return m;
	}
	ResetMerrno();				/* reset CSF error before next call */
	return NULL;
}

static GRIDMAP *write_csf(GRIDMAP * m)
{
	if (m->CSF_MAP == NULL)
		ErrMsg(ER_NULL, "write_csf()");
	if (m->write_row != NULL && m->grid) /* flush last row: */
		m->write_row(m, m->grid[m->current_row], m->current_row);
	else { /* flush complete map: */
		if (RputSomeCells((MAP *) m->CSF_MAP, 0, m->base_size, m->base) != 
						m->base_size)
			ErrMsg(ER_WRITE, m->filename);
	}
	if (m->history)
		MputHistory((MAP *) m->CSF_MAP, m->history);	/* full history */
	if (m->description)
		MputDescription((MAP *) m->CSF_MAP, m->description);

	if (Mclose((MAP *) m->CSF_MAP)) {
		Mperror(m->filename);
		return NULL;
	}
	m->CSF_MAP = NULL;
	return m;
}

void CsfReadRow(GRIDMAP *m, float *buf, unsigned int row)
{
	if (RgetRow(m->CSF_MAP, row, buf) != m->cols)
		ErrMsg(ER_READ, m->filename);
}

void CsfWriteRow(GRIDMAP *m, float *buf, unsigned int row)
{
	if (RputRow(m->CSF_MAP, row, buf) != m->cols)
		ErrMsg(ER_READ, m->filename);
}

static GRIDMAP *dup_csf(GRIDMAP * m, GRIDMAP * dup)
{

	if (m->CSF_MAP != NULL) {	/* use csf duplicate routines: */
		if (dup->celltype == CT_UINT8) {
			if ((dup->CSF_MAP = Rdup(dup->filename, (MAP *) m->CSF_MAP,
					CR_UINT1, VS_NOMINAL)) == NULL) {
				ErrMsg(ER_WRITE, dup->filename);
				return NULL;
			}
		} else {
			if ((dup->CSF_MAP = Rdup(dup->filename, (MAP *) m->CSF_MAP,
					CR_REAL4, VS_SCALAR)) == NULL) {
				ErrMsg(ER_WRITE, dup->filename);
				return NULL;
			}
		}
		if (dup->celltype == CT_UINT8 && RuseAs(dup->CSF_MAP, CR_REAL4)) 
				ErrMsg(ER_IMPOSVAL, "cannot use map as REAL4");
	} else {					/* convert, i.e. create: */
		dup->CSF_MAP = Rcreate(dup->filename, dup->rows, dup->cols, CR_REAL4,
				VS_SCALAR, PT_YDECT2B, dup->x_ul, dup->y_ul, 0.0,
				dup->cellsizex);
		set_mv_float(&(dup->misval));
	}
	if (gl_rowwise)
		RputAllMV(dup->CSF_MAP);
	return dup;
}
#endif

#ifdef HAVE_LIBGIS

static GRIDMAP *read_grass(GRIDMAP * m)
{
/* from: r.out.arc source, by M. Neteler neteler@geog.uni-hannover.de */

	void *raster, *ptr;
	RASTER_MAP_TYPE out_type, map_type;
	char *name;
	char *mapset;
	int fd;
	int row, col;
	struct Cell_head region;

	DUMP(": trying grass map format... ");
	if (! grass()) {
		DUMP("(no grass session recognized) ");
		return NULL;
	}

	name = (char *) m->filename;
	if ((mapset = G_find_cell2(name, "")) == NULL)
		return NULL;
		/*
		G_fatal_error("%s: <%s> cellfile not found\n", G_program_name(), name);
		*/

	map_type = G_raster_map_type(name, mapset);
	out_type = map_type;

	fd = G_open_cell_old(name, mapset);
	if (fd < 0)
		exit(1);

	/* null_row = G_allocate_null_buf(); */
	raster = G_allocate_raster_buf(out_type);

	G_get_window(&region);

	fprintf(stdout, "ncols %d\n", region.cols);
	fprintf(stdout, "nrows %d\n", region.rows);

	m->rows = region.rows;
	m->cols = region.cols;
	m->x_ul = region.west;
	m->y_ul = region.north;
	m->cellsizex = fabs(region.east - region.west) / region.cols;
	m->cellsizey = fabs(region.north - region.south) / region.rows;

	alloc_mv_grid(m);

	if (G_projection() == 3) {	/* Is Projection != LL (3) */
		/* EJP: do something here? */
		/*
		   G_format_easting (region.west, buf, region.proj);
		   fprintf (stdout, "xllcorner %s\n", buf);
		   G_format_northing (region.south, buf, region.proj);
		   fprintf (stdout, "yllcorner %s\n", buf);
		 */
	}

	fd = G_open_cell_old(name, mapset);
	if (fd < 0)
		exit(1);

	for (row = 0; row < m->rows; row++) {
		if (G_get_raster_row(fd, raster, row, out_type) < 0)
			exit(1);
		/*
		   if (G_get_null_value_row(fd, null_row, row) < 0)
		   exit(1);
		 */
		for (col = 0, ptr = raster; col < m->cols; col++,
			 ptr = G_incr_void_ptr(ptr, G_raster_size(out_type))) {
			if (!G_is_null_value(ptr, out_type)) {
				if (out_type == CELL_TYPE)
					map_put_cell(m, row, col, (int) *((CELL *) ptr));
				else if (out_type == FCELL_TYPE)
					map_put_cell(m, row, col, (float) *((FCELL *) ptr));
				else if (out_type == DCELL_TYPE) {
					map_put_cell(m, row, col, (double) *((DCELL *) ptr));
				}
			}
		}						/* for col */
	}							/* for row */
	G_close_cell(fd);
	m->write = write_grass;
	m->type = MT_GRASS;
	m->is_binary = 1;
	m->celltype = CT_IEEE4;
	return m;
}

static GRIDMAP *write_grass(GRIDMAP * m)
{
/* from: r.in.arc source, by M. Neteler neteler@geog.uni-hannover.de */

	char *output;
	char *title;
	int cf;
	struct Cell_head cellhd;
	CELL *cell = NULL;
	FCELL *fcell = NULL;
	DCELL *dcell = NULL;
	int row, col;
	int nrows, ncols;
	int rtype;
	double x;
	char *err;

	output = (char *) m->filename;
	title = m->description;
	if (m->celltype == CT_UINT8)
		rtype = CELL_TYPE;
	else
		rtype = FCELL_TYPE;

	cellhd.zone = G_zone();
	cellhd.proj = G_projection();
	cellhd.cols = m->cols;
	cellhd.rows = m->rows;
	cellhd.ew_res = m->cellsizex;
	cellhd.ns_res = m->cellsizey;
	cellhd.west = m->x_ul;
	cellhd.east = cellhd.west + (cellhd.ew_res * cellhd.cols);
	cellhd.north = m->y_ul;
	cellhd.south = cellhd.north - (cellhd.ns_res * cellhd.rows);

	if ((err = G_adjust_Cell_head(&cellhd, 1, 1)) != 0)
		G_fatal_error("%s: G_adjust_Cell_head failed", G_program_name());

	nrows = m->rows;
	ncols = m->cols;

	if (G_set_window(&cellhd) < 0)
		exit(3);

	if (nrows != G_window_rows())
		G_fatal_error("%s: OOPS -- rows changed from %d to %d\n", nrows, 
				G_program_name(), G_window_rows());

	if (ncols != G_window_cols())
		G_fatal_error("%s: OOPS -- cols changed from %d to %d\n", ncols,
				G_program_name(), G_window_cols());

	switch (rtype) {
	case CELL_TYPE:
		cell = G_allocate_c_raster_buf();
		break;
	case FCELL_TYPE:
		fcell = G_allocate_f_raster_buf();
		break;
	case DCELL_TYPE:
		dcell = G_allocate_d_raster_buf();
		break;
	}
	if ((cf = G_open_raster_new(output, rtype)) < 0)
		G_fatal_error("%s: unable to create raster map %s", 
				G_program_name(), output);

	for (row = 0; row < nrows; row++) {
		G_percent(row, nrows, 5);
		for (col = 0; col < ncols; col++) {
			if (map_cell_is_mv(m, row, col)) {
				switch (rtype) {
				case CELL_TYPE:
					G_set_c_null_value(cell + col, 1);
					break;
				case FCELL_TYPE:
					G_set_f_null_value(fcell + col, 1);
					break;
				case DCELL_TYPE:
					G_set_d_null_value(dcell + col, 1);
					break;
				}
			} else {
				x = map_get_cell(m, row, col);
				switch (rtype) {
				case CELL_TYPE:
					cell[col] = (CELL) x;
					break;
				case FCELL_TYPE:
					fcell[col] = (FCELL) x;
					break;
				case DCELL_TYPE:
					dcell[col] = (DCELL) x;
					break;
				}
			}
		}
		switch (rtype) {
		case CELL_TYPE:
			G_put_c_raster_row(cf, cell);
			break;
		case FCELL_TYPE:
			G_put_f_raster_row(cf, fcell);
			break;
		case DCELL_TYPE:
			G_put_d_raster_row(cf, dcell);
			break;
		}
	}
	fprintf(stderr, "CREATING SUPPORT FILES FOR %s\n", output);
	G_close_cell(cf);
	if (title)
		G_put_cell_title(output, title);
	return m;
}

#endif

#ifdef HAVE_LIBGDAL
static GRIDMAP *read_gdal(GRIDMAP *m) {
	double adfMinMax[2];
	int i, j, bGotMin, bGotMax;
	GDALRasterBandH  hBand;

	/* GDALAllRegister(); -- has been done in main() */
	DUMP(m->filename);
	DUMP(": trying one of GDAL formats... ");

	CPLPushErrorHandler(CPLQuietErrorHandler);
	/* if ((m->hDataset = GDALOpen(m->filename, GA_ReadOnly)) == NULL)
		return NULL; */
	if ((m->hDataset = GDALOpen(m->filename, GA_Update)) == NULL)
		return NULL;
	CPLErrorReset();
	CPLPopErrorHandler();

	DUMP("yes!\n");
	m->hDriver = GDALGetDatasetDriver(m->hDataset);
	if (DEBUG_NORMAL)
		printlog( "GDAL driver: %s [%s]\n", 
			GDALGetDriverShortName(m->hDriver),
			GDALGetDriverLongName(m->hDriver));

	if (GDALGetRasterCount(m->hDataset) > 1)
		pr_warning("only band 1 in %s is used", m->filename);

	if (GDALGetProjectionRef(m->hDataset) != NULL && DEBUG_NORMAL) 
		printlog("Projection is `%s'\n", GDALGetProjectionRef(m->hDataset));

	if (GDALGetGeoTransform(m->hDataset, m->GeoTransform) == CE_None) {
		if (DEBUG_NORMAL) {
			printlog("Origin = (%g,%g), ",
				m->GeoTransform[0], m->GeoTransform[3]);
			printlog("Pixel Size = (%g,%g), ",
				m->GeoTransform[1], m->GeoTransform[5]);
		}
		m->x_ul = m->GeoTransform[0];
		m->y_ul = m->GeoTransform[3];
		m->cellsizex = fabs(m->GeoTransform[1]);
		m->cellsizey = fabs(m->GeoTransform[5]); 
		/* funny, but may be neg.!! */
		/* 
		printf("cellsize: %g x %g\n", m->cellsizex, m->cellsizey);
		*/
	} else
		ErrMsg(ER_IMPOSVAL, "no GDALGetGeoTransform information");

	hBand = GDALGetRasterBand(m->hDataset, 1);
	adfMinMax[0] = GDALGetRasterMinimum(hBand, &bGotMin);
	adfMinMax[1] = GDALGetRasterMaximum(hBand, &bGotMax);
	if(! (bGotMin && bGotMax))
	    GDALComputeRasterMinMax( hBand, TRUE, adfMinMax );

	if (DEBUG_NORMAL)
		printlog( "Min=%g, Max=%g\n", adfMinMax[0], adfMinMax[1] );
	m->cellmin = adfMinMax[0];
	m->cellmax = adfMinMax[1];
	
	if(GDALGetOverviewCount(hBand) > 0 && DEBUG_DUMP)
	    printlog( "Band has %d overviews.\n", GDALGetOverviewCount(hBand));

	if(GDALGetRasterColorTable( hBand ) != NULL && DEBUG_DUMP)
	    printf("Band has a color table with %d entries.\n", 
		     GDALGetColorEntryCount(GDALGetRasterColorTable(hBand)));

	m->rows = GDALGetRasterBandYSize(hBand);
	m->cols = GDALGetRasterBandXSize(hBand);
	alloc_mv_grid(m);
	m->misval = GDALGetRasterNoDataValue(hBand, NULL);
	for (i = 0; i < m->rows; i++) {
		if (GDALRasterIO(hBand, GF_Read, 0, i, m->cols, 1, 
			m->grid[i], m->cols, 1, GDT_Float32, 
			0, 0 ) == CE_Failure)
				pr_warning("error on reading %s\n", m->filename);
		for (j = 0; j < m->cols; j++)
			if (m->grid[i][j] == m->misval)
				set_mv_float(&(m->grid[i][j]));
	}
	m->type = MT_GDAL;
	m->is_binary = 1; /* not necessarily */
	m->celltype = CT_IEEE4;
	/* 
	GDALClose(m->hDataset);
	*/
	m->write = write_gdal;
	return m;
}

static GRIDMAP *write_gdal(GRIDMAP *m) {
	/* NOTE this function still fails on anything */
	GDALRasterBandH  hBand;
	int i, j;
	unsigned char b;

	/*
	if (GDALGetRasterAccess(m->hDataset) != GA_Update)
		printf("NO ACCESS\n");
	*/
	hBand = GDALGetRasterBand(m->hDataset, 1);
	if (m->celltype != CT_UINT8) {
		GDALSetRasterNoDataValue(hBand, m->misval);
		for (i = 0; i < m->rows; i++) { /* one row at a time */
			for (j = 0; j < m->cols; j++) /* { */
				if (is_mv_float(&(m->grid[i][j])))
					m->grid[i][j] = m->misval;
				if (GDALRasterIO(hBand, GF_Write, 0, i, m->cols, 1, 
					m->grid[i], m->cols, 1, GDT_Float32, 0, 0 ) == CE_Failure)
					pr_warning("error on writing %s\n", m->filename);
			/* 
				if (GDALRasterIO(hBand, GF_Write, j, i, 1, 1, 
					&(m->grid[i][j]), 1, 1, GDT_Float32, 
					0, 0 ) == CE_Failure)
						pr_warning("error on writing %s\n", m->filename);
					*/
				/* printf("value written: [%d,%d] %g\n", i, j, m->grid[i][j]); */
			/* } */
		}
	} else {
		GDALSetRasterNoDataValue(hBand, (double) 0xFF);
		for (i = 0; i < m->rows; i++) {
			for (j = 0; j < m->cols; j++) {
				if (is_mv_float(&(m->grid[i][j])))
					b = 0xFF;
				else
					b = (unsigned char) m->grid[i][j];
				if (GDALRasterIO(hBand, GF_Write, j, i, 1, 1, 
					&b, sizeof(unsigned char), 1, GDT_Byte, 
					0, 0) == CE_Failure)
						pr_warning("error on writing %s\n", m->filename);
			}
		}
	}
	/*
	GDALSetRasterStatistics(hBand, m->cellmin, m->cellmax, 0.0, 0.0);
	*/
	GDALClose(m->hDataset);
	return m;
}
#endif

#ifdef HAVE_LIBNETCDF
# include "netcdf.h"
/* Macro to create a NaN */

# ifdef __alpha
# define mknan(x) (((unsigned int *) &x)[0] = 0x7fffffff,\
  ((unsigned int *) &x)[1] = 0xffffffff)
# endif

# ifdef sony
# define mknan(x) (((unsigned long *) &x)[0] = 0x7ff7ffff,\
  ((unsigned long *) &x)[1] = 0x7ff7ffff)
# endif

# if !defined(mknan)
# define mknan(x) (((unsigned long *) &x)[0] = 0x7fffffff,\
  ((unsigned long *) &x)[1] = 0xffffffff)
# endif

/* GMT map support, contributed by Konstantin Malakhanov,
(mailto:kosta@iwwnt.iwwl.rwth-aachen.de) who wrote about this:

And last but not least: I write here all limitations of using GMT 
grids with GSTAT :-

a) GMT grids can be centered at pixels  or at nodes. GSTAT grids are
centered at pixels, so node-wise GSTAT grids will be converted to
pixel-centered.  GMT grids from GSTAT are always pixel-centered.
Convertation could be made in two ways: either decrease the number of
rows and columns by one and set pixel values to mean (or median, or
what you like at most) value of 4 nodes (this changes values, but
preserves boundaries of grid) or extent grid limits to half cell
size to west/east/north/south and use nodes as centers of new pixels
(this preserves values, but slightly changes the limits of grids). I
implemented the second way, so if you have GMT node-centered grid as
mask, then the extensions of result grid from GSTAT will be one cell
size bigger in X- and Y-directions!

b) GMT grids can have multiplication factor and an add offset for
z-values. As GSTAT grid definition does not allow for that, GMT grid
values  with factor different from 1.0 and/or value offset different
from 0.0 will be accordingly transformed during the loading.
(Comment: for reasons I cannot understand GMT grids have sometimes
factor 0.0 which makes no sense. So factor==0.0 will 
be treated as 1.0). GSTAT grids in GMT format always have factor== 
1.0 and value offset==0.0.

c) GMT system allows for rectangular coordinate system or for
geographical projections, but there is no way to detect it from grid
itself (in GMT commands , projection is almost always supplied as
one of arguments by user). So the way GMT grids are treated is
defined by user and not stored in a grid. That means that using
GSTAT for grids which are supposed to have longitude/lattitude
coordinates WILL give results, but these results are useless as 
spheroid of Earth is not taken into account and  it means by no way 
that GSTAT can interpolate/simulate over sphere in geographical 
projections (if one needs such things, take a look at Spherekit at 
http://www.ncgia.ucsb.edu/pubs/spherekit). So I didn't follow this 
branch further.

d) GMT grid definition has fields for names of x-,y-,z-units. These 
fields are ignored at reading and will be set to " " in the result 
grid.

e) GMT grids can have complex z-values. This is neither checked for
nor used!  
*/

/* This are all modified/added things to read/write GMT grids:
   map_read and map_convert are changed just to handle GMT case,
   nothing more else added there. New functions are: read_gmt,
   write_gmt, cdf_read_grd_info.

   You also need header file gmt.h AND netCDF library (by ftp at
   unidata.ucar.edu:/pub/netcdf). You have to compile with additional
   flags -I<path to netcdf #include's> -L<path to netcdf libraries>
   -lnetcdf.

  the copyright part of the README in the GMT 3.0 distribution follows:

  Copyright (c) 1991-1995, P. Wessel & W. H. F. Smith

  Permission to use, copy, modify, and distribute this software and its
  documentation for any purpose without fee is hereby granted, provided
  that the above copyright notice appear in all copies, that both that
  copyright notice and this permission notice appear in supporting
  documentation, and that the name of GMT not be used in
  advertising or publicity pertaining to distribution of the software
  without specific, written prior permission.  The University of Hawaii
  (UH) and the National Oceanic and Atmospheric Administration (NOAA) make no
  representations about the suitability of this software for any purpose.
  It is provided "as is" without expressed or implied warranty.  It is
  provided with no support and without obligation on the part of UH or
  NOAA, to assist in its use, correction, modification, or enhancement.
*/

/*--------------------------------------------------------------------
 *    The GMT-system:   @(#)gmt_cdf.c   2.19  2/23/95
 *
 *    Copyright (c) 1991-1995 by P. Wessel and W. H. F. Smith
 *--------------------------------------------------------------------*/
/*
 *
 *      G M T _ C D F . C   R O U T I N E S
 *
 * Takes care of all grd input/output built on NCAR's netCDF routines (which is
 * an XDR implementation)
 * Most functions return the ncerr error value which will be non-zero if
 * an error occured.
 *
 * Author:      Paul Wessel
 * Date:	9-SEP-1992
 * Modified:    27-JUN-1995
 * Version:     3.0
 *
 * Functions include:
 *
 *      cdf_read_grd_info :     Read header from file
 *      cdf_read_grd :	  Read header and data set from file
 *      cdf_write_grd_info :    Update header in existing file
 *      cdf_write_grd :	 Write header and data set to new file
 *
 * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/* Modified from gmt_cdf.c for GSTAT : */
/* GSTAT grid is stored pixel-wise (in the terms of GMT grids), so */
/* "normal node" GMT grid are converted on the fly. GMT grids written by */
/* GSTAT are always pixel-wise. */

/* GMT grid definition supports also factor and offset for z-values. These */
/* will be use to convert GMT grid values, GMT grids written by GSTAT */
/* have always factor==1.0 and offset==0.0. MISVAL will be set to GMT NaN */

static int cdf_read_grd_info(GRIDMAP * m, double *x_max, double *y_min,
							 double *z_scale_factor, double *z_add_offset,
							 int *node_offset)
{

	long cdfid, nm[2], start[2], edge[2];
	double dummy[2], gmt_misval;
	/* double y_inc; */

	char text[480];
	long x_range_id, y_range_id, z_range_id, inc_id, nm_id, z_id;

	/* Open file and get info */

	if (ncopts)
		ncopts = 0;

	if ((cdfid = ncopen(m->filename, NC_NOWRITE)) == -1)
		return (-1);

	/* Get variable ids */

	x_range_id = ncvarid(cdfid, "x_range");
	y_range_id = ncvarid(cdfid, "y_range");
	z_range_id = ncvarid(cdfid, "z_range");
	inc_id = ncvarid(cdfid, "spacing");
	nm_id = ncvarid(cdfid, "dimension");
	z_id = ncvarid(cdfid, "z");

	/* Get attributes */

/*      ncattget (cdfid, x_range_id, "units", (void *)header->x_units); */
/*      ncattget (cdfid, y_range_id, "units", (void *)header->y_units); */
/*      ncattget (cdfid, z_range_id, "units", (void *)header->z_units); */
/*      ncattget (cdfid, z_id, "long_name", (void *)header->title); */
	ncattget(cdfid, z_id, "scale_factor", (void *) z_scale_factor);
	ncattget(cdfid, z_id, "add_offset", (void *) z_add_offset);
	ncattget(cdfid, z_id, "node_offset", (void *) node_offset);
/*      ncattget (cdfid, NC_GLOBAL, "title", (void *)header->title); */
	ncattget(cdfid, NC_GLOBAL, "source", (void *) text);
	m->description = string_dup(text);
/*      strncpy (m->history, &text[320], 160); */

	/* Get variables */

	start[0] = 0;
	edge[0] = 2;

	ncvarget(cdfid, x_range_id, start, edge, (void *) dummy);
	m->x_ul = dummy[0];
	*x_max = dummy[1];
	ncvarget(cdfid, y_range_id, start, edge, (void *) dummy);
	*y_min = dummy[0];
	m->y_ul = dummy[1];
	ncvarget(cdfid, inc_id, start, edge, (void *) dummy);
	m->cellsizex = dummy[0];
	m->cellsizey = dummy[1];
	ncvarget(cdfid, nm_id, start, edge, (void *) nm);
	m->cols = nm[0];
	m->rows = nm[1];
	ncvarget(cdfid, z_range_id, start, edge, (void *) dummy);
	m->cellmin = dummy[0];
	m->cellmax = dummy[1];

	if (*node_offset == 0) {	/* node grid */
		m->x_ul = m->x_ul - 0.5 * m->cellsizex;
		*x_max = *x_max + 0.5 * m->cellsizex;
		m->y_ul = m->y_ul + 0.5 * m->cellsizey;
		*y_min = *y_min - 0.5 * m->cellsizey;
	}
	m->base_size = m->rows * m->cols;
	mknan(gmt_misval);

	m->misval = gmt_misval;		/* DEFAULT_MISVAL; */
	m->is_binary = 1;
	ncclose(cdfid);

	return (ncerr);
}

static GRIDMAP *read_gmt(GRIDMAP * m)
{
	long cdfid, z_id, /* j2, one_or_zero, *k, */ start[2], edge[2];
	/* long first_col, last_col, first_row, last_row; */
	long i, j, ij /* ,width_in, width_out, height_in, i_0_out, inc = 1 */ ;

/* 	char *file;  */
	float tmp;
	double x_max, y_min, z_scale_factor, z_add_offset;
	int node_offset;
	double /* off, half_or_zero,x,w,e,s,n, */ gmt_misval;

	DUMP(m->filename);
	DUMP(": trying GMT netCDF map format... ");
	if (cdf_read_grd_info
		(m, &x_max, &y_min, &z_scale_factor, &z_add_offset,
		 &node_offset) == -1)
		return NULL;

	m->type = MT_GMT;
	if (sizeof(nclong) != sizeof(long))
		pr_warning("please read comment in source file %s, line %d\n",
				   __FILE__, __LINE__ + 1);
	/* WARNING: the original code used type nclong for all variables
	   in the gmt read/write functions. However, hp-cc seriously complained
	   about this, since they were cast to longs. We changed this to long,
	   but don't know what happens if int's and longs differ in size. So,
	   if you get strange results, please try to compile with all long's
	   changed to nclong.  */

	if (ncopts)
		ncopts = 0;
	if ((cdfid = ncopen(m->filename, NC_NOWRITE)) == -1)
		ErrMsg(ER_IMPOSVAL,
			   "read_gmt: GMT Fatal Error: cannot read file - exiting\n");

	alloc_mv_grid(m);
	/* Get variable id */

	z_id = ncvarid(cdfid, "z");

	m->cellmin = FLT_MAX;
	m->cellmax = -FLT_MAX;

	edge[0] = m->cols;
	mknan(gmt_misval);

	for (i = 0; i < m->rows; i++) {
		ij = start[0] = i * m->cols;
		ncvarget(cdfid, z_id, start, edge, (void *) &m->base[ij]);
		for (j = 0; j < m->cols; j++) {
			tmp = m->grid[i][j];
			if (tmp == gmt_misval) {
				/* set_mv_float(&(m->grid[i][j])); */
			} else {
				if ((z_scale_factor == 1.0 || z_scale_factor == 0.0)
					&& z_add_offset == 0.0) {
				} else
					m->grid[i][j] =
						(float) (tmp * z_scale_factor + z_add_offset);
				m->cellmin = MIN(m->cellmin, m->grid[i][j]);
				m->cellmax = MAX(m->cellmax, m->grid[i][j]);
			}
		}
	}
	ncclose(cdfid);
	m->write = write_gmt;
	return (m);
}

static GRIDMAP *write_gmt(GRIDMAP * m)
{
/*      int pad[4]; */
	int ij, i /* ,j */ ;
	long start[2], edge[2];
	long cdfid, nm[2];
	long /* width_in, */ width_out, height_out /* ,one_or_zero */ ;
	/* long first_col; */
	char *temp = " ";

	double dummy[2];
	double w, e, s, n, tmp;

	char text[480];

	/* dimension ids */
	long side_dim, xysize_dim;

	/* variable ids */
	long x_range_id, y_range_id, z_range_id, inc_id, nm_id, z_id;

	/* variable shapes */
	int dims[1];

	w = m->x_ul;
	e = m->x_ul + m->cols * m->cellsizex;
	s = m->y_ul - m->rows * m->cellsizey;
	n = m->y_ul;

	/* Create file and enter define mode */

	if (ncopts)
		ncopts = 0;
	if ((cdfid = nccreate(m->filename, NC_CLOBBER)) == -1)
		ErrMsg(ER_WRITE, m->filename);

	/* Get dimension of subregion to write */

	width_out = m->cols;		/* rint ((e - w) / m->cellsize) + one_or_zero;  */
	height_out = m->rows;		/* rint ((n - s) / m->cellsize) + one_or_zero;  */

	/* define dimensions */

	side_dim = ncdimdef(cdfid, "side", 2);
	xysize_dim = ncdimdef(cdfid, "xysize", (int) m->base_size);
		/* (width_out * height_out)) */ ;

	/* define variables */

	dims[0] = side_dim;
	x_range_id = ncvardef(cdfid, "x_range", NC_DOUBLE, 1, dims);
	y_range_id = ncvardef(cdfid, "y_range", NC_DOUBLE, 1, dims);
	z_range_id = ncvardef(cdfid, "z_range", NC_DOUBLE, 1, dims);
	inc_id = ncvardef(cdfid, "spacing", NC_DOUBLE, 1, dims);
	nm_id = ncvardef(cdfid, "dimension", NC_LONG, 1, dims);

	dims[0] = xysize_dim;
	z_id = ncvardef(cdfid, "z", NC_FLOAT, 1, dims);

	/* assign attributes */

	ncattput(cdfid, x_range_id, "units", NC_CHAR, 80, (void *) temp);	/* header->x_units); */
	ncattput(cdfid, y_range_id, "units", NC_CHAR, 80, (void *) temp);	/* header->y_units); */
	ncattput(cdfid, z_range_id, "units", NC_CHAR, 80, (void *) temp);	/* header->z_units); */
	ncattput(cdfid, z_id, "long_name", NC_CHAR, 80,
			 (void *) &m->description);
	tmp = 1.0;
	ncattput(cdfid, z_id, "scale_factor", NC_DOUBLE, 1, (void *) &tmp);	/* &header->z_scale_factor); */
	tmp = 0.0;
	i = 1;
	ncattput(cdfid, z_id, "add_offset", NC_DOUBLE, 1, (void *) &tmp);	/* &header->z_add_offset); */
	ncattput(cdfid, z_id, "node_offset", NC_LONG, 1, (void *) &i);	/* &header->node_offset); */
	ncattput(cdfid, NC_GLOBAL, "title", NC_CHAR, 80,
			 (void *) &m->description);

	strncpy(text, m->history ? m->history : "(empty)", 480);
	ncattput(cdfid, NC_GLOBAL, "source", NC_CHAR, 480, (void *) text);

	/* leave define mode */

	ncendef(cdfid);

	/* store header information and array */

	start[0] = 0;
	edge[0] = 2;
	dummy[0] = w;				/* header->x_min; */
	dummy[1] = e;				/* header->x_max; */
	ncvarput(cdfid, x_range_id, start, edge, (void *) dummy);
	dummy[0] = s;				/* header->y_min; */
	dummy[1] = n;				/* header->y_max; */
	ncvarput(cdfid, y_range_id, start, edge, (void *) dummy);
	dummy[0] = m->cellsizex;
	dummy[1] = m->cellsizey;
	ncvarput(cdfid, inc_id, start, edge, (void *) dummy);
	nm[0] = width_out;
	nm[1] = height_out;
	ncvarput(cdfid, nm_id, start, edge, (void *) nm);
	dummy[0] = m->cellmin;
	dummy[1] = m->cellmax;
	ncvarput(cdfid, z_range_id, start, edge, (void *) dummy);
	edge[0] = width_out;

	for (i = 0; i < m->rows; i++) {
		ij = start[0] = i * m->cols;
		ncvarput(cdfid, z_id, start, edge, (void *) &m->base[ij]);
		/* convert to gmt_misval? seems not necessary */
	}
	ncclose(cdfid);

	return (m);
}
#endif

#ifdef HAVE_T2_GRIDFORMAT

GRIDMAP *read_T2(GRIDMAP * m)
{
	FILE *f;
	int line = 1, i, j;
	char *cp;
	float fl;

	DUMP(m->filename);
	DUMP(": trying T2 format... ");
	if ((f = fopen(m->filename, "r")) == NULL)
		return NULL;
	while (line <= 5 && get_line(&line_buf, &line_size, f) != NULL) {
		cp = line_buf + 25;		/* ascii text */
		switch (line) {
		case 1:
			if (strstr(line_buf, "FILETYPE DATATYPE VERNO:") != line_buf) {
				efclose(f);
				return NULL;
			}
			break;
		case 2:
			if (sscanf(cp, "%d %d %lg %lg %lg", &(m->cols), &(m->rows),
					   &(m->cellsizex), &(m->x_ul), &(m->y_ul)) != 5)
				ErrMsg(ER_READ, "T2 header line 2");
			m->cellsizey = m->cellsizex;
			m->y_ul += m->cellsizey * m->rows;
			break;
		case 3:
			if (sscanf(cp, "%g %g %g", &(m->misval), &fl, &fl) != 3)
				ErrMsg(ER_READ, "T2 header line 3");
			break;
		case 4:
			if (sscanf(cp, "%g %g %g %g",
					   &(m->cellmin), &(m->cellmax), &fl, &fl) != 4)
				ErrMsg(ER_READ, "T2 header line 4");
			break;
		case 5:
			m->description = string_dup(line_buf);
			break;
		}
		line++;
	}
	alloc_mv_grid(m);
	for (i = 0; i < m->rows; i++) {
		if (fscanf(f, "%d", &line) != 1 || line != m->rows - i)
			ErrMsg(ER_READ, "T2 file, reading line number");
		for (j = 0; j < m->cols; j++) {
			if (fscanf(f, "%g", &fl) != 1)
				ErrMsg(ER_READ, "T2 file, reading float");
			if (fl != m->misval)
				m->grid[i][j] = fl;
		}
	}
	efclose(f);
	m->type = MT_T2;
	m->write = write_T2;
	DUMP("yes\n");
	return m;
}

static GRIDMAP *write_T2(GRIDMAP * m)
{
	FILE *f;
	int i, j;

	f = efopen(m->filename, "w");
	fprintf(f, "FILETYPE DATATYPE VERNO:  %d %d %d\n", 22, 67, 523);
	fprintf(f, "NX NY DIM XORIG YORIG  :  %d %d  %13E  %13E  %13E\n",
			m->cols, m->rows, m->cellsizex, m->x_ul,
			m->y_ul - m->cellsizey * m->rows);
	fprintf(f, "DELETE UTMZONE ORIENT  :  %13E 0 0\n", m->misval);
	fprintf(f, "MIN MAX MEAN ST.DEV    :  %13E %13E  0  0\n",
			m->cellmin, m->cellmax);
	if (m->description)
		fprintf(f, "%s", m->description);
	else
		fprintf(f, "no description\n");
	for (i = 0; i < m->rows; i++) {
		fprintf(f, "%7d", m->rows - i);
		for (j = 0; j < m->cols; j++)
			fprintf(f, "%s  %13E", j % 5 == 0 ? "\n" : "",
					is_mv_float(&(m->grid[i][j])) ? m->misval : m->
					grid[i][j]);
		fprintf(f, "\n");
	}
	return m;
}
#endif							/* HAVE_T2_GRIDFORMAT */
/************************ ER-Mapper I/O support routines ******************
Contributed by Steve Joyce, who wrote on Fri, 2 Jan 1998:

ER-Mapper support:  I've got it, it works for gstat-2.0, and you are
welcome to it.  The changes are all to mapio.c/h, with the exception of
one line in data.c/h.

Some notes:
ER-Mapper provides a c-function library for programmers to read and
write datafiles in a standard way.  Maybe you remember my first version
used this library, but linking gstat together with ermapper.lib was just
an enormous pain in the ass.  It turns out, ermapper raster files are
fairly straightforward anyway, with a separate ASCII header and binary
data file. So the current version reads and writes the ermapper files
 directly without using ermapper.lib.  This may cause it to fail for
future versions of ER-Mapper, but I can live with that.

ER-Mapper raster files can be multi-channel-  I check the number of
channels on input files and bail out if there is more than one.  Maybe
sometime we can set up the data specification syntax to incluyde a
channel as you talked about before.

ER-Mapper files can specify to skip a number of bytes in the binary
file- I don't take care of this,  but check the file size consistencey
of the binary file.

Ermapper header binary files can have a different root name than the
header- I don't take care of this, and force them to be the same.

ER-Mapper data types can be signed or unsigned chars, ints (16 or 32
bits), 4-byte real or 8-byte double.  I read all formats and cast into
4-byte real for sorage in the gstat gridmap.  I always write 4-byte
real.

Byte order can be specified in ermapper raster files-  I check for it
and reorder as necessary.

The ER-Mapper coordinate origin can be at an arbitrary fractional pixel
position- I correct it to be the upper left corner of the upper left
pixel in the grid.

ER-Mapper coordinates can be either Easting/Northing, Long/Lat, or
Raw(X-Y).  I do a strict check for 'EN' coordinates, because they match
the definition used in gridmaps.  RAW coordinates have positive 'y'
going down (same as cell coordinates) whereas EN coordinates have
positive 'y' going up.
***********************************************************************/

static void swap_multiformat(unsigned char *b, unsigned int m,
							 unsigned int n)
{
/* little <-> big endian swapping of n elements of size m in b */
	unsigned char tmp;
	int i, j;

	/* Order: 0123 => 3210 */
	if (m > 1) {
		for (i = 0; i < n; i++) {
			for (j = 0; j < m / 2; j++) {
				tmp = b[j];
				b[j] = b[m - j - 1];
				b[m - j - 1] = tmp;
			}
			b += m;
		}
	}
}

static unsigned int sizeof_ct(CellType ct)
{
	switch (ct) {
	case CT_UINT8:
		return sizeof(UINT8);
	case CT_UINT16:
		return sizeof(UINT16);
	case CT_UINT32:
		return sizeof(UINT32);
	case CT_INT8:
		return sizeof(INT8);
	case CT_INT16:
		return sizeof(INT16);
	case CT_INT32:
		return sizeof(INT32);
	case CT_IEEE4:
		return sizeof(IEEE4);
	case CT_IEEE8:
		return sizeof(IEEE8);
	case CT_UNKNOWN:
	default:{
			ErrMsg(ER_IMPOSVAL, "unknown celltype");
			return 0;
		}
	}
}

static int read_multiformat_bgrid(GRIDMAP * m, const char *fname, int swap)
{
	unsigned int i, j, size, filesize, bytespercell;
	FILE *f;
	char *buf;
	CellType celltype = CT_UNKNOWN;

	celltype = m->celltype;
	bytespercell = sizeof_ct(celltype);
	size = bytespercell * m->cols * m->rows;
	filesize = file_size(fname);
	if (filesize != size) {
		pr_warning("size of binary file `%s' is %d (should be %d)",
				   fname, filesize, size);
		ErrMsg(ER_READ, fname);
	}
	if ((f = fopen(fname, "rb")) == NULL)
		return 1;
	alloc_mv_grid(m);

	buf = (char *) emalloc(m->cols * bytespercell);

	for (i = 0; i < m->rows; i++) {
		if (fread(buf, bytespercell, m->cols, f) != m->cols) {
			efclose(f);
			return 1;
		}
		if (swap)
			SWAP_M_N(buf, bytespercell, m->cols);

		for (j = 0; j < m->cols; j++) {
			switch (celltype) {
			case CT_UINT8:
			case CT_NONE:
			case CT_UNKNOWN:
				m->grid[i][j] = (float) ((UINT8 *) buf)[j];
				break;
			case CT_UINT16:
				m->grid[i][j] = (float) ((UINT16 *) buf)[j];
				break;
			case CT_UINT32:
				m->grid[i][j] = (float) ((UINT32 *) buf)[j];
				break;
			case CT_INT8:
				m->grid[i][j] = (float) ((INT8 *) buf)[j];
				break;
			case CT_INT16:
				m->grid[i][j] = (float) ((INT16 *) buf)[j];
				break;
			case CT_INT32:
				m->grid[i][j] = (float) ((INT32 *) buf)[j];
				break;
			case CT_IEEE4:
				m->grid[i][j] = (float) ((IEEE4 *) buf)[j];
				break;
			case CT_IEEE8:		/* possible loss of data! */
				m->grid[i][j] = (float) ((IEEE8 *) buf)[j];
				break;
			}
			if (m->grid[i][j] == m->misval)
				set_mv_float(&(m->grid[i][j]));
			else {
				if (m->cellmin == FLT_MAX)
					m->cellmin = m->cellmax = m->grid[i][j];
				else {
					m->cellmin = MIN(m->cellmin, m->grid[i][j]);
					m->cellmax = MAX(m->cellmax, m->grid[i][j]);
				}
			}
		}
	}							/* for i */
	fclose(f);					/* DON'T use efclose on a fopen'ed file */
	efree(buf); /* KS */
	return 0;
}

static GRIDMAP *read_ermapper(GRIDMAP * m)
{
	FILE *f;

	DUMP(m->filename);
	DUMP(": trying ER-Mapper format... ");
	if (file_exists(string_cat(m->filename, ".ers")) && file_exists(m->filename)) {	/* file and file.ers */
		m->is_binary = 1;
		f = efopen(string_cat(m->filename, ".ers"), "r");
		if (read_ermapper_header(m, f)) {
			efclose(f);
			return NULL;
		}
		if (read_multiformat_bgrid(m, m->filename,
								   m->is_binary == BINARY_NON_NATIVE))
				return NULL;

	} else						/* m->filename and m->filename.hdr don't exist both: */
		return NULL;
	m->type = MT_ERMAPPER;
	m->write = write_ermapper;
	efclose(f);
	return m;
}

static int read_ermapper_header(GRIDMAP * m, FILE * f)
{

	char *tok1 = NULL, *tok2 = NULL;
	unsigned int check = 0, ok = 1;
	int i;
	double xdim, ydim, cellx = 0.0, celly = 0.0;
	int bands;

	/* There is potentially a lot more information in ER-Mapper
	   header files but we only extract the keywords that are relevant
	   to creating the GRIDMAP structure. */

	while (ok && get_line(&line_buf, &line_size, f) != NULL) {
		tok1 = strtok(line_buf, "=");	/* first word */
		tok2 = strtok(NULL, "\n\r");	/* second word */
		if (strstr(tok1, "NrOfCellsPerLine") != NULL) {
			if (read_int(tok2, &i))
				ok = 0;
			m->cols = i;
			check = check | CHECK_COLS;
		} else if (strstr(tok1, "NrOfLines") != NULL) {
			if (read_int(tok2, &i))
				ok = 0;
			m->rows = i;
			check = check | CHECK_ROWS;
		} else if (strstr(tok1, "Eastings") != NULL) {
			if (read_double(tok2, &(m->x_ul)))
				ok = 0;
			check = check | CHECK_X_UL;
		} else if (strstr(tok1, "Northings") != NULL) {
			if (read_double(tok2, &(m->y_ul)))
				ok = 0;
			check = check | CHECK_Y_UL;
		} else if (strstr(tok1, "Xdimension") != NULL) {
			if (read_double(tok2, &(xdim)))
				ok = 0;
		} else if (strstr(tok1, "Ydimension") != NULL) {
			if (read_double(tok2, &(ydim)))
				ok = 0;
		} else if (strstr(tok1, "NullCellValue") != NULL) {
			if (read_float(tok2, &(m->misval)))
				ok = 0;
		} else if (strstr(tok1, "NrOfBands") != NULL) {
			if (read_int(tok2, &(bands)))
				ok = 0;
		} else if (strstr(tok1, "RegistrationCellX") != NULL) {
			if (read_double(tok2, &(cellx)))
				ok = 0;
		} else if (strstr(tok1, "RegistrationCellY") != NULL) {
			if (read_double(tok2, &(celly)))
				ok = 0;
		} else if (strstr(tok1, "CellType") != NULL) {
			if (strstr(tok2, "Unsigned8BitInteger") != NULL)
				m->celltype = CT_UINT8;
			else if (strstr(tok2, "Unsigned16BitInteger") != NULL)
				m->celltype = CT_UINT16;
			else if (strstr(tok2, "Unsigned32BitInteger") != NULL)
				m->celltype = CT_UINT32;
			else if (strstr(tok2, "Signed8BitInteger") != NULL)
				m->celltype = CT_INT8;
			else if (strstr(tok2, "Signed16BitInteger") != NULL)
				m->celltype = CT_INT16;
			else if (strstr(tok2, "Signed32BitInteger") != NULL)
				m->celltype = CT_INT32;
			else if (strstr(tok2, "IEEE4ByteReal") != NULL)
				m->celltype = CT_IEEE4;
			else if (strstr(tok2, "IEEE8ByteReal") != NULL)
				m->celltype = CT_IEEE8;
		} else if (strstr(tok1, "Datum") != NULL) {
			er_datum = string_dup(tok2);
		} else if (strstr(tok1, "Projection") != NULL) {
			er_projection = string_dup(tok2);
		} else if (strstr(tok1, "CoordinateType") != NULL) {
			if (strstr(tok2, "EN") == NULL) {
				pr_warning("CoordinateType `%s' is invalid", tok2);
				ErrMsg(ER_IMPOSVAL, "only EN is supported");
			}
		} else if (strstr(tok1, "ByteOrder") != NULL) {
			if (strstr(tok2, "MSBFirst") != NULL) {
				if (cpu_is_little_endian())
					m->is_binary = BINARY_NON_NATIVE;
				else
					m->is_binary = BINARY_NATIVE;
			} else if (strstr(tok2, "LSBFirst") != NULL) {
				if (cpu_is_little_endian())
					m->is_binary = BINARY_NATIVE;
				else
					m->is_binary = BINARY_NON_NATIVE;
			} else {
				pr_warning("byteorder `%s' not recognized", tok2);
				ErrMsg(ER_IMPOSVAL,
					   "only LSBfirst and MSBfirst are supported");
			}
		}
	}

	m->cellsizex = xdim;
	m->cellsizey = ydim;
	/* don't bail out if we don't have square cells */
	check = check | CHECK_CELLSIZE;

	/*KS added -- EJP: is this old or better? */
	/*
	if (xdim == ydim) {
		m->cellsizex = xdim;
        m->cellsizey = xdim; 
		check = check | CHECK_CELLSIZE;
	}
	*/

	/* no support for multi-channel images yet */
	if (bands != 1) {
		pr_warning("Input file has %d channels", bands);
		ErrMsg(ER_IMPOSVAL, "only single channel files are supported");
	}

	/* ER-Mapper registration cell can be arbitrary location */
	m->x_ul -= cellx * m->cellsizex;
	m->y_ul += celly * m->cellsizey;

	if (check != CHECK_SUM) {
		if (check)
			pr_warning("uncomplete ER-Mapper header in %s", m->filename);
		return 1;
	}
	return 0;
}

static GRIDMAP *write_ermapper(GRIDMAP * m)
{
	FILE *f;

	f = efopen(string_cat(m->filename, ".ers"), "w");
	write_ermapper_header(m, f);
	efclose(f);
	write_binary_grid(m, m->filename, 0);
	return m;
}

static void write_ermapper_header(GRIDMAP * m, FILE * fp)
{
/* This writes a functional ER-Mapper header file for v5.0+             */
/* without using functions from the ermaper.lib programming library.	*/
/* If no projection and Datum are cached from reading an input file, 	*/
/* they are set to "RAW". The proper projection info can be added later.*/
	char *d;

	fprintf(fp, "DatasetHeader Begin\n");
	fprintf(fp, "\tDataSetType\t= ERStorage\n");
	fprintf(fp, "\tDataType\t= Raster\n");
	if (cpu_is_little_endian())
		fprintf(fp, "\tByteOrder\t= LSBFirst\n");
	else
		fprintf(fp, "\tByteOrder\t= MSBFirst\n");
	if (m->history != NULL)
		fprintf(fp, "\tComments\t= \"%s\"\n", m->history);
	fprintf(fp, "\tCoordinateSpace Begin\n");
	if (er_datum != NULL)
		fprintf(fp, "\t\tDatum\t=%s\n", er_datum);
	else
		fprintf(fp, "\t\tDatum\t= \"RAW\"\n");
	if (er_projection != NULL)
		fprintf(fp, "\t\tProjection\t=%s\n", er_projection);
	else
		fprintf(fp, "\t\tProjection\t= \"RAW\"\n");
	fprintf(fp, "\t\tCoordinateType\t= EN\n");
	fprintf(fp, "\t\tRotation\t= 0:0:0.0\n");
	fprintf(fp, "\tCoordinateSpace End\n");
	fprintf(fp, "\tRasterInfo Begin\n");
	fprintf(fp, "\t\tCellType\t= IEEE4ByteReal\n");
	fprintf(fp, "\t\tNullCellValue\t= %f\n", m->misval);
	fprintf(fp, "\t\tCellInfo Begin\n");
	fprintf(fp, "\t\t\tXdimension\t= %f\n", m->cellsizex);
	fprintf(fp, "\t\t\tYdimension\t= %f\n", m->cellsizey);
	fprintf(fp, "\t\tCellInfo End\n");
	fprintf(fp, "\t\tNrOfLines\t= %d\n", m->rows);
	fprintf(fp, "\t\tNrOfCellsPerLine\t= %d\n", m->cols);
	fprintf(fp, "\t\tRegistrationCoord Begin\n");
	fprintf(fp, "\t\t\tEastings\t= %f\n", m->x_ul);
	fprintf(fp, "\t\t\tNorthings\t= %f\n", m->y_ul);
	fprintf(fp, "\t\tRegistrationCoord End\n");
	fprintf(fp, "\t\tNrOfBands\t= 1\n");
	fprintf(fp, "\t\tBandId Begin\n");
	/* strip <CR> from m->description */
	d = m->description;
	while (d && *d != '\0') {
		if (*d == '\n')
			*d = ' ';
		d++;
	}
	fprintf(fp, "\t\t\tValue\t= \"%s\"\n",
			m->description ? m->description : "");
	fprintf(fp, "\t\tBandId End\n");
	fprintf(fp, "\tRasterInfo End\n");
	fprintf(fp, "DatasetHeader End\n");
}

static GRIDMAP *read_surfer(GRIDMAP * m) {
	FILE *f = NULL;
	int i, j, size;
	float *buf /* , zlo, zhi */ ;

	DUMP(m->filename);
	DUMP(".grd: trying Surfer DSAA format... ");
	if (file_exists(string_cat(m->filename, ".grd"))) {
		f = efopen(string_cat(m->filename, ".grd"), "r");
		if (read_surfer_header(m, f)) {
			efclose(f);
			return NULL;
		} else { /* use min/max in header to set missing values later on */
			/*
			zlo = m->cellmin;
			zhi = m->cellmax;
			*/
		}
		if (read_ascii_grid(m, f, 0)) {
			efclose(f);
			return NULL;
		}
	} else	/* filename.grd does not exist */
		return NULL;
	efclose(f);
	/* swap row order: */
	size = m->cols * sizeof(float);
	buf = (float *) emalloc(size);
	for (i = 0, j = m->rows - 1; i < m->rows / 2; i++, j--) {
		/* swap rows i and j: */
		memcpy(buf, m->grid[i], size);
		memcpy(m->grid[i], m->grid[j], size);
		memcpy(m->grid[j], buf, size);
	}
	m->is_binary = 0;
	m->misval = SURFER_MISVAL;
	for (i = 0; i < m->rows; i++)
		for (j = 0; j < m->cols; j++)
			if (m->grid[i][j] >= m->misval)
				set_mv_float(&(m->grid[i][j]));
	m->type = MT_SURFER;
	m->write = write_surfer;
	DUMP("yes\n");
	efree(buf);
	return m;
}

static int read_surfer_header(GRIDMAP * m, FILE * f)
{
	char *buf = 0, *tok1, *tok2;
	int i, buf_size = 0, nx, ny;
	float xlo, xhi, ylo, yhi, zlo, zhi;

	/*
	   DSAA
	   nx ny
	   xlo xhi
	   ylo yhi
	   zlo zhi
	   values ... (row order; starting at ylo ending at yhi)
	 */

	DUMP("well...");
	for (i = 0; i < 5; i++) {
		if (get_line(&buf, &buf_size, f) == NULL)	/* i.e., at EOF */
			return 1;
		tok1 = strtok(buf, " \t\n\r");
		tok2 = strtok(NULL, " \t\n\r");
		if (i == 0) {
			/* DSAA */
			if (tok1 == NULL || tok2 != NULL || strcmp(buf, "DSAA"))
				return 1;
			DUMP("line 1 OK\n");
		} else if (i == 1) {	/* nx ny */
			if (tok1 == NULL || tok2 == NULL || read_int(tok1, &nx) ||
				read_int(tok2, &ny))
				return 1;
			DUMP("line 2 OK\n");
		} else if (i == 2) {	/* xlo xhi */
			if (tok1 == NULL || tok2 == NULL || read_float(tok1, &xlo) ||
				read_float(tok2, &xhi))
				return 1;
			DUMP("line 3 OK\n");
		} else if (i == 3) {	/* ylo yhi */
			if (tok1 == NULL || tok2 == NULL || read_float(tok1, &ylo) ||
				read_float(tok2, &yhi))
				return 1;
			DUMP("line 4 OK\n");
		} else if (i == 4) {	/* zlo zhi */
			if (tok1 == NULL || tok2 == NULL || read_float(tok1, &zlo) ||
				read_float(tok2, &zhi))
				return 1;
			DUMP("line 5 OK\n");
		}
	}
	m->cellsizex = (xhi - xlo) / (nx - 1);	/* [xy]hi/lo point to grid centres */
	m->cellsizey = (yhi - ylo) / (ny - 1);
	m->rows = ny;
	m->cols = nx;
	m->x_ul = xlo - 0.5 * m->cellsizex;
	m->y_ul = yhi + 0.5 * m->cellsizey;
	m->cellmin = zlo;
	m->cellmax = zhi;
	if (buf_size)
		efree(buf);
	DUMP("surfer DSAA header ok...");
	return 0;
}

static GRIDMAP *write_surfer(GRIDMAP * m)
{
	FILE *f;
	int i, n, row, col;
	f = efopen(string_cat(m->filename, ".grd"), "w");
	fprintf(f, "DSAA\n");		/* DSAA */
	fprintf(f, "%d %d\n", m->cols, m->rows);	/* nx ny */
	/* xlo xhi: */
	fprintf(f, gl_format, m->x_ul + 0.5 * m->cellsizex);
	fprintf(f, " ");
	fprintf(f, gl_format, m->x_ul + (m->cols - 0.5) * m->cellsizex);	
	fprintf(f, "\n");
	/* ylo yhi */
	fprintf(f, gl_format, m->y_ul - (m->rows - 0.5) * m->cellsizey);
	fprintf(f, " ");
	fprintf(f, gl_format, m->y_ul - 0.5 * m->cellsizey);	
	fprintf(f, "\n");
	/* zlo zhi */
	fprintf(f, gl_format, m->cellmin);
	fprintf(f, " ");
	fprintf(f, gl_format, m->cellmax);	
	fprintf(f, "\n");
	n = m->rows * m->cols;
	for (i = 0; i < n; i++) {
		row = m->rows - (i / m->cols) - 1;	/* from bottom to top */
		col = i % m->cols;
		if (map_cell_is_mv(m, row, col))
			fprintf(f, "%g", m->misval);
		else
			fprintf(f, "%g", m->grid[row][col]);
		fprintf(f, ((i + 1) % 5 == 0) ? "\n" : " ");
	}
	efclose(f);
	return m;
}

#ifdef WITH_GSLIB
static GRIDMAP *read_gslib(GRIDMAP * m)
{
	FILE *f;
	char *buf = NULL, *tok1, *tok2, *tok3;
	double min[3], siz[3];
	float value;
	int n[3], row, col;
	int buf_size = 0, i;


	DUMP(m->filename);
	DUMP(": trying GSLIB grid format... ");

	if (file_exists(m->filename))
		f = efopen(m->filename, "r");
	else
		return NULL;
	for (i = 0; i < 3; i++) {
		if (get_line(&buf, &buf_size, f) == NULL)	/* i.e., at EOF */
			return NULL;
		tok1 = strtok(buf, " \t\n");
		tok2 = strtok(NULL, " \t\n");
		tok3 = strtok(NULL, " \t\n");
		if (tok1 == NULL || tok2 == NULL || tok3 == NULL)
			return NULL;

		switch (i) {
		case 0:
			if (read_double(tok1, &(min[0])) ||
				read_double(tok2, &(min[1])) ||
				read_double(tok3, &(min[2]))) return NULL;
			break;
		case 1:
			if (read_int(tok1, &(n[0])) ||
				read_int(tok2, &(n[1])) || read_int(tok3, &(n[2])))
				return NULL;
			break;
		case 2:
			if (read_double(tok1, &(siz[0])) ||
				read_double(tok2, &(siz[1])) ||
				read_double(tok3, &(siz[2]))) return NULL;
			break;
		}
	}
	if (n[2] > 1)
		ErrMsg(ER_IMPOSVAL, "cannot read grids with  multiple z-layer");
	if (siz[2] > 0.0)
		pr_warning("z-dimension of GSLIB grid (%g) will be lost", siz[2]);

	m = new_map(READ_ONLY);
	m->cols = n[0];
	m->rows = n[1];
	m->x_ul = min[0] - 0.5 * siz[0];
	m->y_ul = min[1] + (n[1] - 0.5) * siz[1];
	m->cellsizex = siz[0];
	m->cellsizey = siz[1];
	alloc_mv_grid(m);

	/* for (z = 0; z < n[2]; z++) */
	for (row = m->rows - 1; row >= 0; row--) {
		for (col = 0; col < m->cols; col++) {
			if (get_line(&buf, &buf_size, f) == NULL) {	/* i.e., at EOF */
				pr_warning("error in gslib grid %s, row %d, col %d",
						   m->filename, row, col);
				map_free(m);
				return NULL;
			}
			read_float(buf, &value);
			map_put_cell(m, row, col, value);
		}
	}
	if (buf_size)
		efree(buf);
	m->type = MT_GSLIB;
	m->write = write_gslib;
	return m;
}
#endif

static GRIDMAP *write_gslib(GRIDMAP * m)
{
	FILE *f;
	int row, col;

	pr_warning("map %s: setting zsize equal to x-size (%g)", m->filename,
			   m->cellsizex);
	f = efopen(m->filename, "w");
	fprintf(f, "%s\n", m->description ? m->description : m->filename);
	fprintf(f, "%g %g %g\n", m->x_ul + 0.5 * m->cellsizex,
			m->y_ul - (m->rows - 0.5) * m->cellsizey, 0.0);
	fprintf(f, "%d %d %d\n", m->cols, m->rows, 1);
	fprintf(f, "%g %g %g\n", m->cellsizex, m->cellsizey, m->cellsizex);
	for (row = m->rows - 1; row >= 0; row--) {
		for (col = 0; col < m->cols; col++) {
			if (map_cell_is_mv(m, row, col))
				fprintf(f, "%g\n", -1.0e30);
			else
				fprintf(f, "%g\n", map_get_cell(m, row, col));
		}
	}
	efclose(f);
	return m;
}

GRIDMAP *map_switch_type(GRIDMAP *in, MAPTYPE type) {

	switch (type) {
	case MT_ARCGRID:
		in->write = write_arcgrid;
		break;
	case MT_IDRISI:
		in->write = write_idrisi;
		break;
	case MT_IDRISI32:
		in->write = write_idrisi32;
		break;
	case MT_CSF:
#ifdef HAVE_LIBCSF
		in->is_binary = 1;
		in->write = write_csf;
#else
		ErrMsg(ER_IMPOSVAL, "csf maps not avaiable");
#endif
		break;
	case MT_ERMAPPER:
		in->is_binary = 1;
		in->write = write_ermapper;
		break;
	case MT_SURFER:
		in->is_binary = 0;
		in->write = write_surfer;
		break;
	case MT_GNUPLOT:
		in->is_binary = 1;
		in->write = write_gnuplot_binary;
		break;
	case MT_GSLIB:
		in->is_binary = 0;
		in->write = write_gslib;
		break;
	case MT_GRASS:
#ifdef HAVE_LIBGRASS
		in->is_binary = 1;
		in->write = write_grass;
#else
		ErrMsg(ER_IMPOSVAL,
			   "Grass maps not avaiable (compile gstat with grass' libgis)");
#endif
	case MT_GMT:
#ifdef HAVE_LIBNETCDF
		in->is_binary = 1;
		in->write = write_gmt;
#else
		ErrMsg(ER_IMPOSVAL,
			   "GMT maps not avaiable (compile gstat with netcdf library)");
#endif
		break;
	case MT_T2:
#ifdef HAVE_T2_GRIDFORMAT
		in->is_binary = 0;
		in->write = write_T2;
#endif
		break;
	case MT_GDAL:
#ifdef HAVE_LIBGDAL
		in->write = write_gdal;
#endif
		break;
	case MT_UNKNOWN:
		ErrMsg(ER_IMPOSVAL, "map type unknown");
	}
	in->type = type;
	return in;
}

void map_name_nr(GRIDMAP *mask, const char *base, char *name, int nr, int max) {
	char *cp;

	if (mask->type != MT_CSF) {
		if (max > 99999)
			sprintf(name, "%s%06d", base, nr);
		else if (max > 9999)
			sprintf(name, "%s%05d", base, nr);
		else if (max > 999)
			sprintf(name, "%s%04d", base, nr);
		else if (max > 99)
			sprintf(name, "%s%03d", base, nr);
		else
			sprintf(name, "%s%02d", base, nr);
		return;
	}

/*
 * CSF/PCRaster convention:
 * base is the name set in the command file with predictions(id): 'base';
 * name is the output, which will look like: 'base0000.011', if nr is 11
 * Big trick: base0001.000 should follow base0000.999, if nr is 1000.
 * reason for this: csf/PCRaster convention (not my idea).
 */
	nr++; /* start at 1 */
	sprintf(name, "%s", base);
	/* get basename: don't know if this is completely portable */
	if ((cp = strrchr(name, '/')) != NULL) /* unix directory delimiter */
		cp++;
	else if ((cp = strrchr(name, '\\')) != NULL) /* DOS directory delimiter */
		cp++;
	else
		cp = name; /* no slashes: this is the base name */
	sprintf(cp, "%011d", nr); /* print the 11-digit number, like 00000000011 */
	memmove(cp + 9, cp + 8, 4); /* move overlapping last 3 digits + '\0' */
	cp[8] = '.'; /* print 8.3 point */
	memcpy(name, base, strlen(base)); /* overprint name excluding '\0' */
}
#endif

#ifdef MAPIO_LIB /* dummy, test main() */
int main(int argc, char *argv[])
{
	printf("all done\n");
}
#endif
