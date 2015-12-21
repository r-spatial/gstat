/*
 * former mapio io functions, now only skeleton things doing row/col <--> x/y
 */

#include <math.h>				/* floor() */
#include <float.h>				/* FLT_MAX */

#include "defs.h"
#include "glvars.h"
#include "utils.h"
#include "debug.h"
#include "userio.h"
#include "mapio.h"

static GRIDMAP *write_error(GRIDMAP * m);

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
	map->write = write_error;
	map->read_row = map->write_row = NULL;
	map->current_row = 0;
	return map;
}

static GRIDMAP *write_error(GRIDMAP * m)
{
	pr_warning("%s: writing this map format is not supported", m->filename);
	assert(0);
	return NULL;
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
	*x = m->x_ul + (col + 0.5) * m->cellsizex;
	*y = m->y_ul - (row + 0.5) * m->cellsizey;
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
				  unsigned int *col /* output value: pointer to column number */ )
{
	assert(m);
	assert(row);
	assert(col);

	if (x < m->x_ul || x > m->x_ul + m->cols * m->cellsizex ||
		y > m->y_ul || y < m->y_ul - m->rows * m->cellsizey)
		return 1;
	*row = (unsigned int) floor((m->y_ul - y) / m->cellsizey);
	*col = (unsigned int) floor((x - m->x_ul) / m->cellsizex);
	if (*row == m->rows) /* on the bottom edge */
		*row = *row - 1;
	if (*col == m->cols) /* on the right edge */
		*col = *col - 1;
	return 0;
}
