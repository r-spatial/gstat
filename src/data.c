/*
 *  data.c: basic i/o routines on DATA structure 
 */
#include <math.h> /* sqrt */
#include <string.h> /* memcpy */

#include "defs.h"
#include "data.h"
#include "mapio.h"
#include "userio.h"
#include "utils.h"
#include "block.h"
#include "debug.h"
#include "glvars.h"
#include "defaults.h"
#include "mtrx.h"
#include "lm.h" /* free_lm() */
#include "gls.h" /* free_glm() */
#include "nsearch.h"
#include "gcdist.h"

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

static void calc_data_mean_std(DATA *d);
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

static void logprint_data_header(const DATA *d);

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

    if (d->point_ids != NULL) {
        for (i = d->n_list - 1; i >= 0; i--)
            efree(d->point_ids[i]);
	}

	if (d->beta != NULL)
		efree(d->beta);

	efree(d);
	return;
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
	data->prob = 1.0;
	data->skip = 0;
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
	printlog("\nidx x:%s;", d->x_coord);
	printlog("y:%s;", d->y_coord);
	printlog("z:%s;", d->z_coord);
	printlog("v:%s;\n", d->variable);
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
    if (d->point_ids) {
        printlog("ID: %s ", d->point_ids[GET_INDEX(p)]);
	}
	printlog("\n");
}

void push_point(DATA *d, const DPOINT *p) {
 	int i;
/*
 * add one point p to the data structure d
 * [counts on the fact that erealloc(NULL,size) calls malloc(size)]
 */

 	if (d->prob < 1.0) {
		ErrMsg(ER_IMPOSVAL, "sample in R, not in gstat");
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

void free_data_gridmap(DATA_GRIDMAP *t) {
	efree(t->grid_base);
	efree(t->dpt);
	efree(t);
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
