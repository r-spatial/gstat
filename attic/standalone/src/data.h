#ifndef DATA_H
#	define DATA_H /* avoid multiple inclusion */
#include <limits.h> /* INT_MAX */
#include <float.h> /* FLT_MAX */

#define ID_OF_VALDATA	INT_MAX
#define ID_OF_AREA		(INT_MAX-1)
#define IS_GLOBAL(d) (d->sel_rad >= DBL_MAX && d->sel_max >= INT_MAX)
#define DELIMITERS		" \t,\n\r"

#define N_POLY     18 
#define POLY_MIN (-19) /* lowest index */

#define POLY_X   (-19)
#define POLY_Y   (-18)
#define POLY_Z   (-17)
#define POLY_X2  (-16)
#define POLY_Y2  (-15)
#define POLY_Z2  (-14)
#define POLY_XY  (-13)
#define POLY_XZ  (-12)
#define POLY_YZ  (-11)
#define POLY_X3  (-10)
#define POLY_Y3   (-9)
#define POLY_Z3   (-8)
#define POLY_X2Y  (-7)
#define POLY_XY2  (-6)
#define POLY_X2Z  (-5)
#define POLY_XZ2  (-4)
#define POLY_Y2Z  (-3)
#define POLY_YZ2  (-2) /* -1 is reserved for "no intercept" */

typedef struct {
	int poly_nr;
	char *name;
	int degree, mode;
} POLY_NM;

typedef struct { 		/* structure to hold one point value: */
	double x, y, z, 	/* x, y and z coordinate */
		variance, 		/* the attribute's variance */
		attr; 			/* attribut value (data value) */
	union {
		float dist2;	/* squared distance to estimate point */
		float weight;	/* weight in block discretization */
		int stratum;	/* stratum of current point in data() list */
	} u;
	double	*X;			/* row entry in X matrix for this DPOINT */
	unsigned int bitfield;
		/* most right bit: IS_POINT (0) or IS_BLOCK (1),
		remaining left bits: index of this point in d->list */
} DPOINT;

/* qtree_search structs (nsearch.c): */
typedef struct {	/* defining a rectangular bounding box */
	double 	x, y, z, size;
	int mode;
} BBOX;

typedef struct qnode {	/* the struct used to define nodes in the search tree */
	int n_node;			/* >= 0: number of data points in this node */
					 	/* negative (-1) if u is a node list */
	union {
		struct qnode **node;/* pointers to 4 or 8 other nodes */			
		DPOINT **list; 	    /* or pointers to data points within this leaf */
	} u;
	BBOX bb;
} QTREE_NODE;

typedef struct {
	double x_ul, y_ul, cellsizex, cellsizey;
	unsigned int rows, cols;
	DPOINT ***dpt,   /* 2d array to (DPOINT *) entries in list */
		**grid_base; /* base of blocked memory allocation */
} DATA_GRIDMAP; /* the management summary */

/* polygon structs: */
typedef struct {
	double		x, y;
} PLOT_POINT;

typedef struct {
	PLOT_POINT	min, max;
} MBR;

typedef struct polygon {
	MBR mbr;
	int lines;
	PLOT_POINT	*p;
    int close; /* 1 - is closed polygon */
} POLYGON;

typedef enum {
	DATA_UNKNOWN = 0,
	DATA_ASCII_TABLE,	/* ascii table */
	DATA_EAS,			/* simplified GeoEAS format */
	DATA_IDRISI_VEC,	/* idrisi .vec */
	DATA_IDRISI32_VEC,	/* idrisi .vct */
	DATA_IDRISI_BIN,	/* Idrisi .img binary */
	DATA_IDRISI_ASCII,	/* Idrisi .img ascii */
	DATA_IDRISI32_BIN,	/* Idrisi32 .rst binary */
	DATA_IDRISI32_ASCII,/* Idrisi32 .rst ascii */
	DATA_GRIDASCII,		/* ArcInfo */
	DATA_GRIDFLOAT,		/* ArcInfo */
	DATA_CSF,			/* PCRaster */
	DATA_T2,			/* Mike-SHE grid format */
	DATA_ERMAPPER,		/* ER-Mapper raster dataset */
	DATA_GNUPLOT,		/* gnuplot binary grid */
	DATA_GMT,			/* GMT netCDF format */
	DATA_SURFER_DSAA,	/* Surfer DSAA ascii grid */
	DATA_GSLIB,			/* GSLIB ascii grid */
	DATA_GRASS,			/* GRASS site list */
	DATA_GRASS_GRID,	/* GRASS raster */
	DATA_GDAL, 			/* GDAL raster */
	DATA_EXT_DBASE  /* CW external database */
} DATA_TYPE_;

typedef struct {
	DATA_TYPE_ type;
	const char *name;
} DATA_TYPE;

extern const DATA_TYPE data_types[];

typedef struct {
	int to_var,  /* merge from current data to this variable */
	col_this_X,  /* merge this column number */
	col_other_X; /* to this column number in the other variable */
} MERGE_TABLE;

typedef struct {
	int size, max_size;
	double *val;
} D_VECTOR;

D_VECTOR *push_d_vector(double d, D_VECTOR *v);
void free_d_vector(D_VECTOR *v);

/* CW added, FTTB copies of DATA members
 * SEARCH_CRITERIA should become part of DATA
 * Dit zijn degene die ik begrijp
 */
typedef struct {
 /* DATA::prob sample later */
 int
		force, 			/* force neighbourhood selection */
		sel_min,
		sel_max,		/* min and max number for neighbourhood selection */
		oct_max,		/* max # pts for each octant; 0: use no octant search */
		oct_filled, /* RETURN VALUE? number of non-empty octants in selection */
		square;     /* use square search neighbourhood, default circular */
	double sel_rad;		/* radius for neighbourhhood selection */
}SEARCH_CRITERIA;

typedef struct {		/* structure that holds data info and lists */
	char *variable,		/* attr name, log(..) */
		*x_coord,		/* name of x coordinate */
		*y_coord,		/* name of y coordinate */
		*z_coord,		/* name of z coordinate */
		*s_coord,   	/* name of stratum variable */
		*V_coord,		/* name of variance */
		*Category,		/* category value, for indicator transform */
		*id_name,		/* name of ID column*/
		*fname,			/* file name */
		**point_ids,	/* IDs of points */
		*var_fn_str, 	/* VarFunction string */
		*nscore_table;	/* normal score table output file name */

	DATA_TYPE type;     /* what is this file? */
	int	id,				/* id of data, number in command file */
		n_list,			/* # points in list */
		n_original,     /* # real data (read from file, not simulated) */
		n_sel,			/* # of points in selection: sel */
		n_max,			/* maximum # in list */
		nsim_at_data,   /* nr. of pts at simulation locations */
		init_max,		/* user-specified maximum n (to save memory) */
		n_sel_max,		/* maximum number in sel */
		n_X, *colX,		/* number and columns in X matrix */
		log, 			/* is attr.value log-transformed ? */
		force, 			/* force neighbourhood selection */
		vdist,			/* use variogram value as distance crit. */
		n_averaged,		/* number of averaged data so far */
		colnx, colny, colnz, /* column-numbers */
		colnvariance,	/* column that holds variance */
		colnvalue,		/* column that holds attribute values */
		colns,          /* column that holds u (if strata) */
		coln_id,        /* column with ID */
		sel_min,
		sel_max,		/* min and max number for neighbourhood selection */
		oct_max,		/* max # pts for each octant; 0: use no octant search */
		oct_filled,     /* number of non-empty octants in selection */
		mode,			/* mode: 1(x),2(y),3(xy),4(z),5(xz).. */
		dummy,			/* is this variable a dummy variable? */
		standard,		/* if standard: data are standardized by 
							dividing all values by data->std */
		calc_residuals,	/* 1: do calculate OLS residuals for vgm est */
		is_residual,	/* attr values are residuals lm X */
		polynomial_degree, /* degree of coordinate polynomial */
		togrid,         /* shift data values to grid centres */
		square,         /* use square search neighbourhood */
		centre,			/* centre area + points() at area centre? */
		region,			/* select data in region? */
		average,		/* average measurements at identical location? */
		every,			/* sample only every x-th observation */
		offset, skip,	/* starting at offset */
		datatype,       
			/*KS added what idrisi data type 0=integer,1=real 1/18/99*/
        filetype;       
			/*KS added what idrisi file type 0=ascii,1=binary 1/18/99*/
	enum { U_UNKNOWN, U_ISDIST, U_ISWEIGHT, U_ISSTRATUM } what_is_u;
	double sel_rad,		/* radius for neighbourhhood selection */
		Icutoff,		/* cutoff value for indicator variable: */
						/* colnvalue = ("real value" <= Icutoff) */
		minX, maxX, minY, maxY, minZ, maxZ, /* min/max of coordinates */
		minvariance, maxvariance, /* min/max variance */
		mv,				/* missing value */
		dX,				/* max. vector norm X-space distance */
		prob,			/* inclusion probability (to sample data file) */
		lambda;			/* lambda value for box-cox transform */
	int  minstratum, maxstratum; /* min/max stratum */
	double mean, std;	/* sample mean and st.dev. of attribute */

/* CW members to hold data in this struct (DATA_TYPE!= DATA_EXT_DBASE)
 */
	DPOINT **list;		/* list of data points, of length n_list */
	DPOINT *P_base;		/* base for pointer array, if allocated blockwise */
#ifdef HAVE_EXT_DBASE
/* CW members when data is held external (DATA_TYPE== DATA_EXT_DBASE) */
	void   *ext_dbase;  /* ptr to EXTDBASE_LINK struct */
/* CW end of edits */
#endif

	DPOINT **sel;			/* list of selection indices, of length n_sel */

	double (*point_norm)(const DPOINT *); /* eucl. vector length */
	double (*pp_norm2)(const DPOINT *, const DPOINT *); /* point-point squared distance */
	double (*pb_norm2)(const DPOINT *, BBOX); /* point-BBOX distance: nsearch.c */
	double (*variance_fn)(double mu); /* variance function */
	double *X_base;		/* base pointer for X arrays, when allocated blockwise */
	void *lm,			/* cast to LM *, see lm.h */
		*glm;			/* remember: several matrices/vecs needed in gls.c */
	/* next 2 entries are for merging regressors across variables */
	/* to avoid double references, each entry var should be less than ->id */
	int n_merge;        /* merge_table size */
	MERGE_TABLE *mtbl;  /* entries in merge table */
	QTREE_NODE *qtree_root; /* a tree-based structure with pointers to list[]
						for fast neighbourhood search */
/* 	POLYGON *poly;	 */	/* where the polygons go */
#ifdef WITH_SPIRAL
	SPIRAL *spiral;     /* spiral search index structure */
#endif
	DATA_GRIDMAP *grid; /* grid map topology if data was read from a map */
	D_VECTOR *beta;
} DATA;

#define X_BIT_SET  1 /* also used in nsearch.c */
#define Y_BIT_SET  2
#define Z_BIT_SET  4
#define V_BIT_SET  8
#define S_BIT_SET 16

#define D_HAS_WEIGHT(d) (d->what_is_u == U_ISWEIGHT)
#define D_HAS_DIST(d)   (d->what_is_u == U_ISDIST)
#define D_HAS_STRATA(d) (d->colnu && d->what_is_u == U_ISSTRATA)

/* following routines do not depend on sizeof(int) (K&R II, p. 48, 49) */
#define left_bits(x) ((x) >> 1) /* shift one, zero most left bit */
#define right_bit(x) ((x) & ~(~0 << 1)) /* zero all except right bit */
#define set_right_bit_on(x) (x = (x) | 1)
#define set_right_bit_off(x) (x = (x) & (~0 << 1))
#define set_left_bits(x,val) (x = (val << 1) | right_bit(x))

#define GET_INDEX(p)     (left_bits((p)->bitfield))
#define SET_INDEX(p,val) (set_left_bits((p)->bitfield,val))
#define IS_BLOCK(p)      (right_bit((p)->bitfield))
#define IS_POINT(p)      (!right_bit((p)->bitfield))
#define SET_BLOCK(p)     (set_right_bit_on((p)->bitfield))
#define SET_POINT(p)     (set_right_bit_off((p)->bitfield))

#if defined(__cplusplus)
extern "C" {
#endif

DATA *read_gstat_data(DATA *d);
DATA *get_area_centre(DATA *area, DATA *valdata);
void centre_area(DATA *area);
void push_point(DATA *d, const DPOINT *p);
void pop_point(DATA *d, int list_nr);
void free_data(DATA *tmp);
void report_data(const DATA *d);
#define print_data_list(d) print_data(d, 1)
#define print_data_selection(d) print_data(d, 0)
void print_data(const DATA *d, int list);
void logprint_point(const DPOINT *p, const DATA *d);
DATA *init_one_data(DATA *data);
int coordinates_are_equal(const DATA *a, const DATA *b);
void init_data_minmax(void);
void setup_data_minmax(DATA *d);
void calc_polynomials(DATA *d);
double calc_polynomial(DPOINT *p, int colX);
char *print_data_line(const DATA *d, char **to);
extern const POLY_NM polynomial[N_POLY];
#define POLY_NAME(i) polynomial[i - POLY_MIN].name
#define POLY_DEGREE(i) polynomial[i - POLY_MIN].degree
void data_add_X(DATA *d, int i);
int push_to_merge_table(DATA *d, int to_var, int col_this_X, int col_other_X);
DATA_GRIDMAP *gsetup_gridmap(double x_ul, double y_ul, double cellsizex, 
			double cellsizey, unsigned int rows, unsigned int cols);
void datagrid_rebuild(DATA *d, int adjust_to_gridcentres);
void set_norm_fns(DATA *d);
double data_block_diagonal(DATA *data);
int intercept_only(const DATA *d);
double v_mu(double mu);
double v_mu2(double mu);
double v_mu3(double mu);
double v_bin(double mu);
double v_identity(double mu);
void setup_polynomial_X(DATA *d);
void calc_polynomial_point(DATA *d, DPOINT *pt);
double pp_norm_gc(const DPOINT *a, const DPOINT *b);

#if defined(__cplusplus)
}
#endif

#endif /* DATA_H */
