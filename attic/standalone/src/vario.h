#ifndef VARIO_H
#	define VARIO_H /* avoid multiple inclusion */

typedef enum {
	NOTSPECIFIED = 0,
	SEMIVARIOGRAM,
	CROSSVARIOGRAM,
	COVARIOGRAM,
	CROSSCOVARIOGRAM,
	PRSEMIVARIOGRAM /* pairwise relative semivariogram */
} SAMPLE_VGM_TYPE;

extern const char *vgm_type_str[];

#define LENGTH_OF_MODEL 100 /* max string length for one variogram model */

typedef enum {
	ZERO_DEFAULT = 0,
	ZERO_INCLUDE,
	ZERO_AVOID,
	ZERO_SPECIAL
} DO_AT_ZERO;

typedef struct {
	double tm[3][3]; /* 3D-transformation matrix */
	double angle[3]; /* angle in <x,y>, ccl from pos x; angle up; rot. angle */
	double ratio[2]; /* ratio axis2:axis1, ratio axis axis3:axis1 */
} ANIS_TM;

typedef enum { 
	NO_FIT = 0, 
	WLS_FIT, 
	WLS_FIT_MOD, 
	WLS_GNUFIT,
	WLS_GNUFIT_MOD, 
	MIVQUE_FIT,
	OLS_FIT,
	WLS_NHH
} FIT_TYPE;

typedef struct {
	int n_est, n_max, cloud, plot_numbers, is_asym;
	int recalc, refit, pseudo, is_directional;
	double *gamma, *dist;
	unsigned long *nh;
	double cutoff, iwidth;
	SAMPLE_VGM_TYPE evt;
	FIT_TYPE fit;
	DO_AT_ZERO zero;
	void *map, /* variogram map structure, i/o using files */
		*S_grid /* variogram map structure, passed from S interface*/ ;
	struct { double x, y, z; } direction;
	DPOINT ***pairs;
	/* optionally, the point pair list -- for j in [ 0, nh[i] >
	((DPOINT ***)pairs)[i][j*2] and its successor are two
	pointers to a pair of data points that were used to calculate gamma[i].
	The length of pairs is (at least) nh[i] * 2. See register_pairs() */
} SAMPLE_VGM; /* a sample variogram */

typedef enum {
	NOT_SP = 0, 
	NUGGET, 
	EXPONENTIAL,
	SPHERICAL, 
	GAUSSIAN, 
	EXCLASS,
#ifdef USING_R
	MATERN,
	STEIN,
#endif
	CIRCULAR,
	LINEAR, 
	BESSEL, 
	PENTASPHERICAL, 
	PERIODIC, 
	HOLE,
	LOGARITHMIC, 
	POWER, 
	SPLINE,
	LEGENDRE,
	MERROR, 
	INTERCEPT
} VGM_MODEL_TYPE;

typedef struct {
	VGM_MODEL_TYPE model;
	const char *name, *name_long, *v_gnuplot, *c_gnuplot;
	double (*fn)(double h, double *r), /* variogram value at h of basic model */
		(*da_fn)(double h, double *r); /* it's derivative to the range parm. */
} V_MODEL;
extern const V_MODEL v_models[];

#define NRANGEPARS 2 /* number of range parameters in variogram models */
typedef struct {
	VGM_MODEL_TYPE model;
	int fit_sill, fit_range, id;
	double range[NRANGEPARS], sill,
	(*fnct)(double h, double *r), /* (partial) unit variogram function */
	(*da_fnct)(double h, double *r); /* (partial) derivative to range of unit variogram */
	ANIS_TM *tm_range;
} VGM_MODEL;

typedef struct {
	long n; /* length */
	double maxdist, *values;
	ANIS_TM *tm_range;
} COV_TABLE;
#define COV_TABLE_VALUE(tablep, dist) \
	(dist >= tablep->maxdist ? tablep->values[tablep->n - 1] : \
	tablep->values[(int) floor(tablep->n * (dist / tablep->maxdist))])
#define SEM_TABLE_VALUE(tablep, dist) \
	(tablep->values[0] - COV_TABLE_VALUE(tablep, dist))

typedef struct {
	char 	*descr, *fname, *fname2; /* descript. and sample variogram (maps) */
	int 	n_models, max_n_models, n_fit, id, id1, id2,
			block_semivariance_set, block_covariance_set, isotropic,
			is_valid_covariance, fit_is_singular;
	VGM_MODEL *part;			/* the basic models */
	COV_TABLE *table;			/* covariance value table */
	double	block_semivariance,	/* average within-block semivariance */
			block_covariance,	/* average within-block covariance */
			max_range,	/* maximum range: where sill is reached */
			sum_sills,	/* sum of partial sill's */
			measurement_error, /* measurement error value--def. zero */
			max_val,	/* maximum value that is ever reached */
			min_val,	/* minimum value that is ever reached */
			SSErr;		/* fit result */
	SAMPLE_VGM *ev;
} VARIOGRAM;

#define dist2(x,y,z) (x*x+y*y+z*z)
#define TM_IS3D(tm) \
	(tm->angle[1] != 0.0 || tm->angle[2] != 0.0 || tm->ratio[1] < 1.0)
#define relative_norm(v,x,y,z) ((v == NULL || dist2(x,y,z) == 0.0) ? 1.0 : \
 (transform_norm(NULL,x,y,z)/transform_norm(v,x,y,z)))
#define UnitCovariance(part,x,y,z) \
 (part.model == INTERCEPT ? (1.0) :\
 (1.0 - part.fnct(transform_norm(part.tm_range,x,y,z), part.range)))
#define UnitSemivariance(part,x,y,z) \
 (part.fnct(transform_norm(part.tm_range,x,y,z), part.range))
#define Covariance(part,x,y,z) (part.sill * UnitCovariance(part,x,y,z))
#define Semivariance(part,x,y,z) (part.sill * UnitSemivariance(part,x,y,z))

#define DA_DELTA 0.001
#define da_Semivariance(part,x,y,z) \
 (part.da_fnct != NULL ? \
 (part.sill * part.da_fnct(transform_norm(part.tm_range,x,y,z), part.range)) : \
 da_general(&(part), transform_norm(part.tm_range, x,y,z)))

/*
 (part.sill * \
 (part.fnct(transform_norm(part.tm_range,x,y,z), part.range * (1.0 + DA_DELTA)) - \
 part.fnct(transform_norm(part.tm_range,x,y,z), part.range * (1.0 - DA_DELTA))) \
 / (2 * part.range * DA_DELTA)))
*/

#define EPSILON 1.0e-30 /* for small, non-zero anisotropy ratios */
#define is_covariogram(v) \
  ((v->ev->evt==COVARIOGRAM||v->ev->evt==CROSSCOVARIOGRAM))
#define is_variogram(v) \
  ((v->ev->evt==SEMIVARIOGRAM||v->ev->evt==CROSSVARIOGRAM))
#define is_direct(v) \
  ((v->ev->evt==COVARIOGRAM||v->ev->evt==SEMIVARIOGRAM))
#define is_cross(v) \
  ((v->ev->evt==CROSSCOVARIOGRAM||v->ev->evt==CROSSVARIOGRAM))

#define MODELHASNORANGE(m) (m == NUGGET || m == INTERCEPT)
#define PARTHASNORANGE(m) (m->model == NUGGET || m->model == INTERCEPT || \
	(m->model == LINEAR && m->range == 0.0))

#if defined(__cplusplus)
extern "C" {
#endif

extern const char *v_model_gnuplot[];
extern const char *v_model_names[];
extern const char *c_model_gnuplot[];

VARIOGRAM *init_variogram(VARIOGRAM *v);
SAMPLE_VGM *init_ev(void);
void vgm_init_block_values(VARIOGRAM *v);
void free_variogram(VARIOGRAM *v);
void logprint_variogram(const VARIOGRAM *v, int verbose);
void fprint_variogram(FILE *f, const VARIOGRAM *v, int verbose);
const char *sprint_variogram(const VARIOGRAM *v, int verbose);
double get_semivariance(const VARIOGRAM *v, double dx, double dy, double dz);
double get_covariance  (const VARIOGRAM *v, double dx, double dy, double dz);
double transform_norm(const ANIS_TM *tm, double dx, double dy, double dz);
void check_variography(const VARIOGRAM **v, int n);
void update_variogram(VARIOGRAM *vp);
double get_max_sills(int n);
double da_general(VGM_MODEL *p, double h);
double effective_range(const VARIOGRAM *v);
int push_variogram_model(VARIOGRAM *v, VGM_MODEL part);
VGM_MODEL_TYPE which_variogram_model(const char *m);
double relative_nugget(VARIOGRAM *v);
int vario(int argc, char **argv);
FIT_TYPE fit_int2enum(int fit);
DO_AT_ZERO zero_int2enum(int zero);
DO_AT_ZERO zero_shift(DO_AT_ZERO now, int next);
FIT_TYPE fit_shift(FIT_TYPE now, int next);
VGM_MODEL_TYPE model_shift(VGM_MODEL_TYPE now, int next);
int get_n_variogram_models(void);
void push_to_v(VARIOGRAM *v, const char *mod, double sill, double *range, 
		int nrangepars, double *d, int fit_sill, int fit_range);
void push_to_v_table(VARIOGRAM *v, double maxdist, int length, double *values,
		double *anis);

#if defined(__cplusplus)
}
#endif

#endif /* VARIO_H */
