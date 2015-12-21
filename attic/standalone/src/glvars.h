#ifndef GLVARS_H
# define GLVARS_H /* avoid multiple inclusion */

typedef enum { 
	NSP = 0,  /* initial value */
	UIF,     /* variogram modelling user interface */
	OKR, UKR, SKR, /* ordinary, universal or simple kriging */
	IDW,  /* inverse distance interpolation */
	MED,  /* (local) sample median or quantile */
	NRS,  /* neighbourhood size */
	LSLM,  /* uncorrelated (or weighted) linear model */
	GSI, ISI, /* Gaussian/indicator (conditional) simulation */
	MAPVALUE,  /* mask map value at data location */
	SEM, COV,  /* sample (cross) semivariance or covariance */
	SPREAD, /* distance to nearest sample */
	XYP,  /* x and y coordinate of location */
	POLY, /* point-in-polygon */
	DIV, /* diversity, range */
	SKEW, /* skewness, kurtosis */
	LSEM, /* locally fitted semivariogram parameters */
	TEST  /* does nothing really */
} METHOD;

typedef struct {
	METHOD m; 
	int is_simulation;
	const char *name;
} METHODS;

extern const METHODS methods[];

typedef enum { 
	MODE_NSP = 0,
	SIMPLE,
	STRATIFY,
	MULTIVARIABLE
} MODE;

#if defined(__cplusplus)
extern "C" {
#endif

int init_global_variables(void);
const char *get_outfile_namei(int i);
const char **get_outfile_name(void);
int dump_all(void);
void check_global_variables(void);
const char *method_string(METHOD i);
int get_n_vars(void);
int get_n_vgms(void);
int get_n_outfile(void);
int get_n_beta_set(void);
int which_identifier(const char *id);
const char *name_identifier(int i);
const char *what_is_outfile(int i);
void push_mask_name(const char *name);
void push_bound(double value);
const char *get_mask_name(int i);
int get_n_masks(void);
void set_method(METHOD);
int is_simulation(METHOD m);
METHOD get_default_method(void);
METHOD get_method(void);
void set_mode(void);
MODE get_mode(void);
double max_block_dimension(int reset);
int n_variograms_set(void);
int decide_on_coincide(void);
int remove_id(const int id);
void remove_all(void);

#ifdef VARIO_H /* vario.h was included before this point: */
VARIOGRAM *get_vgm(int i);
#endif

#ifdef DATA_H /* data.h was included before this point: */
DATA **get_gstat_data(void);
DATA *get_dataval(void);
DATA *get_data_area(void);
DATA *create_data_area(void);
DPOINT *get_block_p(void);
void setup_valdata_X(DATA *d);
#endif

void push_edges_name(const char *name);
const char *get_edges_name(const int i);
int get_n_edges(void);
int *get_n_edges_polys(void);
int *set_n_edges_polys(int *n);

#ifdef DATA_H /* data.h was included before this point: */
POLYGON **get_edges(void);
POLYGON **set_edges(POLYGON **new_edges);
DATA *setup_poly_method(void);
#endif
    
        
#if defined(__cplusplus)
}
#endif

extern int gl_nblockdiscr, gl_seed, gl_n_uk, gl_dots, gl_cressie, gl_zero_est,
	gl_fit, gl_iter, gl_xvalid, gl_gauss, gl_sym_ev, gl_jgraph,
	gl_order, gl_n_intervals, gl_gls_residuals, gl_asym_vgm,
	gl_numbers, gl_nsim, gl_lhs, gl_longlat, gl_n_marginals, gl_sparse, gl_rp,
	gl_coincide, gl_nocheck, gl_spiral, gl_secure, gl_split,
	gl_register_pairs, gl_plotweights, gl_sim_beta, gl_rowwise;
extern double gl_rho, gl_idp, gl_cutoff, gl_iwidth, gl_zmap, gl_gcv,
	gl_quantile, gl_fit_limit, gl_fraction, gl_cn_max, gl_alpha,
	gl_beta, gl_tol_hor, gl_tol_ver, *gl_bounds,
	*gl_marginal_values, gl_zero, gl_zero2;
extern char *gl_delimiters, *gl_plotfile, *gl_format, *gl_gnuplot, 
	*gl_gnuplot35, *gl_mv_string, *gl_pager, *gl_lhsmean, *gl_lhsvar, 
	*gl_gpterm, *logfile_name, *command_file_name, *command_line, *o_filename, 
	*gl_display, *argv0;
extern char **gl_marginal_names; extern const char *method_code[];
extern FILE *plotfile;

#endif /* GLVARS_H */
