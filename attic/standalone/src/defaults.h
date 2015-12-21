#ifndef DEFAULTS_H
#include <limits.h> /* INT_MAX */
#include <float.h> /* DBL_EPSILON */
# define DEFAULTS_H /* avoid multiple inclusion */

#define DEF_alpha           0.0
#define DEF_beta            0.0
#define DEF_bounds         NULL
#define DEF_cn_max         -1.0
#define DEF_coincide         -1
#define DEF_cressie           0
#define DEF_cutoff         -1.0
#define DEF_display    "display"
#define DEF_dots            500
#define DEF_fit               0
#define DEF_fit_limit    1.0E-5
#define DEF_format          "%g"
#define DEF_fraction    0.33333 /* fraction of max_dist for def. cutoff */
#define DEF_gauss             1
#define DEF_gcv             0.0
#define DEF_plotfile "gstatXXX.plt"
#ifndef WIN32
# define DEF_gnuplot   "gnuplot"
#else
# define DEF_gnuplot  "wgnuplot"
#endif
#define DEF_gpterm         NULL
#define DEF_idp             2.0
#define DEF_intervals        15 /* default number of intervals */
#define DEF_is_pdf            0 /* default to cdf indicator simulation */
#define DEF_iter             50
#define DEF_iwidth         -1.0
#define DEF_jgraph            0
#define DEF_lhs               0
#define DEF_longlat           0
#define DEF_n_marginals       0
#define DEF_nocheck           0 /* do check */
#define DEF_marginal_names NULL
#define DEF_marginal_values NULL
#define DEF_mv_string       "NaN"
#define DEF_nblockdiscr       4
#define DEF_n_uk        INT_MAX
#define DEF_numbers           1
#define DEF_nsim              1
#define DEF_ofilename      NULL
#define DEF_order             0
#define DEF_plotweights       0
#define DEF_pager         "more"
#define DEF_pairs             0
#define DEF_quantile        0.5
#define DEF_rowwise           1
#define DEF_rp                1
#define DEF_secure            0
#define DEF_seed              0
#define DEF_sim_beta          0
#define DEF_sparse            0
#define DEF_spiral            0
#define DEF_split             4
#define DEF_sym_ev            0
#define DEF_table_size        0
#define DEF_tol_hor       180.0
#define DEF_tol_ver       180.0
#define DEF_gls_residuals     0
#define DEF_xvalid            0
#define DEF_zero          (DBL_EPSILON * 10.0)
#define DEF_zero_est          0 /* ZERO_DEFAULT */
#define DEF_zmap            0.0
#define GNUFIT_NAME   "gnuplot.fit"

#endif /* DEFAULTS_H */
