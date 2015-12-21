void save_sim(DATA **data, DPOINT *where, int sim, int n_vars,
	const double *value, int *is_pt);
void save_sim_strat(DATA *d, DPOINT *where, int sim, double value, int is_pt);
void restore_data_sel(DATA **data, int sim, int n_vars);
void save_simulations_to_ascii(const char *fname);
void save_simulations_to_maps(GRIDMAP *mask);
void lhs(DATA **d, int n_vars, int stratify);
void init_simulations(DATA **d);
void set_beta(DATA **d, int sim, int n_vars, METHOD method);
void setup_beta(DATA **d, int n_vars, int n_sim);
void print_sim(void);
void free_simulations(void);
#ifndef SIM_DOUBLE
float ***get_msim(void);
#endif
