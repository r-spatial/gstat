/* random.c */
void set_rng_functions(double (*unif)(void), double (*norm)(void),
	const char *name);
void set_seed(unsigned long int i);
unsigned long int get_seed(void);
int e_random(int argc, char *argv[]);

double p_uniform(double v);
double q_uniform(double p);
double r_uniform(void);

double p_normal(double z);
double q_normal(double p);
double r_normal(void);

/*
double p_triangular(double z);
double q_triangular(double p);
double r_triangular(void);
*/
