/* report.c */
void report_xvalid(double *xdata, double *xest, double *xdiff, 
		double *xstd, double *xzscore, int ndata, int var);

void write_points(const char *fname, DATA *d, DPOINT *where, double *est,
	int n_outfl);
