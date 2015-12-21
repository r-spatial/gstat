typedef enum {
	GNUPLOT = 0,
	PSLATEX,
	CGM,
	GIF,
	EPS,
	PNG,
	EEPIC,
	UNKNOWN
} PLOT_TYPE;

int fprint_gnuplot_variogram(FILE *f, const VARIOGRAM *v, char *epsf,
	PLOT_TYPE p, int window_nr);
int fprint_jgraph_variogram(FILE *f, const VARIOGRAM *v);
void fprint_gnuplot_model(FILE *f, const VARIOGRAM *vgm, int fit);
