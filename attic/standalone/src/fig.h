typedef enum { 
	CONTIGUOUS, 
	BLOCKS 
} LEGEND_MODE;

void fig_setup(char *fname, double xmin, double ymin, double xmax, double ymax);
void fig_point(int color, double x, double y, int radius);
void fig_color(int color, float *col);
void fig_legend(char **entries, int n, LEGEND_MODE mode);
void fig_block(int color, int xmin, int ymin, int xmax, int ymax, int line);
void fig_end(void);
void fig_run(int color, double xmin, double ymin, double cellsizex, 
		double cellsizey, int rl);
