#ifndef POLYGON_H
#define POLYGON_H
POLYGON *read_polygons(const char *filename, int *n_polys, double **iso_values);
POLYGON read_n_points(FILE *f, int np);
int point_in_polygon(PLOT_POINT point, POLYGON *p);
void setup_poly_minmax(POLYGON *pl);
void read_edges(void);
void check_edges(DATA *d, const DPOINT *where);
void report_edges(void);
#endif
