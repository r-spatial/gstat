/* sem.c */
#ifndef SEM_H
# define SEM_H /* avoid multiple inclusion */

#if defined(__cplusplus)
extern "C" {
#endif

int calc_variogram(VARIOGRAM *v, const char *fname);
void fill_cutoff_width(DATA *data, VARIOGRAM *v);
int is_directional(VARIOGRAM *v);
void fprint_header_vgm(FILE *f, const DATA *d1, const DATA *d2, 
		const SAMPLE_VGM *ev);
void fprint_sample_vgm(FILE *f, const SAMPLE_VGM *ev);

#if defined(__cplusplus)
}
#endif

#define LONGSIZE (sizeof(unsigned long))
#define MAX_NH (1UL << (4 * LONGSIZE))
#define TO_NH(x,y) (x + (y << (4 * LONGSIZE)))
#define HIGH_NH(x) (x / (1UL << (4 * LONGSIZE)))
#define LOW_NH(x) (x % (1UL << (4 * LONGSIZE)))
#endif /* SEM_H */
