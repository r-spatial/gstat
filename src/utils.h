#ifndef UTILS_H
#define UTILS_H

# include <stddef.h> /* size_t */

/* some famous beware-of-side-effects macro's ! */
#ifndef MAX
# define MAX(a,b) (((a) > (b)) ? (a) : (b))
#endif
#ifndef MIN
# define MIN(a,b) (((a) < (b)) ? (a) : (b))
#endif
#ifndef ABS
#define ABS(a)   (((a) >= 0) ? (a) : (-(a)))
#endif
#ifndef SQR
# define SQR(a) ((a)*(a)) 
#endif
#ifndef PI
# define PI 3.14159265359
#endif
#define NULS(a) (a==NULL ? "<Null>" : a)

/*
 LTI: lower triangular matrix index, stored in an array:
 col | 0 1 2 3
 ----+-------
   0 | 0
 r 1 | 1 2
 o 2 | 3 4 5
 w 3 | 6 7 8 9

 row and col may be interchanged: LTI(a,b)==LTI(b,a)

 LTI2(a,b) is the index of an off-diagonal lower triangular matrix:
 col | 0 1 2 3
 ----+-------
   0 | x
 r 1 | 0 x
 o 2 | 1 2 x
 w 3 | 3 4 5 x
*/

#define LTI(r,c) ((r) >= (c) ? (((r)*(r+1))>>1)+(c) : (((c)*(c+1))>>1)+(r))
#define LTI2(r,c) ((r) >= (c) ? (((r)*(r-1))>>1)+(c) : (((c)*(c-1))>>1)+(r))
/* Note: `>>1' replaced `/2' to circumvent an hp 10.20 optimizer bug */

#if defined(__cplusplus)
extern "C" {
#endif

typedef struct {
	char *str;
	unsigned int max_length;
} STRING_BUFFER;

void set_mv_float(float *f);
void set_mv_double(double *d);
int is_mv_float(const float *f);
int is_mv_double(const double *d);
int almost_equals(const char *tok, const char *str);

void *emalloc(size_t size);
void *ecalloc(size_t nobj, size_t size);
void *erealloc(void *p, size_t size);
void efree(void *p);

#if defined(__cplusplus)
}
#endif

#endif /* UTILS_H */
