#ifndef UTILS_H
#define UTILS_H

# include <stddef.h> /* size_t */

#ifdef NDEBUG
# define assert(x)
#else
# error "assert.h being included" /* remove this to activate assert */
# include <assert.h> /* assert() */
#endif
# include <stdio.h> /* FILE */
# include <string.h> /* FILE */

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

typedef struct {
	double d;
	int index;
} Double_index;

int esystem(char *cmd);
FILE *efopen(const char *filename, const char *mode);
FILE *etmpfile(void);
int   efclose(FILE *stream);
int   eremove(const char *name);
int   file_exists(const char *name);
char *string_prompt(const char *prompt);
char *string_file(const char *name);
#define string_dup(s) strcpy((char *)emalloc((strlen(s) + 1) * sizeof(char)), s)
int string_casecmp(const char *a, const char *b);
const char *string_cat(const char *s, const char *t);
size_t file_size(const char *);
char *get_line(char **s, int *size, FILE *stream);
/* char *ftoa(const char *fmt, float *f); */
char *my_dtoa(const char *fmt, double *d);
int almost_equals(const char *tok, const char *str);
void set_mv_float(float *f);
int is_mv_float(const float *f);
void set_mv_double(double *d);
int is_mv_double(const double *d);
int cpu_is_little_endian(void);
void print_file_record(void);
void elapsed(void);
void gstat_start(void);
char *store_argv(int argc, char *argv[]);
const char *save_string(const char *msg);
void save_strcat(STRING_BUFFER *dest, const char *src);
STRING_BUFFER *resize_strbuf(STRING_BUFFER *b, unsigned int size);
void free_strbuf(STRING_BUFFER *b);
int CDECL double_index_cmp(const Double_index *a, const Double_index *b);
int grass(void);
char *temp_name(void);

#ifdef HAVE_POPEN
FILE *epopen(const char *filename, const char *mode);
int   epclose(FILE *stream);
#endif

#ifdef MEMDEBUG
#	define emalloc(s) edmalloc(s, __FILE__ , __LINE__ )
#	define ecalloc(n, s) edcalloc(n, s, __FILE__, __LINE__)
#	define erealloc(p, s) edrealloc(p, s, __FILE__ , __LINE__ )
#	define efree(s) edfree(s, __FILE__ , __LINE__ )
#	define memtest() check_mem( __FILE__ , __LINE__ )
	void *edmalloc();
	void *edcalloc();
	void *edrealloc();
	void edfree();
#else  /* no MEMDEBUG defined: */
# ifdef DMALLOC
#	define efree free
#	define emalloc malloc
#	define ecalloc calloc
#	define erealloc realloc
# else
	void *emalloc(size_t size);
	void *ecalloc(size_t nobj, size_t size);
	void *erealloc(void *p, size_t size);
	void efree(void *p);
# endif
#	define memtest()  /* empty */
#endif

#if defined(__cplusplus)
}
#endif

#ifdef DMALLOC
# include <dmalloc.h>
#endif

#endif /* UTILS_H */
