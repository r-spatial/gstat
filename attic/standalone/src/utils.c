/*
    Gstat, a program for geostatistical modelling, prediction and simulation
    Copyright 1992, 2011 (C) Edzer Pebesma

    Edzer Pebesma, edzer.pebesma@uni-muenster.de
	Institute for Geoinformatics (ifgi), University of Münster 
	Weseler Straße 253, 48151 Münster, Germany. Phone: +49 251 
	8333081, Fax: +49 251 8339763  http://ifgi.uni-muenster.de 

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version. As a special exception, linking 
    this program with the Qt library is permitted.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

    (read also the files COPYING and Copyright)
*/

/*
 * utils.c: error checking functions for file, memory and string handling
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <ctype.h> /* tolower(), isspace() */
#include <math.h> /* floor() */
#include <string.h> /* strlen(), memcmp() */

#include "defs.h"
#ifdef HAVE_STAT_H
# include <sys/types.h> 
# include <sys/stat.h>
#endif /* HAVE_STAT_H */

#ifdef HAVE_UNISTD_H
# include <unistd.h>
#endif

#ifdef HAVE_LIBREADLINE
# include <readline/readline.h>
# include <readline/history.h>
#endif

#include "userio.h"
#ifdef HAVE_LIBCSF
# include "csftypes.h"
#endif
#include "utils.h"
#include "glvars.h"
#include "debug.h"

#ifdef HAVE_LIBGIS
# include <grass/gis.h>
# include <grass/gisdefs.h>
#endif 

#ifdef USING_R
# define NO_STD_IN_OUT
#endif

static void convert_null_to_space(char *cp, const char *name, const FILE *f);

typedef enum { IS_FILE, IS_PIPE } FILE_TYPE;
typedef enum { IS_OPEN, IS_CLOSED, IS_REMOVED } FILE_STATUS;

typedef struct {
	char *name, *mode;
	int nr, remove_at_exit;
	FILE_TYPE type;
	FILE_STATUS status;
	const FILE *f;
} FILE_RECORD;

static FILE_RECORD *file_record = NULL;
static int file_record_size = 0;
static time_t start;

static void record_open(const FILE *f, const char *name, const char *mode, 
	FILE_TYPE t);
static void record_closed(const FILE *f);
static void record_removed(const char *name);
static FILE_TYPE what_is_file(const FILE *f);
static const char *stream_name(const FILE *f);

#ifdef MEMDEBUG
/* #define RECORD_MALLOC */
/* #define RECORD_FREE  */
#endif

#ifndef MEMDEBUG
#ifndef DMALLOC
void efree(void *p) {
	if (p == NULL)
		pr_warning("efree(): NULL pointer as argument");
	else /* there's little point in calling free(NULL) */
		free(p);
}

void *emalloc(size_t size) {
	void *p = NULL;
	if (size == 0) {
		pr_warning("emalloc(): size 0 requested");
		return NULL;
	}
	p = (void *) malloc(size);
	if (p == NULL) {
		if (DEBUG_DUMP)
			message("malloc(%u) returned NULL", size);
		ErrMsg(ER_MEMORY, "");
	}
	return p;
}

void *ecalloc(size_t nobj, size_t size) {
	void *p = NULL;

	if (size == 0) {
		pr_warning("ecalloc(): size 0 requested");
		return NULL;
	}
	p = (void *) calloc(nobj, size);
	if (p == NULL) {
		if (DEBUG_DUMP)
			message("calloc(%u,%u) returned NULL", nobj, size);
		ErrMsg(ER_MEMORY, "");
	}
	return p;
}

void *erealloc(void *p, size_t size) {
	if (size == 0) {
		pr_warning("erealloc(): size 0 requested");
		return NULL;
	}
	if (p == NULL)
		p = (void *) malloc(size);
	else
		p = (void *) realloc(p, size);
	if (p == NULL) {
		if (DEBUG_DUMP)
			message("realloc(%u) returned NULL\n", size);
		ErrMsg(ER_MEMORY, "");
	}
	return p;
}
#endif

#else  /* MEMDEBUG defined: */
unsigned int n_mallocs = 0;

void print_n(void) {
	printf("# malloc's: %u\n", n_mallocs);
}

void edfree(void *p, char *file, int line) {
#ifdef RECORD_FREE
	printlog("%s:%d: free()\n", file, line);
#endif
	if (p == NULL)
		pr_warning("%s%s%s%d", "free(): NULL pointer as argument. File:", file,
				" Line: ", line);
	else /* little point in calling free(NULL) */
		free(p);
}

void *edmalloc(size_t size, char *file, int line) {
	void *p = NULL;
#ifdef RECORD_MALLOC
	printlog("%s:%d: malloc(%u)\n", file, line, size);
#endif
#ifndef USING_R
	if (n_mallocs == 0)
		atexit(print_n);
#endif
	n_mallocs++;
	p = (void *) malloc(size);
	if (p == NULL) {
		message("\n%s%s%s%d%s%u\n", "malloc(): out of memory in FILE: ", file,
				" LINE: ", line, " SIZE : ", size);
		ErrMsg(ER_MEMORY, "");
	}
	return p;
}

void *edcalloc(size_t nobj, size_t size, char *file, int line) {
	void *p = NULL;

#ifdef RECORD_MALLOC
	printlog("%s:%d: calloc(%u)\n", file, line, size);
#endif
	p = (void *) calloc(nobj, size);
	if (p == NULL) {
		message("\n%s%s%s%d%s%u\n", "calloc(): out of memory in FILE: ", file,
				" LINE: ", line, " SIZE : ", size);
		ErrMsg(ER_MEMORY, "");
	}
	return p;
}

void *edrealloc(void *p, size_t size, char *file, int line) {
#ifdef RECORD_MALLOC
	printlog("%s:%d: realloc(%u)\n", file, line, size);
#endif
	if (p == NULL)
		p = (void *) malloc(size);
	else
		p = (void *) realloc(p, size);
	if (p == NULL) {
		message("\n%s%s%s%d%s%u\n", "realloc(): out of memory in FILE: ", file,
				" LINE: ", line, " SIZE : ", size);
		ErrMsg(ER_MEMORY, "");
	}
	return p;
}

void check_mem(char *f, int *l) {
	char *p;

	p = (char *)emalloc(100);
	efree(p);
	message("check_mem FILE %s LINE %d\n", f, *l);
	return;
}
#endif /* MEMDEBUG */

FILE *efopen(const char *filename, const char *mode) {
/*
 * open file filename, warning with some diagnostics on error
 * return FILE *, error message on error.
 */
	FILE *tmp = NULL;
	int isread, error = 0;
#ifdef HAVE_STAT_H
	struct stat statbuf;
#endif /* HAVE_STAT_H */

	if (filename == NULL)
		ErrMsg(ER_NULL, "efopen()");
	if (filename[0] == '\0')
		ErrMsg(ER_NOFILE, "in function efopen()");
	if ((strchr(mode, '+')))
		ErrMsg(ER_IMPOSVAL, "efopen(): + mode not supported");
	isread = ((strchr(mode, 'r')) != NULL); /* read */

#ifndef NO_STD_IN_OUT
	if (strcmp(filename, "-") == 0) /* stdin/stdout */
		return isread ? stdin : stdout;
#endif

	switch (*filename) {
#ifdef HAVE_POPEN
		case '|': /* pipe */
			return (tmp = epopen(++filename, mode));
#endif
		case '>': /* append */
			if (isread) {
				pr_warning("file: %s", filename);
				ErrMsg(ER_IMPOSVAL, "efopen(): cannot read an append file");
			}
			while (isspace(*(++filename))) /* avoid '> file' to become ' file' */
				;
			mode = "a";
			/* BREAKTHROUGH: */
		default:
			tmp = fopen(filename, mode);
			if (tmp == NULL) {
				error = 1;
				if (! isread) /* try read-opening the thing, to fstat it */
					tmp = fopen(filename, "r");
			}
#ifdef HAVE_STAT_H
			if (tmp && fstat(fileno(tmp), &statbuf) < 0) {
				message("cannot fstat `%s'\n", filename);
				error = 1;
			}
			if (tmp && (statbuf.st_mode & S_IFMT) == S_IFDIR) {
				message("`%s' is a directory\n", filename);
				error = 1;
			}
#endif /* HAVE_STAT_H */
			break;
	}
	if (error == 1) {
		if (isread == 1) 
			ErrMsg(ER_READ, filename);
		else 
			ErrMsg(ER_WRITE, filename);
	}
	record_open(tmp, filename, mode, IS_FILE);
	return tmp;
}

int efclose(FILE *stream) {
/* 
 * close file stream, warning on error
 */
	int i;

#ifdef HAVE_POPEN
	if (what_is_file(stream) == IS_PIPE)
		i = pclose(stream);
	else
#endif
		i = fclose(stream);
	if (i == EOF)
		pr_warning("error on closing file");
	record_closed(stream);
	return i;
}

int esystem(char *cmd) {
	if (gl_secure) {
		pr_warning("prevented: system(\"%s\"):", cmd);
		ErrMsg(ER_SECURE, "");
	}
	return system(cmd);
}

FILE *etmpfile(void) {
	FILE *f;

	f = tmpfile();
	if (f == NULL)
		ErrMsg(ER_WRITE, "could not obtain tmpfile()");
	return f;
}

#ifdef HAVE_POPEN
FILE *epopen(const char *filename, const char *mode) {
	FILE *tmp = NULL;
	int isread;

	if (gl_secure) {
		pr_warning("prevented: popen(\"%s\"):", filename);
		ErrMsg(ER_SECURE, "");
	}
	isread = ((strchr(mode, 'r')) != NULL);
	tmp = (FILE *) popen((char *) filename, (char *) mode);
	record_open(tmp, filename, mode, IS_PIPE);
	if (tmp == NULL) {
		if (isread) 
			ErrMsg(ER_PREAD, filename);
		else 
			ErrMsg(ER_PWRITE, filename);
	}
	return tmp;
}

int epclose(FILE *stream) {
	int i = 0;

	i = pclose(stream);
	record_closed(stream);
	if (i == -1)
		pr_warning("pclose() on invalid stream");
	return i;
}
#endif

int eremove(const char *name) {
/*
 * remove file name, warning on error
 */
	int i;

	if (gl_secure) {
		pr_warning("secure mode prevented remove(\"%s\")", name);
		return 0;
	} 
	i = remove(name);
	record_removed(name);
	if (i != 0)
		pr_warning("error on removing file `%s'", name);
	return i;
}

int file_exists(const char *name) {
	FILE *f = NULL;

	if ((f = fopen(name, "r")) != NULL) {
		fclose(f);
		return 1;
	} else
		return 0;
}

#ifndef USING_R
char *get_line(char **s, int *size, FILE *stream) {
/* 
 * read line in *s, return number of chars read;
 * resize s and adjust *size if neccesary;
 * PRE: *s is a char *, pointing to NULL or dynamically allocated memory
 * return NULL on EOF and empty string;
 * after last line read.
 */
#define INCR 64
	int c;
	char cr;
	int n = 0;

	if (s == NULL || size == NULL || stream == NULL)
		ErrMsg(ER_NULL, "get_line()");
	if (*size == 0 || *s == (char *) NULL) {
		*s = (char *) emalloc(INCR * sizeof(char));
		*size = INCR;
	}
	while ((c = fgetc(stream)) != EOF) {
		cr = c;
		convert_null_to_space(&cr, NULL, stream);
		/* printf("char:[%c],int[%d]\n", c, c); */
		(*s)[n] = c;
		n++;
		if (n == *size - 1) { /* resize: leave space for '\0' */
			*size += INCR;
			*s = erealloc(*s, *size);
		}
		if (c == '\n') { /* end-of-line */
			(*s)[n] = '\0'; /* terminate string */
			return *s;
		}
	}
	/* at EOF: */
	(*s)[n] = '\0';
	if (n > 0) /* we've had character(s): */
		return *s;
	return NULL; /* EOF */
}

char *string_prompt(const char *prompt) {
	char *buf = NULL, *line = NULL;
	int buf_size = 4096, line_size = 0, i = 0;

	buf = (char *) emalloc(buf_size);
	buf[0] = '\0';
	printf("Enter commands, end with `e' or EOF\n");
	do {
		if (line != NULL)
			efree(line);
#ifdef HAVE_LIBREADLINE
		if ((line = readline(prompt)) != NULL && strlen(line))
			add_history(line);
#else
		line = NULL;
		line_size = 0;
		fprintf(stdout, "%s", prompt);
		line = get_line(&line, &line_size, stdin);
#endif /* else HAVE_LIBREADLINE */
		if (almost_equals(line, "e$\n") || almost_equals(line, "q$\n")) {
			efree(line);
			line = NULL;
		}
		if (line && strlen(line)) { /* non-empty string: add to buf */
			i += strlen(line) + 1; /* + trailing \n, from readline() */
			if (i >= buf_size - 1)
				buf = (char *) erealloc(buf, buf_size *= 2);
			strcat(buf, line);
#ifdef HAVE_LIBREADLINE
			strcat(buf, "\n");
#endif
		}
	} while (line);
	return buf;
}
#endif

char *string_file(const char *fname) {
/*
 * read file in dynamically allocated character string
 */
	FILE *in = NULL;
	char *buf = NULL, cr;
	int c;
	int buf_size = 1000, i = 0;

	/* read as ascii file */
	in = (fname == NULL) ? (FILE *) stdin : efopen(fname, "r");
	buf = (char *) emalloc(buf_size * sizeof(char));
	while ((c = fgetc(in)) != EOF) {
		cr = c;
		convert_null_to_space(&cr, fname ? fname : "stdin", NULL);
		if (i == buf_size) {
			buf_size += 1000;
			buf = (char *) erealloc(buf, buf_size * sizeof(char));
		}
		buf[i] = c;
		i++;
	}
	buf[i] = '\0'; /* close string */
	if (fname != NULL)
		efclose(in); /* close file */
	/* free unnecesary memory: */
	return (char *) erealloc(buf, (i + 1) * sizeof(char));
}

static void convert_null_to_space(char *cp, const char *fname, const FILE *stream) {
	static const char *fn = NULL;

	/* convert null characters */
	if (*cp == '\0') {
		*cp = ' ';
		if (fname == NULL)
			fname = stream_name(stream);
		if (fn != fname) { /* print only once: */
			pr_warning("converted null-character(s) in `%s' to space", fname);
			fn = fname;
		}
	}
	return;
}

int string_casecmp(const char *a, const char *b) {
	/* after strcmp(), K&RII p. 106: */
	int i;

	for (i = 0; tolower(a[i]) == tolower(b[i]); i++)
		if (a[i] == '\0')
			return 0;
	return tolower(a[i]) - tolower(b[i]);
}

const char *string_cat(const char *s, const char *t) {
/*
 * beware of strtok-like side effect: each second call overwrites the
 * first' call return value.
 */
	char *cp = NULL;

	if (cp == NULL)
		cp = emalloc((strlen(s) + strlen(t) + 1) * sizeof(char));
	else
		cp = erealloc(cp, (strlen(s) + strlen(t) + 1) * sizeof(char));
	strcpy(cp, s);
	strcat(cp, t);
	return cp;
}

size_t file_size(const char *fname) {
	FILE *f = NULL;
	size_t size;

	f = fopen(fname, "rb");
	if (f == NULL)
		return 0;
	fseek(f, 0L, SEEK_END); /* jump to end-of-file */
	size = ftell(f); /* get size of f at end-of-file */
	fclose(f);
	return size;
}

char *ftoa(const char *fmt, float *a) {
/*
 * BEWARE of the sideffect:
 * NEVER use printf("%10s %10s", ftoa("%g", 1.0), ftoa("%g", 2.0));
 * instead: printf("%10s", ftoa("%g", 1.0)); printf(" %10s", ftoa("%g", 2.0));
 */
    static char *s = NULL;

	if (s == NULL) /* first time: */
		s = (char *) emalloc(MAX(50, 1 + strlen(gl_mv_string)));
    s[0] = '\0';
    if (is_mv_float(a))
        sprintf(s, "%s", gl_mv_string);
    else
        sprintf(s, fmt, *a);
    return s;
}

char *my_dtoa(const char *fmt, double *a) {
/*
 * BEWARE of the sideffect:
 * NEVER use printf("%10s %10s", my_dtoa("%g", 1.0), my_dtoa("%g", 2.0));
 * instead: printf("%10s", my_dtoa("%g", 1.0)); printf(" %10s", my_dtoa("%g", 2.0));
 */
    static char *s = NULL;

	if (s == NULL) /* first time: */
		s = (char *) emalloc(MAX(50, 1 + strlen(gl_mv_string)));
    s[0] = '\0';
    if (is_mv_double(a))
        sprintf(s, "%s", gl_mv_string);
    else
        sprintf(s, fmt, *a);
    return s;
}

/*
 * almost_equals() compares string value of token tok with str[], and
 *   returns TRUE if they are identical up to the first $ in str[].
 * (admitted, this was stolen from gnuplot)
 */
int almost_equals(const char *tok, const char *str) {
	int i, after = 0, start = 0, len;

	if (tok == NULL) 
		return 0;	/* must be a value--can't be equal */
	len = strlen(tok);
	for (i = 0; i < len + after; i++) {
		if (str[i] != tok[start + i]) {
			if (str[i] != '$')
				return 0;
			else {
				after = 1;
				start--;
			}
		}
	}
	/* i now beyond end of token string */
	return(after || str[i] == '$' || str[i] == '\0');
}

void set_mv_float(float *f) {
#ifdef HAVE_LIBCSF /* csftypes.h was included: */
	SET_MV_REAL4(f);
#else
	memset(f, 0xFF, sizeof(float));
#endif
}

void set_mv_double(double *d) {
#ifdef HAVE_LIBCSF
	SET_MV_REAL8(d);
#else
	memset(d, 0xFF, sizeof(double));
#endif
}

int is_mv_float(const float *f) {
#ifdef HAVE_LIBCSF
	return IS_MV_REAL4(f);
#else
	const unsigned char u[sizeof(float)] = { 0xFF, 0xFF, 0xFF, 0xFF };
	/* will choke if sizeof(float) != 4 */
	return (memcmp(f, u, sizeof(float)) == 0);
#endif
}

int is_mv_double(const double *d) {
#ifdef HAVE_LIBCSF
	return IS_MV_REAL8(d);
#else
	const unsigned char u[sizeof(double)] =
		{ 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF };
	/* will choke if sizeof(double) != 8 */
	return (memcmp(d, &u, sizeof(double)) == 0);
#endif
}

int cpu_is_little_endian(void) {
/*
 * returns 0 if the current cpu is not little-endian,
 * returns 1 if it is. NOTE: VMS-order is not evaluated (don't know it!)
 */
#ifdef HAVE_CONFIG_H /* configure did the work... */
# ifdef WORDS_BIGENDIAN 
	return 0;
# else
	return 1;
# endif
#else /* do it self: */
	unsigned long u = 1;
	char *cp;

	cp = (char *) (&u);
	return (cp[0] == 1); /* are we little-endian? */
#endif
}

static void record_open(const FILE *f, const char *name, const char *mode,
		FILE_TYPE t) {
	int i;

	for (i = 0; i < file_record_size; i++) {
		if (strcmp(file_record[i].name, name) == 0) {
			if (file_record[i].status != IS_OPEN) /* file can be read twice */
				break; /* take this place for recording only if not open */
		}
	}
	if (i == file_record_size) { /* no matches: increase list */
		file_record_size += 1;
		file_record = (FILE_RECORD *) erealloc(file_record,
			file_record_size * sizeof(FILE_RECORD));
	} else {
		efree(file_record[i].name);
		efree(file_record[i].mode);
	}
	file_record[i].f = f;
	file_record[i].name = string_dup(name);
	file_record[i].mode = string_dup(mode);
	file_record[i].type = t;
	file_record[i].status = IS_OPEN;
	file_record[i].nr = i;
}

static void record_closed(const FILE *f) {
	int i;

	for (i = 0; i < file_record_size; i++)
		if (file_record[i].f == f) {
			if (file_record[i].status == IS_CLOSED)
				pr_warning("file %s closed twice", file_record[i].name);
			file_record[i].status = IS_CLOSED;
			file_record[i].f = NULL;
		}
}

static void record_removed(const char *name) {
	int i;

	for (i = 0; i < file_record_size; i++) {
		if (strcmp(file_record[i].name, name) == 0) {
			file_record[i].status = IS_REMOVED;
			file_record[i].f = NULL;
		}
	}
}

void print_file_record(void) {
	int i, open;

	for (i = open = 0; i < file_record_size; i++) {
		printlog("%d: %s `%s' (mode %s) ",
			file_record[i].nr,
			file_record[i].type == IS_PIPE ? "pipe" : "file",
			file_record[i].name,
			file_record[i].mode);
		switch (file_record[i].status) {
			case IS_OPEN:
				printlog("is open\n");
				break;
			case IS_CLOSED:
				printlog("was closed\n");
				break;
			case IS_REMOVED:
				printlog("was removed\n");
				break;
		}
	}
}

static FILE_TYPE what_is_file(const FILE *f) {
	int i;
	for (i = 0; i < file_record_size; i++)
		if (file_record[i].f == f)
			return file_record[i].type;
	assert(0);
	return IS_OPEN; /* never reached */
}

static const char *stream_name(const FILE *f) {
	int i;
	for (i = 0; i < file_record_size; i++)
		if (file_record[i].f == f)
			return file_record[i].name;
	assert(0);
	return "bogus"; /* never reached */
}

void gstat_start(void) {
	start = time(NULL);
}

void elapsed(void) {
	int hrs, mns, sec;
	double diff;

	diff = difftime(time(NULL), start); /* difference in seconds */
	if (diff < 10)
		return;
	hrs = floor(diff/3600.0);
	mns = floor((diff - hrs * 3600.0)/60.0);
	sec = floor(diff - hrs * 3600.0 - mns * 60.0);
	if (hrs == 0)
		printlog("elapsed time %d:%02d\n", mns, sec);
	else
		printlog("elapsed time %d:%02d:%02d\n", hrs, mns, sec);
	return;
}

char *store_argv(int argc, char *argv[]) {
	int i, len = 0;
	char *cp;

	for (i = 0, len = argc; i < argc; i++)
		len += strlen(argv[i]);
	cp = (char *) emalloc(len * sizeof(char));
	for (i = 0, cp[0] = '\0'; i < argc; i++) {
		strcat(cp, argv[i]);
		if (i < argc - 1)
			strcat(cp, " ");
	}
	return cp;
}

const char *save_string(const char *msg) {
#define MAX_SIZE (ERROR_BUFFER_SIZE/2)
	static char *s, *empty = "";

	if (msg == NULL)
		return empty;
	if (strlen(msg) > MAX_SIZE) { /* will usually never happen... */
		s = (char *) emalloc(MAX_SIZE * sizeof(char));
		strncpy(s, msg, MAX_SIZE-5);
		s[MAX_SIZE-5] = '\0';

		strcat(s, "...");
		return s;
	}
	return msg;
}

STRING_BUFFER *resize_strbuf(STRING_BUFFER *b, unsigned int size) {
	if (b == NULL) {
		b = (STRING_BUFFER *) emalloc(sizeof(STRING_BUFFER));
		b->str = (char *) emalloc(size * sizeof(char));
		b->str[0] = '\0';
	} else
		b->str = (char *) erealloc(b->str, size * sizeof(char));
	b->max_length = size;
	return b;
}

void free_strbuf(STRING_BUFFER *b) {
	if (b == NULL)
		return;
	efree(b->str);
	efree(b);
}

void save_strcat(STRING_BUFFER *dest, const char *src) {
	int len;
	
	assert(dest != NULL);
	assert(src != NULL);
	len = strlen(src) + 1;
	if (strlen(dest->str) + strlen(src) > dest->max_length)
		resize_strbuf(dest, dest->max_length + MAX(len,ERROR_BUFFER_SIZE));
	dest->str = strcat(dest->str, src);
}

int CDECL double_index_cmp(const Double_index *a, const Double_index *b) {
/* ANSI-qsort() conformant Double_index comparison function: sort on field d */
	if (a->d < b->d)
		return -1;
	if (a->d > b->d)
		return 1;
	return 0;
}

int grass(void) {
	static int gisinit = 0;
	int env, lock;
	char *str, *home /* , *gisrc */ ;

	if (gisinit == 1) /* been here before... */
		return gisinit;

	env = ((getenv("LOCATION") || 
		getenv("LOCATION_NAME"))&& getenv("GISDBASE") && getenv("MAPSET") ||
		getenv("GISRC"));
	if ((home = getenv("HOME")) == NULL)
		home = "";
	str = (char *) emalloc(strlen(home) + 20);
	str[0] = '\0';
	strcat(str, home);
	strcat(str, "/.gislock5");
	lock = file_exists(str);
	if (env || lock) {
#ifdef HAVE_LIBGIS
		if (gisinit == 0) {
			G_gisinit("gstat");
			/*
			gisrc = G__get_gisrc_file();
			printf("gisrc: [%s]\n", gisrc ? gisrc : "NULL");
			*/
			gisinit = 1;
		}
#else
		pr_warning("this version of gstat was not compiled with grass support");
#endif
	}
	efree(str);
	return gisinit;
}
