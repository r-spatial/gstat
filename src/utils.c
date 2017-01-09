/*
 * utils.c: error checking functions for file, memory and string handling
 */
#include <stdlib.h> /* free(), malloc() etc */
#include <ctype.h> /* tolower(), isspace() */
#include <string.h> /* strlen(), memcmp() */

#include "defs.h"
#include "userio.h"
#include "utils.h"
#include "glvars.h"
#include "debug.h"

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

void set_mv_float(float *f) {
	memset(f, 0xFF, sizeof(float));
}

void set_mv_double(double *d) {
	memset(d, 0xFF, sizeof(double));
}

int is_mv_float(const float *f) {
	const unsigned char u[sizeof(float)] = { 0xFF, 0xFF, 0xFF, 0xFF };
	/* will choke if sizeof(float) != 4 */
	return (memcmp(f, u, sizeof(float)) == 0);
}

int is_mv_double(const double *d) {
	const unsigned char u[sizeof(double)] =
		{ 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF };
	/* will choke if sizeof(double) != 8 */
	return (memcmp(d, &u, sizeof(double)) == 0);
}

/*
 * almost_equals() compares string value of token tok with str[], and
 *   returns TRUE if they are identical up to the first $ in str[].
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
