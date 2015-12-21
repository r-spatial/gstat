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
 * read.c: read int or float from a string, missing value and error handling
 */
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <string.h>

#include "defs.h"
#include "utils.h"
#include "userio.h"
#include "glvars.h"
#include "read.h"

/*
 * functions for save reading and error checking
 * on reading int, long, float and double
 */

#ifndef INT_MAX
#	define INT_MAX 32767
#	define INT_MIN -32767
#endif
#ifndef LONG_MAX
#	define LONG_MAX +2147483647L
#	define LONG_MIN -2147483647L
#endif
#ifndef FLT_MAX
#	define FLT_MAX 1.0E+37F
#	define FLT_MIN 1.0E-37F
#endif
#ifndef HUGE_VAL
#  define HUGE_VAL    1.7976931348623157e+308
#endif

int read_float(const char *s, float *f) {
/* return 1 on error, 0 on no error */
	double d = 0, min = (double) FLT_MIN, max = (double) FLT_MAX;
	int warning = 0;

	warning = read_double(s, &d);
	if (is_mv_double(&d)) {
		set_mv_float(f);
		return 0;
	}
	if (fabs(d) > max || (fabs(d) < min && d != 0.0)) {
		message("value outside valid range +/-[%g, %g]\n", min, max);
		ErrMsg(ER_RANGE, s);
	}
	*f = (float) d;
	return warning;
}

int read_double(const char *s, double *d) {
/* return 1 on error, 0 on no error */
	char *cp;

	if (s == NULL)
		ErrMsg(ER_NULL, "read_double()");
	if (s[0] == '\0') 
		ErrMsg(ER_IMPOSVAL, "read_double(): empty string");
	if (strcmp(s, gl_mv_string) == 0) {
		set_mv_double(d);
		return 0;
	}
	errno = 0;
	*d = strtod(s, &cp);
	if (errno == ERANGE) {
		message("value outside valid range +/-[%g, %g]\n", DBL_MIN, DBL_MAX);
		ErrMsg(ER_RANGE, s);
	}
	if (*cp == '\0')
		return 0;
	else {
#ifdef READ_WARNING
		pr_warning("read_double(): unconverted suffix: `%s'", cp);
#endif
		return 1;
	}
}

int read_int(const char *s, int *i) {
/* return 1 on error, 0 on no error */
	long int l;
	int warning = 0;

	warning = read_long(s, &l);
	if (warning) 
		return warning;
	if (l > INT_MAX || l < INT_MIN) {
		message("value outside valid range [%d, %d]\n", INT_MIN, INT_MAX);
		ErrMsg(ER_RANGE, s);
	}
	*i = (int) l;
	return warning;
}

int read_uint(const char *s, unsigned int *u) {
/* return 1 on error, 0 on no error */
	unsigned long ul;
	int warning = 0;

	warning = read_ulong(s, &ul);
	if (warning) 
		return warning;
	if (ul > UINT_MAX) {
		message("value outside valid range [0, %u]\n", UINT_MAX);
		ErrMsg(ER_RANGE, s);
	}
	*u = (unsigned int) ul;
	return warning;
}

int read_long(const char *s, long int *l) {
/* return 1 on error, 0 on no error */
	char *cp;

	if (s == NULL)
		ErrMsg(ER_NULL, "read_long()");
	if (s[0] == '\0') 
		ErrMsg(ER_IMPOSVAL, "read_long(): empty string");
	errno = 0;
	*l = strtol(s, &cp, 10);
	if (errno == ERANGE) {
		message("value outside valid range [%ld, %ld]\n", LONG_MIN, LONG_MAX);
		ErrMsg(ER_RANGE, s);
	}
	if (*cp == '\0')
		return 0;
	else  {
#ifdef READ_WARNING
		pr_warning("read_long(): unconverted suffix: `%s'", cp);
#endif
		return 1;
	}
}

int read_ulong(const char *s, unsigned long *u) {
/* return 1 on error, 0 on no error */
	char *cp;

	if (s == NULL)
		ErrMsg(ER_NULL, "read_long()");
	if (s[0] == '\0') 
		ErrMsg(ER_IMPOSVAL, "read_long(): empty string");
	errno = 0;
	*u = strtoul(s, &cp, 10);
	if (errno == ERANGE) {
		message("value outside valid range [0, %lu]\n", ULONG_MAX);
		ErrMsg(ER_RANGE, s);
	}
	if (*cp == '\0')
		return 0;
	else  {
#ifdef READ_WARNING
		pr_warning("read_long(): unconverted suffix: `%s'", cp);
#endif
		return 1;
	}
}
