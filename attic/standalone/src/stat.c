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
 * stat.c: simple basic statistic functions
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "defs.h"

#ifdef HAVE_LIBGSL
# include "gsl/gsl_statistics.h"
#endif

#include "userio.h"
#include "data.h"
#include "utils.h"
#include "debug.h"
#include "lex.h"
#include "read.h"
#include "glvars.h"
#include "stat.h"

double sample_mean(double *list, int n) {
	int i;
	double mn = 0.0;

	if (list == NULL)
		ErrMsg(ER_NULL, "sample_mean()");
	if (n == 0) 
		ErrMsg(ER_RANGE, "sample_mean(): no values");
	for (i = 0; i < n; i++)
		mn += list[i];
	return mn/(1.0 * n);
}

double sample_var(double *list, double mean, int n) {
	int i;
	double var = 0.0;

	if (list == NULL)
		ErrMsg(ER_NULL, "sample_var()");
	if (n <= 1 || list == NULL) 
		ErrMsg(ER_RANGE, "sample_var(): <= 1 values");
	for (i = 0; i < n; i++)
		var += SQR(list[i] - mean);
	return (var / (n - 1.0));
}

double sample_std(double *list, double mean, int n) {
	return sqrt(sample_var(list, mean, n));
}

double est_quant(double *list, double p, int n) {
/*
* function returns the value of the p-th quantile
* of the ordered list *list, row exists from list[0]..list[length-1],
* p is in [0..1];
* a missing value is generated when the quantile lies not within
* the valid data range, and is undetermined therefore
*/
#ifndef HAVE_LIBGSL
	double order, where;
	int below, above;

	if (n < 2)
		ErrMsg(ER_RANGE, "est_quant(): < 2 obs.");
	if (p < 0.0 || p > 1.0)
		ErrMsg(ER_RANGE, "can't calculate quantile outside [0,1]");
	order = p * (n - 1);
	/* order = n * (p * n)/(n + 1); */
	below = (int) floor(order); /* the index below order */
	if (below < 0) 
		return list[0];
	above = below + 1;			/* the index above order */
	if (above >= n) 
		return list[n - 1];
	where = order - below;
	return (1 - where) * list[below] + where * list[above]; 
#else
	return gsl_stats_quantile_from_sorted_data(list, 1, n, p);
#endif
}

void calc_r(double *a, double *b, int n, double *corr) {
	double mean[2], sp = 0.0, ss1 = 0.0, ss2 = 0.0;
	int i, j;

	set_mv_double(corr);
	/* calc r: */
	mean[0] = 0.0; 
	mean[1] = 0.0;
	for (i = 0, j = 0; i < n; i++) 
		if (!is_mv_double(&(b[i])) && !is_mv_double(&(a[i]))) {
			mean[0] += a[i];
			mean[1] += b[i];
			j++;
		}
	if (j == 0)
		return;
	mean[0] /= j;
	mean[1] /= j;
	for (i = 0; i < n; i++) 
		if (! is_mv_double(&(b[i])) && !is_mv_double(&(a[i]))) {
			sp += (a[i] - mean[0]) * (b[i] - mean[1]);
			ss1 += SQR(b[i]-mean[1]);	
			ss2 += SQR(a[i]-mean[0]);	
		}
	if (ss1 > 0.0 && ss2 > 0.0)
		*corr = sp/(sqrt(ss1) * sqrt(ss2));
	return;
}

#ifndef USING_R
int stats(char *name, int silent, double q) {
	static D_VECTOR *dv = NULL;
	double mean = 0.0;

	if (dv == NULL) {
		dv = (D_VECTOR *) emalloc(sizeof(D_VECTOR));
		dv->max_size = dv->size = 0;
		dv->val = NULL;
	}
	dv->size = 0;

	read_vector(dv, name);

	assert(dv->size > 0);

	qsort(dv->val, (size_t) dv->size, sizeof(double),
		(int CDECL (*)(const void *, const void *)) d_cmp);

	mean = sample_mean(dv->val, dv->size);

	if (dv->size <= 1) {
		pr_warning("calc_stats(): n <= 1");
		return 0;
	}

	if (! silent) {
		if (name)
			printf("          ");
		if (q == 0.25)
			printf("     min    1st Q   median    3rd Q      max     mean      std      n\n");
		else
			printf("     min   Q(%.2f)  median   Q(%.2f)     max     mean      std      n\n",
				q, 1.0-q);
	}
	/* printf("%8.3g %8.3g %8.3g %8.3g %8.3g %8.3g %8.3g %6d\n", */
	if (name)
		printf("%-10s", name);
	printf("%8g %8g %8g %8g %8g %8g %8g %6d\n",
		dv->val[0], 
		est_quant(dv->val, q, dv->size), 
		est_quant(dv->val, .5, dv->size), 
		est_quant(dv->val, 1.0-q, dv->size), 
		dv->val[dv->size-1], 
		mean, 
		sample_std(dv->val, mean, dv->size), 
		dv->size);
	return 0;
}
#endif

int CDECL d_cmp(const double *a, const double *b) {
	if (*a < *b)
		return -1;
	if (*a > *b)
		return 1;
	return 0;
}
