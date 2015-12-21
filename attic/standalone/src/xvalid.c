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
 * xvalid.c: cross validation
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "defs.h"
#include "userio.h"
#include "data.h"
#include "utils.h"
#include "debug.h"
#include "glvars.h"
#include "getest.h"
#include "select.h"
#include "stat.h"
#include "report.h"
#include "xvalid.h"
 
#define IS_VAR(m) (m == SKR || m == OKR || m == UKR || m == LSLM)

static int cross = 0;

static void remove_where_from_selection(DATA *data, DPOINT *pt, int at0);

void cross_valid(DATA **data) {
/*
 * DATE: 		10 sept 1992
 * BY: 			Edzer J. Pebesma.
 * RETURNS:		none
 * PURPOSE:		perform cross validation on data[0]->list[i], 
 * given data[0]->list[j] (j = 0..i-1,i+1..data[0]->n) and data[j] (j > 0) 
 * PARAMETERS: data is a data array, data[0]->n_list > 0; 
 * SIDE EFFECT:	some memory allocation and free_ing
 */
	int i, j, k, var, sel_max_ori, n_Xtot;
	double *xdata, *xpred, *xdiff, *xstd, *xzscore;
	static double *est;
	double *X = NULL;

	DPOINT where, *block;

	cross = 1;
/* get global variables: */
	if (o_filename != NULL)
		write_points(o_filename, data[0], NULL, NULL, get_n_outfile());
	for (i = n_Xtot = 0; i < get_n_vars(); i++)
		n_Xtot += data[i]->n_X;
	X = (double *) emalloc(n_Xtot * sizeof(double));
	for (i = 0; i < n_Xtot; i++)
		X[i] = 1.0; /* will lead to correct results on ordinary cokriging */
	est = (double *) emalloc (get_n_outfile() * sizeof(double));
	sel_max_ori = data[0]->sel_max;
	data[0]->sel_max = MIN(data[0]->n_list, sel_max_ori) + 1; /* prevent "global" */
					 
	block = get_block_p();
	if (block->x > 0.0 || block->y > 0.0 || block->z > 0.0 || 
			get_data_area() != NULL)
		pr_warning("xvalid.c: cross validating with non-zero block/area size");

/* initialise: */
	var = IS_VAR(get_method());
/* var = TRUE if variance is returned, otherwise FALSE */
	xdata = (double *) emalloc(data[0]->n_list * sizeof(double));
	xpred = (double *) emalloc(data[0]->n_list * sizeof(double));
	xdiff = (double *) emalloc(data[0]->n_list * sizeof(double));
	xstd = (double *) emalloc(data[0]->n_list * sizeof(double));
	xzscore = (double *) emalloc(data[0]->n_list * sizeof(double));
	for (i = 0; i < data[0]->n_list; i++) { /* initialise all arrays: */
		set_mv_double(&(xdata[i]));
		set_mv_double(&(xpred[i]));
		set_mv_double(&(xdiff[i]));
		set_mv_double(&(xstd[i]));
		set_mv_double(&(xzscore[i]));
	}
	for (i = 0; i < data[0]->n_list; i++) { 
		xdata[i] = data[0]->list[i]->attr;
		where = *(data[0]->list[i]);
		for (j = 0; j < data[0]->n_X; j++)
			X[j] = where.X[j];
		where.X = X;
		/* qtree_rebuild(data); */
		for (j = 0; j < get_n_outfile(); j++)
			set_mv_double(&est[j]);
		for (j = 0; j < get_n_vars(); j++) {
			select_at(data[j], &where);
			remove_where_from_selection(data[j], &where, j == 0);
		}
		get_est(data, get_method(), &where, est);
		if (o_filename != NULL) {
			for (j = 1; j < get_n_vars(); j++) {
				if (!(data[j]->n_X == 1 && data[j]->colX[0] == 0)) {
					set_mv_double(&(est[2 * j]));
					set_mv_double(&(est[2 * j + 1]));
					for (k = 0; k < get_n_vars(); k++)
						if (k != j)
							set_mv_double(&(est[2 * get_n_vars() + LTI2(j,k)]));
				}
			}
			write_points(o_filename, data[0], &where, est, get_n_outfile());
		}
		if (! is_mv_double(&est[0])) {
			xpred[i] = est[0]; 
			if (! is_mv_double(&xdata[i]))
				xdiff[i] = xpred[i] - xdata[i];
			if (var) {
				if (! is_mv_double(&est[1])) 
					if (est[1] >= 0.0) {
						xstd[i] = sqrt(est[1]); 
						if (xstd[i] != 0.0 && !is_mv_double(&xdata[i]))
							xzscore[i] = xdiff[i]/xstd[i];
					}
			} 
		} else {
			set_mv_double(&xpred[i]); 
			set_mv_double(&xdiff[i]); 
			set_mv_double(&xstd[i]);
			set_mv_double(&xzscore[i]);
		} /* if .... */
		print_progress(i, data[0]->n_list);
	} /* for i, all records */
	/* 
	 * once and for all: copy original data pointers back: 
	 */
	data[0]->sel_max = sel_max_ori;
	print_progress(data[0]->n_list, data[0]->n_list);
	report_xvalid(xdata, xpred, xdiff, xstd, xzscore, data[0]->n_list, var);
	if (o_filename != NULL)
		write_points(NULL, NULL, NULL, NULL, 0); /* closes file */
	efree(xdata); 
	efree(xpred); 
	efree(xdiff); 
	efree(xstd); 
	efree(xzscore);
	efree(X);
	cross = 0;
	return;
}

static void remove_where_from_selection(DATA *data, DPOINT *pt, int at0) {
	/* remove point pt from data->sel */
	int n_sel, i;

	n_sel = data->n_sel;
	if (n_sel == 0)
		return;

	for (i = 0; i < data->n_sel; i++) {
		if (data->pp_norm2(data->sel[i], pt) == 0) {
			data->sel[i] = data->sel[data->n_sel-1]; /* copy last to i */
			data->n_sel--; /* cut 1 from list */
		}
	}
	if (at0 && n_sel == data->n_sel)
		ErrMsg(ER_NULL, "point(s) not found in remove_where_from_selection");
}
