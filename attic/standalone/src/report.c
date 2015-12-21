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
 * report.c: write (cross) validation report
 */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "defs.h"

#ifdef HAVE_LIBGIS
#include "gis.h"
#include "site.h"
#endif

#include "userio.h"
#include "data.h"
#include "utils.h" /* my_dtoa() */
#include "debug.h"
#include "glvars.h"
#include "stat.h"
#include "predict.h"
#include "report.h"

static void write_ascii_header(FILE *out_file, DATA *data, 
	int n_outfl);

void report_xvalid(double *xdata, double *xpred, double *xdiff, double *xstd,
		double *xzscore, int ndata, int var) {	
/*
 * DATE: Tue Oct  6 11:55:44 MET 1992
 * BY  : Edzer J. Pebesma
 * PURPOSE: report summary statistics of these five lists
 * SIDE EFFECTS: none
 */
	int i, nXdata = 0, nXpred = 0, nXdiff = 0, n_std = 0, nZscore = 0,
		compare(const double *a, const double *b);
	double min[5], max[5], p25[5], p75[5], p50[5], mean[5], std[5];
	double corr = 0.0;

	set_mv_double(&corr);
	calc_r(xdata, xpred, ndata, &corr);
	for (i = 0; i < 5; i ++) {
		set_mv_double(&(min[i])); 
		set_mv_double(&(p25[i])); 
		set_mv_double(&(p50[i])); 
		set_mv_double(&(p75[i])); 
		set_mv_double(&(max[i])); 
		set_mv_double(&(mean[i])); 
		set_mv_double(&(std[i])); 
	}
	/* select not missing values, put mv's at the end: */
	/* sorting arrays: */
	qsort(xdata, (size_t) ndata, sizeof(double), 
			(int CDECL (*)(const void *,const void *)) compare);
	while (!is_mv_double(&(xdata[nXdata])) && nXdata < ndata)
		nXdata++;
	qsort(xpred, (size_t) ndata, sizeof(double), 
			(int CDECL (*)(const void *,const void *)) compare);
	while (!is_mv_double(&(xpred[nXpred])) && nXpred < ndata)
		nXpred++;
	qsort(xdiff, (size_t) ndata, sizeof(double), 
			(int CDECL (*)(const void *,const void *)) compare);
	while (!is_mv_double(&(xdiff[nXdiff])) && nXdiff < ndata)
		nXdiff++;
	if (var) { /* do everything for xstd and xzscore */
		qsort(xstd, (size_t) ndata, sizeof(double), 
			(int CDECL (*)(const void *,const void *)) compare);
		while ((! is_mv_double(&(xstd[n_std]))) && (n_std < ndata))
			n_std++;
		qsort(xzscore, (size_t) ndata, sizeof(double), 
			(int CDECL (*)(const void *,const void *)) compare);
		while ((! is_mv_double(&(xzscore[nZscore]))) && (nZscore < ndata))
			nZscore++;
	}
	/* calculate statistics: */
	if (nXdata) {
		min[0]=xdata[0];
		max[0]=xdata[nXdata-1];
		mean[0] = sample_mean(xdata, nXdata); 
		if (nXdata > 1) {
			p25[0]=est_quant(xdata, 0.25, nXdata);
			p50[0]=est_quant(xdata, 0.5, nXdata);
			p75[0]=est_quant(xdata, 0.75, nXdata);
			std[0] = sample_std(xdata, mean[0], nXdata);
		}
	}

	if (nXpred) {
		min[1]=xpred[0];
		max[1]=xpred[nXpred-1];
		mean[1] = sample_mean(xpred, nXpred);
		if (nXpred > 1) {
			p25[1]=est_quant(xpred, 0.25, nXpred);
			p50[1]=est_quant(xpred, 0.5, nXpred);
			p75[1]=est_quant(xpred, 0.75, nXpred);
			std[1] = sample_std(xpred, mean[1], nXpred);
		}
	}

	if (nXdiff) {
		min[2]=xdiff[0];
		max[2]=xdiff[nXdiff-1];
		mean[2] = sample_mean(xdiff, nXdiff);
		if (nXdiff > 1) {
			p25[2]=est_quant(xdiff, 0.25, nXdiff);
			p50[2]=est_quant(xdiff, 0.5, nXdiff);
			p75[2]=est_quant(xdiff, 0.75, nXdiff);
			std[2] = sample_std(xdiff, mean[2], nXdiff);
		}
	}

	if (var) {
		if (n_std) {
			min[3]=xstd[0];
			max[3]=xstd[n_std-1];
			mean[3] = sample_mean(xstd, n_std);
			if (n_std > 1) {
				p25[3]=est_quant(xstd, 0.25, n_std);
				p50[3]=est_quant(xstd, 0.5, n_std);
				p75[3]=est_quant(xstd, 0.75, n_std);
				std[3] = sample_std(xstd, mean[3], n_std);
			}
		}
		if (nZscore) {
			min[4]=xzscore[0];
			max[4]=xzscore[nZscore-1];
			mean[4] = sample_mean(xzscore, nZscore);
			if (nZscore > 1) {
				p25[4]=est_quant(xzscore, 0.25, nZscore);
				p50[4]=est_quant(xzscore, 0.5, nZscore);
				p75[4]=est_quant(xzscore, 0.75, nZscore);
				std[4] = sample_std(xzscore, mean[4], nZscore);
			}
		}
	}

	/* output: */
	printlog("corr(Obs, Pred): %s  [%s]\n\n",
		my_dtoa("%6.4g", &corr),
		method_string(get_method()));
	printlog("              observed   predicted   pred.-obs.   pred.std.     zscore\n");
	printlog("======================================================================\n");
	printlog("%-10s%12s", "minimum", my_dtoa("%6.4g", &(min[0])));
	printlog("%12s", my_dtoa("%6.4g", &(min[1]))); 
	printlog("%12s", my_dtoa("%6.4g", &(min[2]))); 
	printlog("%12s", my_dtoa("%6.4g", &(min[3]))); 
	printlog("%12s\n", my_dtoa("%6.4g", &(min[4]))); 
	printlog("%-10s%12s", "1st q.", my_dtoa("%6.4g", &(p25[0])));
	printlog("%12s", my_dtoa("%6.4g", &(p25[1]))); 
	printlog("%12s", my_dtoa("%6.4g", &(p25[2]))); 
	printlog("%12s", my_dtoa("%6.4g", &(p25[3]))); 
	printlog("%12s\n", my_dtoa("%6.4g", &(p25[4]))); 
	printlog("%-10s%12s", "median", my_dtoa("%6.4g", &(p50[0])));
	printlog("%12s", my_dtoa("%6.4g", &(p50[1]))); 
	printlog("%12s", my_dtoa("%6.4g", &(p50[2]))); 
	printlog("%12s", my_dtoa("%6.4g", &(p50[3]))); 
	printlog("%12s\n", my_dtoa("%6.4g", &(p50[4]))); 
	printlog("%-10s%12s", "3rd q.", my_dtoa("%6.4g", &(p75[0])));
	printlog("%12s", my_dtoa("%6.4g", &(p75[1]))); 
	printlog("%12s", my_dtoa("%6.4g", &(p75[2]))); 
	printlog("%12s", my_dtoa("%6.4g", &(p75[3]))); 
	printlog("%12s\n", my_dtoa("%6.4g", &(p75[4]))); 
	printlog("%-10s%12s", "maximum", my_dtoa("%6.4g", &(max[0])));
	printlog("%12s", my_dtoa("%6.4g", &(max[1]))); 
	printlog("%12s", my_dtoa("%6.4g", &(max[2]))); 
	printlog("%12s", my_dtoa("%6.4g", &(max[3]))); 
	printlog("%12s\n\n", my_dtoa("%6.4g", &(max[4]))); 
	printlog("%-10s%12d%12d%12d%12d%12d\n", "n",
		nXdata, nXpred, nXdiff, n_std, nZscore);
	printlog("%-10s%12s", "mean", my_dtoa("%6.4g", &(mean[0])));
	printlog("%12s", my_dtoa("%6.4g", &(mean[1]))); 
	printlog("%12s", my_dtoa("%6.4g", &(mean[2]))); 
	printlog("%12s", my_dtoa("%6.4g", &(mean[3]))); 
	printlog("%12s\n", my_dtoa("%6.4g", &(mean[4]))); 
	printlog("%-10s%12s", "std.dev.", my_dtoa("%6.4g", &(std[0])));
	printlog("%12s", my_dtoa("%6.4g", &(std[1]))); 
	printlog("%12s", my_dtoa("%6.4g", &(std[2]))); 
	printlog("%12s", my_dtoa("%6.4g", &(std[3]))); 
	printlog("%12s\n", my_dtoa("%6.4g", &(std[4]))); 
	return;
}

int CDECL compare(const double *a, const double *b) 
/* ansi conformant qsort cmp, puts mv's at the end */
{
	if (is_mv_double(a)) /* a is bigger */
		return 1;
	if (is_mv_double(b)) /* b is bigger */
		return -1;

	if (*a < *b) 
		return -1; 
	if (*a > *b) 
		return 1; 
	return 0;
}

#ifndef USING_R
static void write_ascii_header(FILE *out_file, DATA *data, int n_outfl) {
	char *lf = "\n";
	int i = 0;

	if (data->mode & X_BIT_SET) i++; 		/* xcoord */
	if (data->mode & Y_BIT_SET) i++; 		/* ycoord */
	if (data->mode & Z_BIT_SET) i++; 		/* zcoord */
	if (data->mode & S_BIT_SET) i++;
	if (data->colnvalue > 0) i++;	/* obs */
    if (data->point_ids) i++;
    
	if (get_mode() == STRATIFY)
		i+= 2;
	else
		i+= n_outfl; 					/* est[0], .. , est[n_outfl-1]  */
	if (data->type.type == DATA_EAS) {
		fprintf(out_file, "prediction results from file: %s variable: %s\n",
			data->fname, data->variable); 
		fprintf(out_file, "%d\n", i);
	} else {
		fprintf(out_file, "#");
		lf = " ";
	}
	if (data->mode & X_BIT_SET) fprintf(out_file, "%s%s", NULS(data->x_coord), lf);
	if (data->mode & Y_BIT_SET) fprintf(out_file, "%s%s", NULS(data->y_coord), lf);
	if (data->mode & Z_BIT_SET) fprintf(out_file, "%s%s", NULS(data->z_coord), lf);
	if (data->mode & S_BIT_SET) fprintf(out_file, "%s%s", NULS(data->s_coord), lf);
	if (data->colnvalue > 0) fprintf(out_file, "%s%s", NULS(data->variable), lf);
	if (get_mode() == STRATIFY)
		fprintf(out_file, "predictions%svariances%s", lf, lf);
	else
		for (i = 0; i < n_outfl; i++)
			fprintf(out_file, "%s%s", what_is_outfile(i), lf);

    if (data->point_ids)
        fprintf(out_file, "%s%s",data->id_name,lf);

    if (data->type.type != DATA_EAS) fprintf(out_file, "\n");
    
}

static void output_line(FILE *out_file, DATA *data, DPOINT *where, 
		double *est, int n_outfl) {
	int i;

	assert(out_file != NULL);
	assert(out_file != stdout);

	if (data->mode & X_BIT_SET) {
		fprintf(out_file, " ");
		fprintf(out_file, gl_format, where->x);
	}
	if (data->mode & Y_BIT_SET) {
		fprintf(out_file, " ");
		fprintf(out_file, gl_format, where->y);
	}
	if (data->mode & Z_BIT_SET) {
		fprintf(out_file, " ");
		fprintf(out_file, gl_format, where->z);
	}
	if (data->mode & S_BIT_SET)
		fprintf(out_file, " %d", where->u.stratum + strata_min);
	if (data->colnvalue > 0) {
		fprintf(out_file, " ");
		fprintf(out_file, gl_format, where->attr);
	}
	for (i = 0; i < n_outfl; i++)
		fprintf(out_file, " %s", my_dtoa(gl_format, &(est[i])));

    if (data->point_ids)
        fprintf(out_file, " %s",data->point_ids[GET_INDEX(where)]);
    
	fprintf(out_file, "\n");
	return;
}


void write_points(const char *fname, DATA *d, DPOINT *where, double *est,
	int n_outfl) {

	static FILE *f = NULL;
#ifdef HAVE_LIBGIS
	static Site *site = NULL;
	static int dim = 2;
	int i;
#endif 

	if (! grass()) {
		if (where == NULL) {
			if (fname != NULL) {
				f = efopen(fname, "w");
				write_ascii_header(f, d, n_outfl);
			} else
				efclose(f);
		} else {
			if (f == NULL)
				ErrMsg(ER_NULL, "write_points(): f");
			output_line(f, d, where, est, n_outfl);
		}
	} else {
#ifdef HAVE_LIBGIS
		if (where == NULL) {
			if (fname != NULL) { /* initialize: */
				DUMP("opening grass sites list\n");
				if (d->mode & Z_BIT_SET)
					dim++;
				if ((f = G_sites_open_new((char *) fname)) == NULL)
					G_fatal_error("%s: cannot open sites file %f for writing\n",
						G_program_name());
				site = G_site_new_struct(CELL_TYPE, dim, 0, n_outfl);
			} else { /* close: */
				DUMP("closing grass sites list\n");
				fclose(f);
				dim = 2;
				G_site_free_struct(site);
				site = NULL;
			}
		} else {
			assert(site != NULL);
			assert(d != NULL);
			/* fill site: */
			site->east = where->x;
			site->north = where->y;
			if (d->mode & Z_BIT_SET)
				site->dim[0] = where->z;
			if (d->mode & S_BIT_SET)
				site->ccat = where->u.stratum + strata_min;
			else
				site->ccat = GET_INDEX(where) + 1;
			for (i = 0; i < n_outfl; i++) {
				if (is_mv_double(&(est[i]))) {
					site->dbl_att[i] = -9999.0;
					if (DEBUG_DUMP)
						printlog(" [%d]:mv ", i);
				} else {
					site->dbl_att[i] = est[i];
					if (DEBUG_DUMP)
						printlog(" value[%d]: %g ", i, site->dbl_att[i]);
				}
			}
			if (DEBUG_DUMP)
				printlog("\n");
			G_site_put(f, site);
		}
#else 
		ErrMsg(ER_IMPOSVAL, "gstat/grass error: libgis() not linked");
#endif 
	}
}
#endif /* USING_R */
