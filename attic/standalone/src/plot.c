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
 * plot.c: kriging weights/variogram plot drivers to gnuplot and jgraph
 */

#include <stdio.h>
#include <stdlib.h> /* getenv() */
#include <string.h>
#include <math.h>

#include "defs.h"
#include "userio.h"
#include "debug.h"
#include "data.h"
#include "utils.h"
#include "vario.h"
#include "glvars.h"
#include "sem.h"
#include "version.h"
#include "plot.h"

/* jgraph defaults: */
#define SIZE 0.013 /* rel distance (in) of number from cross */
#define N_PTS 200 /* number of semivariogram model points for jgraph */
#define PLOT_WINDOW_NR(v) \
		(v->id1==v->id2?v->id1:(get_n_vars()+LTI2(v->id1,v->id2)))

typedef struct {
	PLOT_TYPE p;
	int pt, lt;
	char *term_name, *term_opt, *gamma_str, *distance_str;
} GNUPLOT_TERM;

GNUPLOT_TERM gnuplot_terms[] = {
 { GNUPLOT, 2, 3, NULL, NULL, NULL, NULL }, /* use default option */
 { PSLATEX, 1, 1, "pslatex", "", "$\\gamma(h)$", "$h$" },
 { CGM, 1, 1, "cgm", "", NULL },
 { GIF, 2, 2, "gif", 
 	"transparent size 640, 480 xffffff x000000 x404040 xff0000 x0000ff x00ff00",
	NULL, NULL },
 { EPS, 1, 1, "postscript", "eps solid 17", NULL, NULL },
 { PNG, 2, 3, 
 		"png", "small transparent xffffff x000000 x404040 xff0000",
		NULL, NULL },
 { EEPIC, 2, 3, "eepic", "", "{$\\gamma(h)$}", "{$h$}" },
 { UNKNOWN, 0, 0, NULL, NULL, NULL, NULL }
};

static int set_key(FILE *f, const VARIOGRAM *v, double min, double max);
static void get_minmax_y(const VARIOGRAM *v, double max_x, double *min_y, double *max_y);

#ifndef USING_R
int fprint_gnuplot_variogram(FILE *stream, const VARIOGRAM *v,
		char *fname, PLOT_TYPE p, int window_nr) {
	double min_y = 0.0, max_y = 0.0, max_x = 0.0, range;
	int i, key_outside = 1;
	GNUPLOT_TERM t = gnuplot_terms[0];
	const char *est = NULL; 
	int inline_data = 0;

	assert(v->ev != NULL);
	assert(stream != NULL);
	/* plot sample vgm and/or model: */
	assert(v->ev->n_est > 0 || v->n_models > 0); 

	for (i = 0; gnuplot_terms[i].p != UNKNOWN; i++) /* find the entry */
		if (gnuplot_terms[i].p == p) {
			t = gnuplot_terms[i];
			break;
		}
	assert(t.p != UNKNOWN);

	if (gl_gpterm != NULL) {
		t.term_name = gl_gpterm;
		t.term_opt = "";
	}

	fprintf(stream, "#\n# gnuplot file, created by gstat %s\n#\n\n", VERSION);
	fprintf(stream, "set tics out\n");
	if (fname) {
		fprintf(stream, "## for %s term:\n", t.term_name);
		if (t.term_name != NULL)
			fprintf(stream, "set term %s %s\n", t.term_name, t.term_opt);

		if (p == PSLATEX || p == EPS || p == EEPIC) {
			fprintf(stream, "set border 3\n");
			fprintf(stream, "set xtics nomirror\nset ytics nomirror\n");
		}

		/* output file: */
		if (*fname != '\0')
			fprintf(stream, "%sset output '%s'\n\n", 
				p == GNUPLOT ? "# " : "", fname);
	} else { /* plot to screen terminal */
#ifdef WIN32
	    fprintf(stream, "set term windows\n");
#else
		if (gl_gnuplot35 == NULL && getenv("DISPLAY") != NULL
				&& getenv("GNUTERM") == NULL) 
			/* will default to x11; set nr for 3.7: */
			fprintf(stream, "set term x11 %d\n", 
					window_nr >= 0 ? window_nr : PLOT_WINDOW_NR(v));
#endif
	}

	/* print basic models: */
	for (i = 1; v_models[i].name != NULL; i++)
		fprintf(stream, "%s\n", 
			is_variogram(v) ? v_models[i].v_gnuplot : v_models[i].c_gnuplot);

	fprintf(stream, "\n# erase previous label settings:\nset nolabel\n");
	if (p != EEPIC)
		fprintf(stream, "\n# set dottet lines at y=0:\nset zeroaxis\n");
	fprintf(stream, "\n# set x and y labels:\n");
	fprintf(stream, "set xlabel '%s'\n", 
		t.distance_str == NULL ? "distance" : t.distance_str);
	fprintf(stream, "set ylabel '");
	switch (v->ev->evt) {
		case SEMIVARIOGRAM: 
			fprintf(stream, "%s", 
				t.gamma_str == NULL ? "semivariance" : t.gamma_str);
			break;
		case CROSSVARIOGRAM: 
			if (v->ev->pseudo)
				fprintf(stream, "%s", "pseudo cross semivariance");
			else
				fprintf(stream, "%s", "cross semivariance");
			break;
		case COVARIOGRAM: 
			fprintf(stream, "%s", "C(h)");
			break;
		case CROSSCOVARIOGRAM: 
			fprintf(stream, "%s", "cross covariance");
			break;
		case NOTSPECIFIED: 
			ErrMsg(ER_IMPOSVAL, "evt unspecified");
			break;
	} 
	if (v->ev->cloud)
		fprintf(stream, " cloud");
	if (v->ev->is_directional) {
		fprintf(stream, " (dir. ");
		if (gl_tol_hor < 90) {
			fprintf(stream, "<x,y> %g +/- %g", gl_alpha, gl_tol_hor);
			if (gl_tol_ver < 90)
				fprintf(stream, "; ");
		}
		if (gl_tol_ver < 90)
			fprintf(stream, "<z> %g +/- %g", gl_beta, gl_tol_ver);
		fprintf(stream, ")");
	}
	fprintf(stream, "'\n\n");
	if (v->ev->plot_numbers)
		fprintf(stream, "\n# set number of point pairs as labels:\n");
	for (i = 0; i < (v->ev->zero == ZERO_AVOID ? 
			v->ev->n_est - 1 : v->ev->n_est); i++) {
		if (v->ev->nh[i] > 0) {
			if (v->ev->plot_numbers && v->ev->n_est <= gl_dots) {
				if (v->ev->cloud)
					fprintf(stream, "set label '%ld,%ld' at %e,%e left\n", 
						HIGH_NH(v->ev->nh[i]) + 1, LOW_NH(v->ev->nh[i]) + 1,
						v->ev->dist[i] + SIZE * v->ev->cutoff, v->ev->gamma[i]);
				else
					fprintf(stream, "set label '%ld' at %e,%e left\n", 
						v->ev->nh[i], v->ev->dist[i] + SIZE * v->ev->cutoff,
						v->ev->gamma[i]);
			}
		}
	}
	if (v->ev->n_est == 0) 
		max_x = 2.5 * effective_range(v);
	else
		max_x = 1.05 * v->ev->cutoff;
	get_minmax_y(v, max_x, &min_y, &max_y);
	if (gl_gnuplot35 == NULL)
		key_outside = set_key(stream, v, min_y, max_y);
	fprintf(stream, "\n# set axis limits:\n");
	fprintf(stream, "xmin = 0; xmax = %e\n", max_x);
	fprintf(stream, "set xrange [xmin:xmax]\n"); 
	if (max_y > min_y) {
		range = max_y - min_y;
		fprintf(stream, "ymin = %e; ymax = %e\n", 
			min_y < 0.0 ? min_y - 0.05 * range : 0.0,
			max_y + 0.05 * range + 0.15 * key_outside * range);
		fprintf(stream, "set yrange [ymin:ymax]\n"); 
	}
	/* plot file: */

	fprintf(stream, "\n# do the data & function plotting:\n");

	fprintf(stream, "plot ");

	if (v->ev->n_est > 0) {	
		if (v->fname == NULL && o_filename == NULL) {
			inline_data = 1;
			est = "-";
		} else
			est = v->fname ? v->fname : o_filename;
	
		fprintf(stream, "'%s'", est);
		/* using: */
		if (v->ev->cloud)
			fprintf(stream, "%s", " using 3:4 ");
		else
			fprintf(stream, "%s", " using 4:5 ");
		/* title: */
		if (v->id1 == v->id2)
			fprintf(stream, "title '%s' ", name_identifier(v->id1));
		else
			fprintf(stream, "title '%s x %s' ", name_identifier(v->id1), 
				name_identifier(v->id2));
		if (v->ev->n_est > gl_dots)
			fprintf(stream, "with dots 1");
		else
			fprintf(stream, "with points pt %d", t.pt);
	}
	if (v->n_models > 0) {
		fprintf(stream, "%s\\\n ", v->ev->n_est > 0 ? "," : "");
		fprint_gnuplot_model(stream, v, 0);
		fprintf(stream, "\\\n title '%s' with lines lt %d", v->descr, t.lt);
		/*
		if (p == PSLATEX)
			fprintf(stream, " with lines 1 1");
		*/
	}
	fprintf(stream, "\n"); /* do it ! */
	if (inline_data) {
		fprint_sample_vgm(stream, v->ev);
		fprintf(stream, "e\n");
	}

#ifdef WIN32
	/* prevent wgnuplot from disappearing; kms */
   	fprintf(stream, "pause -1 \"Hit Enter to continue\"\n");
#endif

	fflush(stream);
	return 0;
}

void fprint_gnuplot_model(FILE *f, const VARIOGRAM *vgm, int fit) {
	int i;
	char a[50], c[50];

	for (i = 0; i < vgm->n_models; i++) {
		if (fit && vgm->part[i].fit_sill) /* print parameter to string: */
			sprintf(c, "c%d", i);
		else
			sprintf(c, "%.2e", vgm->part[i].sill);

		if (fit && vgm->part[i].fit_range) /* print parameter to string: */
			sprintf(a, "a%d", i);
		else
			sprintf(a, "%.2e", vgm->part[i].range[0] *
				relative_norm(vgm->part[i].tm_range, vgm->ev->direction.x,
				vgm->ev->direction.y, vgm->ev->direction.z));

		fprintf(f, "%s * %s(%s,x)", c,
			v_models[vgm->part[i].model].name, a);

		if (i < vgm->n_models - 1)
			fprintf(f, " + ");
	}
}

int fprint_jgraph_variogram(FILE *f, const VARIOGRAM *v) {
/*
 *  print variogram to a jgraph input file 
 *  (jgraph is a program that converts such a file to PS/EPS)
 *  URL: file://ftp.princeton.edu/pub/jgraph.Z
 *
 *  don't play tricks: f might be a pipe.
 */
	int i, sill_reached = 0;
	double max_y = 0.0, min_y = 0.0, max_x = 0.0, x, y, z, s;

	if (v->ev->n_est <= 0)
		return 1;
	fprintf(f, "(*\n  jgraph file, created by gstat %s\n*)\n",
		VERSION);
	fprintf(f, "newgraph\n");
	fprintf(f, "(* draw estimates as crosses: *)\n");
	fprintf(f, "newcurve marktype cross label : ");
	if (v->id1 == v->id2)
		fprintf(f, "%s\n pts\n", name_identifier(v->id1));
	else
		fprintf(f, "%s x %s\n pts\n", name_identifier(v->id1),
			name_identifier(v->id2));
	for (i = 0; i < ((v->ev->zero == ZERO_AVOID) ?
			v->ev->n_est - 1 : v->ev->n_est); i++)
		if (v->ev->nh[i] > 0)
			fprintf(f, "%g %g\n", v->ev->dist[i], v->ev->gamma[i]);
	if (v->ev->plot_numbers) {
		fprintf(f, "(* draw numbers *)\n");
		for (i = 0; i < v->ev->n_est; i++) {
			if (v->ev->nh[i] > 0)
				fprintf(f, "newstring x %g y %g hjl vjc fontsize 8 : %ld\n",
					v->ev->dist[i] + SIZE * v->ev->cutoff,
					v->ev->gamma[i], v->ev->nh[i]);
		} /* for */
		max_x = (1.0 + 2.0 * SIZE) * v->ev->cutoff;
	} else
		max_x = v->ev->cutoff;
	get_minmax_y(v, max_x, &min_y, &max_y);
/*
 * LEGEND, TITLES: 
 */
	fprintf(f, "(* draw legend, title and axis *)\nlegend top\n");
	fprintf(f, "(* title : (fill in title) *)\n");
/*
 * XAXIS:
 */
	fprintf(f, "xaxis min 0 max %g size 4.0 label : Distance\n", max_x);
	if (v->ev->cutoff > 99999.0 && v->ev->cutoff < 100001.0)
		fprintf(f, "hash 25000 mhash 4\n");
	if (v->ev->plot_numbers)
		fprintf(f, "(* Note: remove `max ...' if not correct *)\n");
/*
 * YAXIS: 
 */
	min_y *= 1.05;
	max_y *= 1.05;
	if (min_y == max_y) {
		min_y = -1.0;
		max_y =  1.0;
	}
	fprintf(f, "yaxis min %g ", min_y);
	if (max_y > 0)
		fprintf(f, "max %g ", max_y);
	else
		fprintf(f, "(* max 0 *) hash_format g ");
	fprintf(f, "size 2.5 label : ");
	switch (v->ev->evt) {
		case SEMIVARIOGRAM: 
			fprintf(f, "Semivariance\n");
			break;
		case CROSSVARIOGRAM: 
			if (v->ev->pseudo)
				fprintf(f, "Pseudo Cross Semivariance\n");
			else
				fprintf(f, "Cross Semivariance\n");
			break;
		case COVARIOGRAM: 
			fprintf(f, "Covariance\n");
			break;
		case CROSSCOVARIOGRAM: 
			fprintf(f, "Cross Covariance\n");
			break;
		case NOTSPECIFIED: 
			ErrMsg(ER_IMPOSVAL, "evt unspecified");
			break;
	} 
/* 
 * process model: 
 */
	if (v->n_models > 0) {
		fprintf(f, "(* draw variogram model as line: *)\n");
		fprintf(f, "newline label : %s\npts\n", v->descr);
		for (i = 0, x = 0; i < N_PTS; i++) {
			if (sill_reached)
				i = N_PTS - 1; /* jump to last one */
			x = (i * v->ev->direction.x) / (N_PTS-1.0) * max_x;
			y = (i * v->ev->direction.y) / (N_PTS-1.0) * max_x;
			z = (i * v->ev->direction.z) / (N_PTS-1.0) * max_x;
			if (x < 1.0e-30 && y < 1.0e-30 && z < 1.0e-30)
				x = 1.0e-30; /* avoid line from 0 to nugget */
			s = (is_variogram(v) ?
				get_semivariance(v, x, y, z) :
				get_covariance(v, x, y, z));
			fprintf(f, "\t%g %g\n", transform_norm(NULL, x, y, z), s);
			if (x > v->max_range)
				sill_reached = 1;
		} /* for i */
	} 
	return 0;
}
#endif

static void get_minmax_y(const VARIOGRAM *v, double max_x, double *min_y, double *max_y) {
	int i;
	double x, y, z, s;

	for (i = 0; i < v->ev->n_est; i++) { 
		if (v->ev->nh[i] > 0) {
			*max_y = MAX(*max_y, v->ev->gamma[i]); /* max_y >= 0 */
			*min_y = MIN(*min_y, v->ev->gamma[i]); /* min_y <= 0 */
		}
	}
	if (v->n_models > 0) { /* get max_y, min_y: */
		for (i = 0, x = 0; i < N_PTS; i++) {
			x = (i * v->ev->direction.x) / (N_PTS-1.0) * max_x;
			y = (i * v->ev->direction.y) / (N_PTS-1.0) * max_x;
			z = (i * v->ev->direction.z) / (N_PTS-1.0) * max_x;
			if (x < 1.0e-30 && y < 1.0e-30 && z < 1.0e-30)
				x = 1.0e-30; /* avoid line from 0 to nugget */
			s = (is_variogram(v) ?
				get_semivariance(v, x, y, z) :
				get_covariance(v, x, y, z));
			*max_y = MAX(*max_y, s);
			*min_y = MIN(*min_y, s);
		}
	}
}

static int set_key(FILE *f, const VARIOGRAM *v, double min, double max) {
/*
 * returns 1 if key should be set outside plot region
 * returns 0 if it fits inside plot region (and sets key top or bottom)
 */
	int i;
	enum { POS_INIT, POS_TOP, POS_BOTTOM, POS_OUTSIDE } pos = POS_INIT;

	min = min + 0.25 * (max - min);
	max = max - 0.25 * (max - min);
	for (i = 0; i < ((v->ev->zero == ZERO_AVOID) ? 
			v->ev->n_est - 1 : v->ev->n_est); i++) {
		if (v->ev->nh[i] > 0 && v->ev->dist[i] > 0.5 * v->ev->cutoff) {
			if (pos == POS_INIT) {
				pos = POS_OUTSIDE;
				if (v->ev->gamma[i] < max)
					pos = POS_TOP;
				if (v->ev->gamma[i] > min)
					pos = POS_BOTTOM;
			} else {
				if ((pos == POS_TOP && v->ev->gamma[i] > max) ||
					(pos == POS_BOTTOM && v->ev->gamma[i] < min))
					pos = POS_OUTSIDE;
			}
		}
	}
	switch (pos) {
		case POS_INIT:
		case POS_OUTSIDE:
			fprintf(f, "\n# set key outside plotting region:\nset key\n");
			return 1;
		case POS_TOP:
			fprintf(f, "\n# set key in plotting region:\nset key top\n");
			break;
		case POS_BOTTOM:
			fprintf(f, "\n# set key in plotting region:\nset key bottom\n");
			break;
	}
 	return 0;
}
