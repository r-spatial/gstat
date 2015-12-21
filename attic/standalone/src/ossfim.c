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
#include <stdio.h>
#include <stdlib.h> 
#include <string.h> /* memcpy */
#include <math.h>

#include "defs.h"

#ifdef HAVE_GETOPT_H
# include <getopt.h>
#endif
#ifdef HAVE_UNISTD_H
# include <unistd.h>
#endif
#ifndef HAVE_GETOPT
# include "getopt.h"
#endif

#include "userio.h"
#include "debug.h"
#include "data.h"
#include "vario.h"
#include "utils.h"
#include "glvars.h"
#include "select.h"
#include "nsearch.h"
#include "parse.h"
#include "read.h"
#include "lex.h"
#include "block.h"
#include "mapio.h"
#include "map2gd.h"
#include "plot.h"
#include "gls.h"

static void generate_grid(DATA *d, double samplespacing, int neighbourhood);
void ossfim2map(double **table, const char *name, double s, double S,
	double b, double B, int dx, int dy);

static void generate_grid(DATA *d, double samplespacing, int neighbourhood) {
	DPOINT pt;
	int i, j, n;

	d->n_list = d->n_sel = 0;
	qtree_free(d->qtree_root);
	d->qtree_root = NULL;
	pt.z = 0.0;
	pt.attr = 0.0;
	n = floor(sqrt(neighbourhood)) / 2 + 1;
	d->maxX = d->maxY = n * samplespacing;
	d->minX = d->minY = -n * samplespacing;
	gl_split = 4 * n * n;
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			pt.x = (i + 0.5) * samplespacing;
			pt.y = (j + 0.5) * samplespacing;
			push_point(d, &pt); /* NE */
			pt.y *= -1;
			push_point(d, &pt); /* SE */
			pt.x *= -1;
			push_point(d, &pt); /* SW */
			pt.y *= -1;
			push_point(d, &pt); /* NE */
		}
	}
	d->sel_max = neighbourhood;
	return;
}

/*
	n:m:B:b:S:s:V:v:x:y:"
	-V s variogram model to print (GIF, on stdout!)\n\ 
*/

#define OSSFIM_HELPSTRING "\
gstat/ossfim options:\n\
	-n # number of nearest neighbours\n\
	-s $ minimum sample spacing\n\
	-S $ maximum sample spacing\n\
	-x # number of sample spacings to evaluate\n\
	-b $ minimum block size\n\
	-B $ maximum block size\n\
	-y # number of block sizes to evaluate\n\
	-v s variogram model\n\
	-m s 'map' name (PNG file written to)\n\
	-h   print this help\n\
[ # is a number, $ a floating point value, s a string ]\n"

int ossfim(int argc, char *argv[]) {
	int c, n = 25, dx = 9, dy = 9, i, j, plot_vgm = 0;
	double b = 1, B = 10, s = 1, S = 10, blocksize, samplespacing, est[2],
		**table;
	DATA **d = NULL;
	DPOINT *block = NULL, where;
	char *vgm_str = "1 Exp(10)", *map_name = NULL;
	VARIOGRAM *vgm;

	while ((c = getopt(argc, argv, "n:m:B:b:S:s:V:v:x:y:h")) != EOF) {
		switch (c) {
			case 'h':
				printf("%s", OSSFIM_HELPSTRING);
				exit(0);
				break;
			case 'n':
				if (read_int(optarg, &n) || n <= 0)
					ErrMsg(ER_ARGOPT, "n");
				break;
			case 'b':
				if (read_double(optarg, &b) || b < 0)
					ErrMsg(ER_ARGOPT, "b");
				break;
			case 'B':
				if (read_double(optarg, &B) || B <= 0)
					ErrMsg(ER_ARGOPT, "B");
				break;
			case 's':
				if (read_double(optarg, &s) || s <= 0)
					ErrMsg(ER_ARGOPT, "s");
				break;
			case 'S':
				if (read_double(optarg, &S) || S <= 0)
					ErrMsg(ER_ARGOPT, "S");
				break;
			case 'x':
				if (read_int(optarg, &dx) || dx <= 0)
					ErrMsg(ER_ARGOPT, "x");
				break;
			case 'y':
				if (read_int(optarg, &dy) || dy <= 0)
					ErrMsg(ER_ARGOPT, "y");
				break;
			case 'v':
				vgm_str = optarg;
				break;
			case 'V':
				plot_vgm = 1;
				vgm_str = optarg;
				break;
			case 'm':
				map_name = optarg;
				break;
			default:
				ErrClo(optopt);
				break;
		}
	}

	which_identifier("dummy grid");
	d = get_gstat_data();
	init_one_data(d[0]);
	d[0]->id = 0;
	d[0]->n_list = d[0]->n_max = 0;
	d[0]->mode = X_BIT_SET | Y_BIT_SET | V_BIT_SET;
	set_norm_fns(d[0]);
	vgm = get_vgm(0);
	if (read_variogram(vgm, vgm_str))
		ErrMsg(ER_SYNTAX, vgm_str);
	vgm->ev->evt = SEMIVARIOGRAM;
	vgm->id1 = vgm->id2 = d[0]->id;
	block = get_block_p();
	block->z = 0.0;
	block->x = block->y = -1.0;
	est[0] = 0.0;
	est[1] = -1.0;
	where.x = where.y = where.z = 0.0;
	where.X = (double *) emalloc(sizeof(double));
	where.X[0] = 1.0;

	if (plot_vgm)
		return fprint_gnuplot_variogram(stdout, vgm, "", GIF, 0);

	table = (double **) emalloc((dy + 1) * sizeof(double *));
	for (i = 0; i <= dy; i++)
		table[i] = (double *) emalloc((dx + 1) * sizeof(double));

	/* do it: */
	for (i = 0; i <= dx; i++) { /* sample spacing loop */
		samplespacing = s + (i / (1.0 * dx)) * (S - s);
		generate_grid(d[0], samplespacing, n);
		select_at(d[0], &where);
		for (j = 0; j <= dy; j++) { /* block sizes loop */
			reset_block_discr();
			vgm_init_block_values(vgm);
			blocksize = b + (j / (1.0 * dy)) * (B - b);
			block->x = block->y = blocksize;
			if (blocksize == 0.0)
				SET_POINT(&where);
			else
				SET_BLOCK(&where);
			gls(d, 1, GLS_BLUP, &where, est);
			if (map_name)
				table[i][j] = sqrt(est[1]);
			else
				printlog("%g %g %g\n", samplespacing, blocksize, sqrt(est[1]));
		}
	}
	if (map_name)
		ossfim2map(table, map_name, s, S, b, B, dx, dy);
	return 0;
}

void ossfim2map(double **table, const char *name, double s, double S,
	double b, double B, int dx, int dy) {
	
	GRIDMAP *m;
	TICKS x, y;
	char str[100];
	int i, j;
	extern int nice_legend;
	double bs;
	
	m = new_map(WRITE_ONLY);
	m->filename = name;
	m->rows = dy + 1;
	m->cols = dx + 1;
	m->cellsizex = m->cellsizey = 1.0;
	m->x_ul = m->y_ul = 1.0;
	for (i = 0; i <= dx; i++) /* row */
		for (j = 0; j <= dy; j++) /* col */
			map_put_cell(m, dy - j, i, table[i][j]); /* flips */

	/* fill ticks */
	x.every = y.every = 1;
	x.n = dx + 1;
	y.n = dy + 1;
	x.entries = (char **) emalloc(x.n * sizeof(char *));
	for (i = 0; i <= dx; i++) { /* sample spacings: */
		sprintf(str, "%3g", s + i * (S - s) / dx);
		x.entries[i] = string_dup((char *) str);
	}
	y.entries = (char **) emalloc(y.n * sizeof(char *));
	for (i = 0; i <= dy; i++) { /* sample spacings: */
		bs = b + i * (B - b) / dy;
		sprintf(str, "%g x %g", bs, bs);
		y.entries[dy - i] = string_dup((char *) str);
	}
	nice_legend = 1;
	one_map2gd(m, name, &y, &x, 1);
	return;
}
