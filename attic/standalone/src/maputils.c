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
#include <math.h> 
#include <float.h>
#include <string.h>
#include <stdlib.h>

#include "defs.h"

#ifdef HAVE_GETOPT_H
# include <getopt.h>
#endif
#ifdef HAVE_UNISTD_H
# include <unistd.h> /* isatty() */
#endif
#ifndef HAVE_GETOPT
# include "getopt.h"
#endif

#ifdef HAVE_LIBCSF
#include "csf.h"
#endif

#include "defs.h"
#include "glvars.h"
#include "utils.h"
#include "debug.h"
#include "read.h"
#include "userio.h"
#include "mapio.h"
#include "stat.h"
#include "maputils.h"

static void convert_help(char *prog);

int map_nominal(int argc, char *argv[]) {
	GRIDMAP *in = NULL, *out = NULL;
	int i, j, k;
	float value;

	if (argc < 4) {
		printf("usage: %s out_map in_map1 in_map2 ...\n", argv[0]);
		exit(0);
	}
	for (k = 2; k < argc; k++) {
		in = new_map(READ_ONLY);
		in->filename = argv[k];
		if ((in = map_read(in)) == NULL) {
			printf("cannot read map %s\n", argv[k]);
			exit(1);
		}
		if (out == NULL)
			out = map_dup(argv[1], in);
		for (i = 0; i < in->rows; i++) {
			for (j = 0; j < in->cols; j++) {
				if (!map_cell_is_mv(in, i, j)) {
					if (map_cell_is_mv(out, i, j)) /* initialize */
						map_put_cell(out, i, j, 1.0 * (argc - 2));
					value = map_get_cell(in, i, j);
					if (value > 0.0) 
						map_put_cell(out, i, j, 1.0 * (k - 2)); /* starting at 0 */
				}
			}
		}
		map_free(in);
	}
	out->write(out); /* does the error handling */
	return 0;
}

int map_cover(int argc, char *argv[]) {

	GRIDMAP *in = NULL, *out = NULL;
	int i, j, k;
	float value;

	if (argc < 4) {
		printf("usage: %s out_map in_map1 in_map2 ...\n", argv[0]);
		exit(0);
	}
	for (k = 2; k < argc; k++) {
		in = new_map(READ_ONLY);
		in->filename = argv[k];
		if ((in = map_read(in)) == NULL) {
			printf("cannot read map %s\n", argv[k]);
			exit(1);
		}
		if (out == NULL)
			out = map_dup(argv[1], in);
		for (i = 0; i < in->rows; i++) {
			for (j = 0; j < in->cols; j++) {
				if (map_cell_is_mv(out, i, j) && !map_cell_is_mv(in, i, j)) {
					value = map_get_cell(in, i, j);
					map_put_cell(out, i, j, value);
				}
			}
		}
		map_free(in);
	}
	out->write(out); /* does the error handling */
	return 0;
}

int map_cut(int argc, char *argv[]) {
/* cut a square region from a map */

	GRIDMAP *in = NULL, *out = NULL;
	int i, j, from_row, to_row, from_col, to_col;
	float value;

	if (argc != 7) {
		printf("usage: %s from_row to_row from_col to_col in_map out_map\n", argv[0]);
		exit(0);
	}
	if (read_int(argv[1], &from_row))
		ErrMsg(ER_RDINT, argv[1]);
	if (read_int(argv[2], &to_row))
		ErrMsg(ER_RDINT, argv[2]);
	if (read_int(argv[3], &from_col))
		ErrMsg(ER_RDINT, argv[3]);
	if (read_int(argv[4], &to_col))
		ErrMsg(ER_RDINT, argv[4]);
	from_row--; from_col--; to_row--; to_col--; /* start counting at 0 */

	in = new_map(READ_ONLY);
	in->filename = argv[5];
	if ((in = map_read(in)) == NULL) {
		printf("cannot read map %s\n", argv[5]);
		exit(1);
	}
	if (from_col < 0 || from_col > in->cols-1)
		ErrMsg(ER_IMPOSVAL, "from_col");
	if (to_col < 0 || to_col > in->cols-1)
		ErrMsg(ER_IMPOSVAL, "to_col");
	if (from_row < 0 || from_row > in->rows-1)
		ErrMsg(ER_IMPOSVAL, "from_row");
	if (to_row < 0 || to_row > in->rows-1)
		ErrMsg(ER_IMPOSVAL, "to_row");

	if (in->CSF_MAP) {
#ifdef HAVE_LIBCSF
		out = new_map(READ_ONLY);
		out->CSF_MAP = Rcreate(argv[6],
			to_row - from_row + 1, to_col - from_col + 1,
			RgetCellRepr(in->CSF_MAP), RgetValueScale(in->CSF_MAP), 
			MgetProjection(in->CSF_MAP), in->x_ul + from_col * in->cellsizex,
			MgetProjection(in->CSF_MAP) == PT_YDECT2B ?
				in->y_ul - from_row * in->cellsizey :
				in->y_ul + from_row * in->cellsizey,
			RgetAngle(in->CSF_MAP), in->cellsizex);
		out->type = in->type;
#endif
	} else
		out = map_dup(argv[6], in); /* large enough... */
	out->rows = to_row - from_row + 1;
	out->cols = to_col - from_col + 1;
	out->x_ul = in->x_ul + from_col * in->cellsizex;
	out->y_ul = in->y_ul - from_row * in->cellsizey;
	for (i = from_row; i <= to_row; i++) {
		for (j = from_col; j <= to_col; j++) {
			if (!map_cell_is_mv(in, i, j)) {
				value = map_get_cell(in, i, j);
				map_put_cell(out, i - from_row, j - from_col, value);
			}
		}
	}
	map_free(in);
	out->write(out); /* does the error handling */
	return 0;
}

int map_diff(int argc, char *argv[]) {
	GRIDMAP *a = NULL, *b = NULL;
	int i, j, n_diff = 0;
	float val_a, val_b;

	if (argc != 3) {
		printf("usage: %s map_a map_b\n", argv[0]);
		return(1);
	}
	a = new_map(READ_ONLY);
	a->filename = argv[1];
	if ((a = map_read(a)) == NULL) {
		printf("cannot read map %s\n", argv[1]);
		return(1);
	}
	b = new_map(READ_ONLY);
	b->filename = argv[2];
	if ((b = map_read(b)) == NULL) {
		printf("cannot read map %s\n", argv[2]);
		map_free(a);
		return(1);
	}
	printf("%s and %s ", argv[1], argv[2]);
	if (! map_equal(a,b)) {
		printf("differ in map topology\n");
		map_free(a); map_free(b);
		return(1);
	}
	for (i = 0; i < a->rows; i++) {
		for (j = 0; j < a->cols; j++) {
			if (map_cell_is_mv(a, i, j) != map_cell_is_mv(b, i, j)) {
				if (n_diff == 0)
					printf("[%d,%d]: missing values do not match ", i+1, j+1);
				n_diff++;
			} else if (!map_cell_is_mv(a, i, j)) {
				val_a = map_get_cell(a, i, j);
				val_b = map_get_cell(b, i, j);
				if (val_a != val_b) { 
					if ((fabs(val_a - val_b)) > gl_zero) {
						if (n_diff == 0)
							printf("[%d,%d]: %g != %g ", i+1, j+1, val_a, val_b);
						n_diff++;
					}
				}
			}
		}
	}
	if (n_diff)
		printf("(%d cells differ)\n", n_diff);
	else 
		printf("are equal\n");
	map_free(a);
	map_free(b);
	return n_diff;
}

int map_lnh(int argc, char *argv[]) {

	GRIDMAP *m2s, *p2s, *out;
	float m, p, level;
	int i, j;
	int ans;

#ifdef HAVE_LIBCSF
	CSF_LEGEND *buf;
	char *cp;
#endif
	
	if (argc != 5) {
		printf("usage: %s m2s p2s out level\n", argv[0]);
		exit(0);
	}
	m2s = new_map(READ_ONLY);
	m2s->filename = argv[1];
    if ((m2s = map_read(m2s)) == NULL) {
        ErrMsg(ER_READ, argv[1]);
        exit(1);
    }
	p2s = new_map(READ_ONLY);
	p2s->filename = argv[2];
    if ((p2s = map_read(p2s)) == NULL) {
        ErrMsg(ER_READ, argv[2]);
        exit(1);
    }
	m2s->celltype = CT_UINT8;
	if ((out = map_dup(argv[3], m2s)) == NULL) {
		ErrMsg(ER_READ, argv[3]);
		exit(1);
	}

	level = atof(argv[4]);

    for (i = 0; i < m2s->rows; i++) {
        for (j = 0; j < m2s->cols; j++) {
            if (!(map_cell_is_mv(m2s, i, j) || map_cell_is_mv(p2s, i, j))) {
				m = map_get_cell(m2s, i, j);
				p = map_get_cell(p2s, i, j);
				if (p < m)
					printf("%s: %g < %g -- a bit funny that is!\n", 
						argv[0], p, m);
				if (p < level)
					ans = 0.0;
				else if (m > level)
					ans = 2.0;
				else
					ans = 1.0;
				map_put_cell(out, i, j, ans);
			} 
		}
	}

#ifdef HAVE_LIBCSF
	if (out->CSF_MAP) {
    	buf = (CSF_LEGEND *) malloc (4 * sizeof(CSF_LEGEND));

		/* buf[0].nr = -1; */
		/* name: */
		cp = strrchr(argv[2], '.');
		if (cp)
			*cp = '\0';
		sprintf(buf[0].descr, "%s [%g]", argv[2], level);
		buf[1].nr = 0;
		sprintf(buf[1].descr, "lower");
		buf[2].nr = 1;
		sprintf(buf[2].descr, "not distinguishable");
		buf[3].nr = 2;
		sprintf(buf[3].descr, "higher");

		MputLegend(out->CSF_MAP, buf, 4);
		RputValueScale(out->CSF_MAP, VS_NOMINAL);
	}
#endif
	out->write(out);
	map_free(out);
	return 0;
}

int map_q(int argc, char *argv[]) {
	int i, j, k, n;
	double *stack, q;
	GRIDMAP **in, *out;

	if (argc < 4) {
		printf("usage: %s <q> out_map in_map ...\n", argv[0]);
		exit(0);
	}
	if (read_double(argv[1], &q))
		ErrMsg(ER_RDFLT, argv[1]);
	printlog("q: %g\n", q);
	n = argc - 3;
	in = (GRIDMAP **) emalloc(n * sizeof(GRIDMAP *));
	stack = (double *) emalloc(n * sizeof(double));
	for (i = 0; i < n; i++) {
		in[i] = new_map(READ_ONLY);
		in[i]->filename = argv[i+3];
		in[i] = map_read(in[i]);
	}
	out = map_dup(argv[2], in[0]);
	for (i = 0; i < out->rows; i++) {
		for (j = 0; j < out->cols; j++) {
			if (! map_cell_is_mv(in[0], i, j)) {
				for (k = 0; k < n; k++)
					stack[k] = map_get_cell(in[k], i, j);
    			qsort(stack, (size_t) n, sizeof(double),
	        		(int CDECL (*)(const void *, const void *)) d_cmp);
				map_put_cell(out, i, j, (float) est_quant(stack, q, n));
			}
		}
	}
	out->write(out);
	map_free(out);
	for (i = 0; i < n; i++)
		map_free(in[i]);
	efree(in);
	efree(stack);
	return 0;
}

int map_convert(int argc, char *argv[]) {

	GRIDMAP *in = NULL, *out = NULL;
	MAPTYPE t = MT_UNKNOWN;
	unsigned int i, j;
	int c;
	char format = '-';

	while ((c = getopt(argc, argv, "f:h")) != EOF) {
		switch (c) {
			case 'f': format = *optarg; break;
			case 'h': convert_help(argv[0]); break;
		}
	}

	if (argc - optind < 2) { /* no two arguments left */
		convert_help(argv[0]);
		exit(0);
	}
	in = new_map(READ_ONLY);
	in->filename = argv[optind];
	if ((in = map_read(in)) == NULL)
		ErrMsg(ER_READ, argv[optind]);

#ifdef HAVE_LIBCSF
	if (in->type == MT_CSF)
		Mclose(in->CSF_MAP);
#endif

	switch (format) {
		case 'a': t = MT_ARCGRID; in->is_binary = 0; break;
		case 'b': t = MT_IDRISI; in->is_binary = 1; break;
		case 'c': t = MT_CSF; break;
		case 'e': t = MT_ERMAPPER; break;
		case 'd': t = MT_SURFER; break;
		case 'f': 
			t = MT_ARCGRID; 
			in->is_binary = 1;
			in->celltype = CT_IEEE4;
			break;
		case 'G': t = MT_GNUPLOT; break;
		case 'g': t = MT_GSLIB; break;
		case 'i': t = MT_IDRISI; in->is_binary = 0; break;
		case 'n': t = MT_GMT; break;
		case 't': t = MT_T2; break;
		default:
			printf("%c: unknown target type\n", format);
			exit(1);
			break;
	}
	in = map_switch_type(in, t);
	out = map_dup(argv[optind + 1], in);
	/* copy all cells */
	for (i = 0; i < in->rows; i++)
		for (j = 0; j < in->cols; j++)
			if (! map_cell_is_mv(in, i, j))
				map_put_cell(out, i, j, map_get_cell(in, i, j));
				/* internal missing values should be kept NaN! */
	out->write(out); /* does error handling and NaN -> misval conversion */
	return 0;
}

static void convert_help(char *prog) {
	printf("Minimal map converter (topology and cell values only)\n");
	printf("usage: %s [options] in_map out_map\n", prog);
	printf("options:\n\t-f [type]\n\t-h   print help\n");
	printf("target types supported:\n");
	printf("\ta gridascii\n");
	printf("\tb idrisi binary\n");
#ifdef HAVE_LIBCSF
	printf("\tc PCRaster/csf\n");
#endif
	printf("\td Surfer DSAA\n");
	printf("\te er-mapper\n");
	printf("\tf gridfloat\n");
	printf("\tg GSLIB\n");
	printf("\ti idrisi ascii\n");
#ifdef HAVE_LIBNETCDF
	printf("\tn GMT/netcdf\n");
#endif
#ifdef HAVE_T2_GRIDFORMAT
	printf("\tt T2\n");
#endif
}
