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
 * map2fig: generate .fig figure from grid map or point data,
 * for interactive editing with xfig (or printing with fig2dev)
 */

#include <stdio.h>
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
#include "utils.h"
#include "palet.h"
#include "mapio.h"
#include "data.h"
#include "read.h"
#include "fig.h"
#include "map2fig.h"

#define DATA2FIG_OPTIONS "-c classes\n\
-l log-transform\n\
-L # legend mode (0: contiguous, 1: blocks)\n\
-n # number of classes\n\
-r f radius (?)\n\
-m s map name\n\
-p # palette number\n\
-x # x column number\n\
-y # y column number\n\
-v # attribute column number\n\
# : number; f : floating point; s : string\n\
"
#define DEFAULT_RADIUS 100 /* units 1/1200 inch */
#define MINMAX -99999.0
static void draw_map(GRIDMAP *m, int bg);
D_VECTOR *find_classes(GRIDMAP *m, DATA *d, double min, double max, 
		int nclass, int lt);
static int find_class(D_VECTOR *classes, double attr);

int map2fig(int argc, char *argv[]) {
	int c, nclass = 8, i, j, ipal, radius = DEFAULT_RADIUS, 
		logtransform = 0;
	DATA *d = NULL;
	D_VECTOR *classes = NULL;
	PALETTE pal = GYR;
	double class, min = MINMAX, max = MINMAX;
	char *cp = NULL, **entries, *map_name = NULL, *fname = NULL;
	GRIDMAP *m;
	LEGEND_MODE leg_mode = CONTIGUOUS;

	d = init_one_data(d);
	d->colnx = 1; 
	d->colny = 2; 
	d->colnvalue = 3;
	while ((c = getopt(argc, argv, "c:lL:n:m:p:r:v:x:y:z:")) != EOF) {
		switch (c) {
			case 'c':
				nclass = -1;
				while ((cp = strtok(nclass == -1 ? optarg : NULL, " \t,"))) {
					if (read_double(cp, &class))
						ErrClo(c);
					classes = push_d_vector(class, classes);
					nclass++;
				}
				break;
			case 'l': logtransform = 1; break;
			case 'L': 
				if (read_int(optarg, &(i)))
					ErrClo(c);
				switch (i) {
					case 0: leg_mode = CONTIGUOUS; break;
					case 1: leg_mode = BLOCKS; break;
					default: ErrClo(c); break;
				}
				break;
			case 'n':
				if (read_int(optarg, &(nclass)) || nclass <= 1)
					ErrClo(c);
				break;
			case 'm': map_name = optarg; break;
			case 'r':
				if (read_int(optarg, &(radius)) || radius < 1)
					ErrClo(c);
				break;
			case 'p':
				if (read_int(optarg, &(ipal)) || ipal < 0)
					ErrClo(c);
				pal = int2pal(ipal);
				break;
			case 'x': 
				if (read_int(optarg, &(d->colnx)))
					ErrClo(c);
				break;
			case 'y':
				if (read_int(optarg, &(d->colny)))
					ErrClo(c);
				break;
			case 'v':
				if (read_int(optarg, &(d->colnvalue)))
					ErrClo(c);
				break;
		}
	}

    if (argc - optind < 1) {
        printf("map2fig: Copyright 2000 Edzer J. Pebesma\n");
        printf("usage: map2fig [options] {map|data} fig_file\n");
        printf("options:\n%s", DATA2FIG_OPTIONS);
        return 0;
    }

	m = new_map(READ_ONLY);
	m->filename = argv[optind];
	if (map_read(m) == NULL) { /* fig data: */
		map_free(m);
		m = NULL;
		d->fname = argv[optind];
		d->id = 0;
		d = read_gstat_data(d);
	} 

	if (classes == NULL) /* find linear classes from min to max */
		classes = find_classes(m, d, min, max, nclass, logtransform);

	/*
	for (i = 0; i < classes->size; i++)
		printf("%g%s", classes->val[i], i < classes->size - 1 ? "," : "\n");
	*/

	/* classify data: */
	if (m == NULL) {
		for (i = 0; i < d->n_list; i++)
			d->list[i]->attr = find_class(classes, d->list[i]->attr);
	} else {
		for (i = 0; i < m->rows; i++) {
			for (j = 0; j < m->cols; j++) {
				if (!map_cell_is_mv(m, i, j))
					map_put_cell(m, i, j, 
							find_class(classes, map_get_cell(m, i, j)));
			}
		}
	}

    if (argc - optind == 2)
		fname = argv[optind + 1];

	if (m == NULL)
		fig_setup(fname, d->minX, d->minY, d->maxX, d->maxY);
	else
		fig_setup(fname, m->x_ul, m->y_ul - m->rows * m->cellsizey, 
			m->x_ul + m->cols * m->cellsizex, m->y_ul);

	for (i = 0; i < nclass; i++)
		fig_color(i, GetRGB(pal, (1.0 * i)/(nclass - 1), 255));

	if (m == NULL) {
		for (i = 0; i < d->n_list; i++)
			fig_point((int) d->list[i]->attr, 
					d->list[i]->x, d->list[i]->y, radius);
	} else
		draw_map(m, 0);

	if (map_name) {
		m = new_map(READ_ONLY);
		m->filename = map_name;
		draw_map(map_read(m), 1);
	} 

	entries = (char **) emalloc((nclass + 1) * sizeof(char *));
	if (leg_mode > 0) {
		for (i = 0; i < nclass; i++) {
			entries[i] = emalloc(50 * sizeof(char));
			sprintf(entries[i], "%g - %g", classes->val[i], classes->val[i+1]);
			printf("## [%s]\n", entries[i]);
		}
	} else {
		for (i = 0; i < nclass + 1; i++) {
			entries[i] = emalloc(50 * sizeof(char));
			sprintf(entries[i], "%g", classes->val[i]);
			printf("## [%s]\n", entries[i]);
		}
	}

	fig_legend(entries, nclass, leg_mode);
	fig_end();
	return 0;
}

D_VECTOR *find_classes(GRIDMAP *m, DATA *d, double min, double max, 
		int nclass, int lt) {
	int i, j;
	D_VECTOR *classes = NULL;

	if (m == NULL) {
		if (min == MINMAX)
			min = d->list[0]->attr;
		if (max == MINMAX)
			max = d->list[0]->attr;
		for (i = 0; i < d->n_list; i++) {
			min = MIN(min, d->list[i]->attr);
			max = MAX(max, d->list[i]->attr);
		}
	} else {
		for (i = 0; i < m->rows; i++) {
			for (j = 0; j < m->cols; j++) {
				if (!map_cell_is_mv(m, i, j)) {
					if (min == MINMAX)
						min = map_get_cell(m, i, j);
					if (max == MINMAX)
						max = map_get_cell(m, i, j);
					min = MIN(min, map_get_cell(m, i, j));
					max = MAX(max, map_get_cell(m, i, j));
				}
			}
		}
	}
	if (lt) {
		min = log(min);
		max = log(max);
		for (i = 0; i <= nclass; i++)
			classes = push_d_vector(
				exp(min + ((max - min)* i)/nclass), classes);
	} else
		for (i = 0; i <= nclass; i++)
			classes = push_d_vector(
					min + ((max - min) * i)/nclass, classes);
	return classes;
}

static int find_class(D_VECTOR *classes, double attr) {
/* returns -1 if outside class boundaries; else class entry */
	int i, sz;

	assert(classes != NULL && classes->size > 0);

	sz = classes->size;

	if (attr < classes->val[0] || attr > classes->val[sz-1])
		return -1;
		
	for (i = 1; i < sz; i++) {
		if (attr < classes->val[i] || (i == sz - 1 && attr == classes->val[i]))
			return (i - 1);
	}
	assert(0);
}

static void draw_map(GRIDMAP *m, int bg) {
	int i, j, rl;
	double x, y, val;

	assert(m);

	for (i = 0; i < m->rows; i++) {
		for (j = 0; j < m->cols; j++) {
			if (! map_cell_is_mv(m, i, j)) {
				map_rowcol2xy(m, i, j, &x, &y);
				val = map_get_cell(m, i, j);
				rl = 1;
				while (j + rl < m->cols && 
						!map_cell_is_mv(m, i, j+rl) &&
						map_get_cell(m, i, j+rl) == val)
					rl++;
				fig_run(bg ? -1 : (int) val, x, y, m->cellsizex, 
						m->cellsizey, rl);
				j += rl - 1;
			}
		}
	}
}
