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
 * fig.c: drivers for generating fig (xfig) figures (for map2fig.c)
 * xfig: ftp://ftp.x.org/contrib/applications/drawing_tools/xfig/
 */

#include <stdio.h>
#include <math.h>

#include "defs.h"
#include "userio.h"
#include "utils.h"
#include "palet.h"
#include "mapio.h"
#include "glvars.h"
#include "data.h"
#include "fig.h"

#define FIG_HEADER "#FIG 3.2\n\
Landscape\n\
Center\n\
Metric\n\
A4\n\
100.00\n\
Single\n\
-2\n\
1200 2\n\
0 32 #c0c0c0\n"

static FILE *out = NULL;
static double Xmin, Xmax, Ymin, Ymax;
static int nesting_level = 0;

#define FIGXMIN 2250
#define FIGXMAX 6750 
#define FIGYMIN 2250
#define FIGYMAX 6750 

#define LEGXMIN 6975
#define LEGXMAX 9000
#define LEGYMIN 2250
#define LEGYMAX 5400

#define LEGSIZX 495
#define LEGSIZY 225

#define LEGSKIP 135 /* between y blocks */
#define TEXTSKIP 90 /* between block and text */

#define FONT_POINT 12
#define BG_DEPTH   900
#define COLOR_OFFSET 33

#define TOFIGX(x) \
	((int)(FIGXMIN + (FIGXMAX - FIGXMIN) * (x - Xmin)/(Xmax - Xmin)))
#define TOFIGY(y) \
	((int)(FIGYMIN + (FIGYMAX - FIGYMIN) * (Ymax - y)/(Ymax - Ymin)))

void fig_setup(char *fname, double xmin, double ymin, double xmax, double ymax) {
	
	if (fname)
		out = efopen(fname, "w");
	else
		out = stdout;
	fprintf(out, "%s", FIG_HEADER);
	/* save: */
	if ((xmax - xmin) > (ymax - ymin)) { /* fit x coordinates to FIG: */
		Xmin = xmin;
		Xmax = xmax;
		Ymin = ymin;
		Ymax = ymin + (xmax - xmin);
	} else {
		Ymin = ymin;
		Ymax = ymax;
		Xmin = xmin;
		Xmax = xmin + (ymax - ymin);
	}
	printf("## x[%g,%g],y[%g,%g]\n",Xmin,Xmax,Ymin,Ymax);
}

void fig_point(int color, double x, double y, int radius) {
/*
1 3 0 1 0 33 100 0 20 0.000 1 0.0000 3780 3915 186 186 3780 3915 3960 3960
1 3 0 1 0 4  100 0 20 0.000 1 0.000     0    0 150 150    0    0  150    0
4th 1: line thinckness 
*/
	if (nesting_level == 0) {
		fprintf(out, "6 %d %d %d %d\n", FIGXMAX, FIGYMIN, FIGXMIN, FIGYMAX);
		nesting_level = 1;
	}
	fprintf(out,"1 3 0 1 0 %d 100 0 20 0.000 1 0.000 %d %d %d %d %d %d %d %d\n",
		color + COLOR_OFFSET, 
		TOFIGX(x), TOFIGY(y), radius, radius, TOFIGX(x), TOFIGY(y),
		TOFIGX(x)+radius, TOFIGY(y));
}

void fig_legend(char *entries[], int n, LEGEND_MODE leg_mode) {
/*
4 0 0 100 0 16 12 0.0000 4 135 1335 1755 5175 Some Text Here\001
16 is font; 
12 font_size; 
4->ps font (0: latex); 
135,1335: height/lenght
x,y, text
*/
	int i, yoffset = 0, leg_skip = 0;

	if (leg_mode == BLOCKS)
		leg_skip = LEGSKIP; /* the area between legend blocks */

	if (nesting_level > 0) {
		fprintf(out, "-6 # end map compound\n"); /* ends map compound */
		nesting_level--;
	}

	fprintf(out, "6 %d %d %d %d # start legend block\n",
		LEGXMAX, LEGYMIN, LEGXMIN, LEGYMAX);
	nesting_level++;

	fprintf(out, "6 %d %d %d %d # start colored legend entries\n", 
		LEGXMIN+LEGSIZX, LEGYMIN, LEGXMIN, LEGYMAX);
	nesting_level++;

	for (i = yoffset = 0; i < n; i++) {
		/* draw block: */
		fig_block(i + COLOR_OFFSET, 
			LEGXMIN, LEGYMIN + yoffset, LEGXMIN + LEGSIZX,
			LEGYMIN + yoffset + LEGSIZY, 1);
		yoffset += (LEGSIZY + leg_skip);
	}
	fprintf(out, "-6 # end colored legend entries\n");
	nesting_level--;

	fprintf(out, "6 %d %d %d %d # start legend text\n", 
		LEGXMAX, LEGYMIN, LEGXMIN + LEGSIZX, LEGYMAX); /* start text compound */
	nesting_level++;

	switch (leg_mode) {
		case BLOCKS:
			for (i = yoffset = 0; i < n; i++) {
				/* draw text: */
				fprintf(out, 
					"4 0 0 100 0 16 %d 0.0000 4 135 1335 %d %d %s\\001\n",
					FONT_POINT, LEGXMIN + LEGSIZX + TEXTSKIP, 
					LEGYMIN + LEGSIZY + yoffset, entries[i]);
				yoffset += (LEGSIZY + leg_skip);
			}
		break;
		case CONTIGUOUS:
			for (i = yoffset = 0; i <= n; i++) {
				/* draw text: */
				fprintf(out, 
					"4 0 0 100 0 16 %d 0.0000 4 135 1335 %d %d %s\\001\n",
					FONT_POINT, LEGXMIN + LEGSIZX + TEXTSKIP, 
					LEGYMIN + yoffset + 70,  
					/* 70 to align text with legend block top line */
					entries[i]);
				yoffset += LEGSIZY;
		}
		break;
	}
	fprintf(out, "-6 # end legend text\n");
	nesting_level--;

	fprintf(out, "-6 # end legend compound\n");
	nesting_level--;
}

void fig_run(int color, double x0, double y0, double cellsizex, double cellsizey, int rl) {
	assert(rl >= 1); /* rl == 1 means one cell */
	if (nesting_level == 0) {
		fprintf(out, "6 %d %d %d %d\n", FIGXMAX, FIGYMIN, FIGXMIN, FIGYMAX);
		nesting_level = 1;
	}
	fig_block(color + COLOR_OFFSET, 
		TOFIGX(x0 - 0.5 * cellsizex),
		TOFIGY(y0 - 0.5 * cellsizey),
		TOFIGX(x0 + (rl - 0.5) * cellsizex),
		TOFIGY(y0 + 0.5 * cellsizey), 
		0);
}

void fig_block(int color, int xmin, int ymin, int xmax, int ymax, int line) {
/*
2 2 0 0 0 32 100 0 20 0.000 0 0 -1 0 0 5
	 1710 3555 3195 3555 3195 4590 1710 4590 1710 3555
*/
	if (color == 32)
		fprintf(out, "2 2 0 %d 0 %d %d 0 20 0.000 0 0 -1 0 0 5\n", 
			line, color, BG_DEPTH);
	else
		fprintf(out, "2 2 0 %d 0 %d 100 0 20 0.000 0 0 -1 0 0 5\n", 
			line, color);
	fprintf(out, "\t  %d %d %d %d %d %d %d %d %d %d\n",
		xmin, ymin, xmax, ymin, xmax, ymax, xmin, ymax, xmin, ymin);
}

void fig_color(int color, float *col) {
	fprintf(out, "0 %d #%02x%02x%02x\n", 
		color + COLOR_OFFSET,
		(int) floor(col[0]), 
		(int) floor(col[1]),
		(int) floor(col[2]));
}

void fig_end(void) {
	while (nesting_level != 0)
		fprintf(out, "-6 # end <unknown> compound\n");
	fprintf(out, "## %s\n", command_line);
	if (out != stdout)
		efclose(out);
}
