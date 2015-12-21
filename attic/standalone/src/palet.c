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
for the BPY palet, thanks to it's creator, look at
http://www-ihe.etec.uni-karslruhe.de/mitarbeiter/vonhagen/palette.en.html 
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "defs.h"

#ifdef HAVE_GETOPT_H
# include <getopt.h>
#endif
#ifdef HAVE_UNISTD_H
# include<unistd.h>
#endif
#ifndef HAVE_GETOPT
# include "getopt.h"
#endif

#include "utils.h"
#include "userio.h"
#include "palet.h"

#define MAXCOL 255.0

#define PALETUSAGE "\
usage: palet -p <palet> [options]\n\
options are:\n\
-p #   palet #: 1 gyr; 2 ryg; 3 rainbow; 4 blue-pink-yellow\n\
-i     print as integers (default floats)\n\
-m #   maxcolor # (default 255)\n\
-g     print as gnuplot (png term) palet\n\
-n #   max nr of colours\n\
-c #   use CMYK scheme\n\
-u s   use s as palette (e.g., 0xff0000,0x00ff00,0x0000ff )\n\
"

#define GET_RED(x) (x >> 16)
#define GET_GRE(x) ((x >> 8) % 256)
#define GET_BLU(x) (x & 0xff)

typedef struct {
	int n;
	int *RGB;
} PAL;

static PAL *push_to_palet(PAL *p, unsigned int rgb);
static void interpolate_pal(double pos, PAL *pal, float colors[]);
static PAL *user_palette = NULL;

int palet(int argc, char *argv[]) {
	int i = 0, ipal = 0, nr = 20, c, asint = 0, asgnuplot = 0, CMYK = 0;
	unsigned int color;
	PALETTE pal = NOTSET;
	float pos, *col, maxcol = MAXCOL;
	char *cp;

	while ((c = getopt(argc, argv, "cign:m:p:u:")) != EOF) {
		switch (c) {
			case 'c': CMYK = 1; break;
			case 'i': asint = 1; break;
			case 'n': nr = atoi(optarg); break;
			case 'm': maxcol = atof(optarg); break;
			case 'g': asgnuplot = 1; maxcol = 255; break;
			case 'p': ipal = atoi(optarg); break;
			case 'u': 
				while ((cp = strtok(i == 0 ? optarg : NULL, ", ")) != NULL) {
					if (sscanf(cp, "0x%06x", &color) != 1)
						ErrClo(c);
					user_palette = push_to_palet(user_palette, color);
					/* printf("%s : %06x\n", cp, color); */
					i++;
				}
				pal = USER_DEFINED;
				break;
			default:
				ErrClo(c);
				break;
		}
	}

	if (pal == NOTSET)
		pal = int2pal(ipal);

	if (nr <= 1) {
		fprintf(stderr, "number %d too small\n", nr); 
		exit(1);
	}

	if (nr <= 1) {
		printf("nr should be > 1\n");
		exit(1);
	}

	for (i = 0; i < nr; i++) {
		pos = (1.0 * i)/(nr - 1);
		col = GetRGB(pal, pos, maxcol);
		if (asgnuplot) {
			if (i == 0)
				printf("xffffff x000000 x000000 \\\n");
			printf("x%02x%02x%02x", (int) floor(col[0]), (int) floor(col[1]),
				(int) floor(col[2]));
			if ((i + 1) % 10 == 0)
				printf("%s\n", i != (nr - 1) ? " \\" : "");
			else
				printf(" ");
		} else if (asint) {
			printf("%3d %3d %3d\n", 
				(int) col[0], (int) col[1], (int) col[2]);
		} else if (CMYK == 0)
			printf("%g %g %g\n", col[0], col[1], col[2]);
		else
			printf("%g %g %g %g\n", maxcol - col[0], 
				maxcol - col[1], maxcol - col[2], 0.0);
	}

	return 0;
}


float *GetRGB(PALETTE palet, float pos, float maxcol) {
/*
 * Return R,G and B value as int in [0,maxcol]
 * given some colorset
 */
	static float colors[3];
	int i;
	static PALETTE p = NOTSET;
	static PAL *pal = NULL;

	if (pos < 0.0 || pos > 1.0)
		printf("cannot find colour for pos %g\n", pos);
		
	if (pal == NULL || p != palet) { /* set or change pal: */
		switch (palet) {
			case RYG: 
				pal = push_to_palet(NULL, 0xff0000);
				pal = push_to_palet(pal, 0x00ff00);
				break;
			case GYR: 
				pal = push_to_palet(NULL, 0x00ff00);
				pal = push_to_palet(pal, 0xff0000);
				break;
			case UNIRAS: 
				pal = push_to_palet(NULL, 0xff00ff);
				pal = push_to_palet(pal, 0x0000ff);
				pal = push_to_palet(pal, 0x00ffff);
				pal = push_to_palet(pal, 0x00ff00);
				pal = push_to_palet(pal, 0xffff00);
				pal = push_to_palet(pal, 0xff0000);
				break;
			case BPY: 
				pal = push_to_palet(NULL, 0x000066);
				pal = push_to_palet(pal, 0x0000a3);
				pal = push_to_palet(pal, 0x0000e1);
				pal = push_to_palet(pal, 0x1900ff);
				pal = push_to_palet(pal, 0x4900ff);
				pal = push_to_palet(pal, 0x7a00ff);
				pal = push_to_palet(pal, 0xaa16e8);
				pal = push_to_palet(pal, 0xdb35c9);
				pal = push_to_palet(pal, 0xff54aa);
				pal = push_to_palet(pal, 0xff738b);
				pal = push_to_palet(pal, 0xff926c);
				pal = push_to_palet(pal, 0xffb14d);
				pal = push_to_palet(pal, 0xffd02e);
				pal = push_to_palet(pal, 0xffef0f);
				pal = push_to_palet(pal, 0xffff5f);
				break;
			case GYR0: 
				pal = push_to_palet(NULL,0x197f00);
				pal = push_to_palet(pal, 0x339900);
				pal = push_to_palet(pal, 0x59ad00);
				pal = push_to_palet(pal, 0x7fbf00);
				pal = push_to_palet(pal, 0x99cc00);
				pal = push_to_palet(pal, 0xb2d800);
				pal = push_to_palet(pal, 0xd8ed00);
				pal = push_to_palet(pal, 0xf2ff00);
				pal = push_to_palet(pal, 0xfad600);
				pal = push_to_palet(pal, 0xffbf00);
				pal = push_to_palet(pal, 0xffa800);
				pal = push_to_palet(pal, 0xff7f00);
				pal = push_to_palet(pal, 0xff6600);
				pal = push_to_palet(pal, 0xff5400);
				pal = push_to_palet(pal, 0xff0000);
				break;
			case BW:
			case WB:
				/* pal = push_to_palet(NULL, 0xffffff); */
				pal = push_to_palet(NULL, 0xacacac);
				pal = push_to_palet(pal, 0x8c8c8c);
				pal = push_to_palet(pal, 0x6c6c6c);
				pal = push_to_palet(pal, 0x4c4c4c);
				pal = push_to_palet(pal, 0x323232);
				pal = push_to_palet(pal, 0x191919);
				pal = push_to_palet(pal, 0x0a0a0a);
				/* pal = push_to_palet(pal, 0x000000); */
				break;
			case USER_DEFINED:
				pal = user_palette;
				break;
			default:
				printf("unknown color set: %d\n", palet);
				exit(1);
				break;
		}
		p = palet;
		/*
		for (i = 0; i < pal->n; i++)
			printf("%06x\n", pal->RGB[i]);
		*/
	}

	interpolate_pal(pos, pal, colors);

	for (i = 0; i < 3; i++)
		colors[i] = maxcol * colors[i]/255.0;

	if (colors[0] < 0.0)
		colors[0] = 0.0;
	if (colors[1] < 0.0)
		colors[1] = 0.0;
	if (colors[2] < 0.0)
		colors[2] = 0.0;
	return colors;
}

PALETTE int2pal(int ipal) {
	PALETTE pal;

	switch (ipal) {
		case 1: pal = GYR; break;
		case 2: pal = RYG; break;
		case 3: pal = UNIRAS; break;
		case 4: pal = BPY; break;
		case 5: pal = GYR0; break;
		case 6: pal = BW; break;
		case 7: pal = WB; break;
		default: 
			fprintf(stderr, "usage: %s\n", PALETUSAGE); 
			exit(0);
			break;
	}
	return pal;
}

static void interpolate_pal(double pos, PAL *pal, float colors[3]) {
	int nc, offset, r1, r2, g1, g2, b1, b2;

	nc = pal->n;
	if (pos == 1.0) {
		colors[0] = GET_RED(pal->RGB[nc-1]);
		colors[1] = GET_GRE(pal->RGB[nc-1]);
		colors[2] = GET_BLU(pal->RGB[nc-1]);
		return;
	}
	/*
	if (pos == 0.0) {
		colors[0] = GET_RED(pal->RGB[0]);
		colors[1] = GET_GRE(pal->RGB[0]);
		colors[2] = GET_BLU(pal->RGB[0]);
		return;
	}
	*/
	offset = (int) floor(pos * (nc - 1));
	assert(offset <= nc - 1);
	pos = pos * (nc - 1) - offset; 
	/* now pos measures between two colour indexes */
	r1 = GET_RED(pal->RGB[offset]);
	g1 = GET_GRE(pal->RGB[offset]);
	b1 = GET_BLU(pal->RGB[offset]);
	r2 = GET_RED(pal->RGB[offset+1]);
	g2 = GET_GRE(pal->RGB[offset+1]);
	b2 = GET_BLU(pal->RGB[offset+1]);
	/* printf("%d %d %d %d %d %d : %g\n", r1, r2, g1, g2, b1, b2, pos); */
	colors[0] = pos * r2 + (1.0 - pos) * r1;
	colors[1] = pos * g2 + (1.0 - pos) * g1;
	colors[2] = pos * b2 + (1.0 - pos) * b1;
	return;
}

static PAL *push_to_palet(PAL *p, unsigned int rgb) {
	if (p == NULL) {
		p = (PAL *) emalloc(sizeof(PAL));
		p->n = 0;
		p->RGB = (int *) emalloc(10 * sizeof(int));
	}
	p->RGB = (int *) erealloc(p->RGB, (p->n+1) * sizeof(unsigned int));
	p->RGB[p->n] = rgb;
	p->n++;
	return p;
}
