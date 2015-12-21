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
#include <string.h>
#include <stdlib.h>
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

#include "utils.h"
#include "mapio.h"
#include "userio.h"
#include "palet.h"
#include "read.h"
#include "map2gd.h"

/*
External (wrt. gstat) uses:
libgd graphics library by Thomas Boutell, www.boutell.com.
   v. 1.5 and lower support gif; v. 1.6 and up only support png
   Get libraray and documentation from http://www.boutell.com/gd/
png support requires additionally:
- libpng (png library) http://www.cdrom.com/pub/png
- libz (compression lib used by libpng) http://www.cdrom.com/pub/infozip/zlib
*/

#define OPTIONS "\
-c f  colour table file f (%r %g %b on each line)\n\
-i    make interlaced gif/png\n\
-m $  set minimum value to $\n\
-M $  set maximum value to $\n\
-t    no transparent background\n\
-T s  use string s as legend title (-T """" to suppress title)\n\
-p n  use n pixels per map grid cell\n\
-P #  use palet #\n\
-l n  print legend, n: 0 no legend; 1 right; 2 bottom; 3 classes\n\
-n #  number of legend classes (cont.)\n\
-a    automatic legend boundaries\n\
-b    draw box\n\
-C 1,2,3,4 use class boundaries 1,2,3,4\n\
"

#define LEGSIZEX 90
#define LEGSIZEY 45
#define TICKSHEIGHT 15
#define TICKSWIDTH 35

#define MAXCOL 255

#define TABLE_SIZE 1000

#ifdef HAVE_LIBCSF
# include "csf.h"
# define IS_CLASSIFIED(m) (legend == 3 || \
	(m->CSF_MAP && RgetValueScale(m->CSF_MAP) != VS_SCALAR && \
	RgetValueScale(m->CSF_MAP) != VS_CONTINUOUS))
#else
 #define IS_CLASSIFIED(m) (legend == 3)
#endif
#define COLOUR_TSIZE 16

int px = 0, interlace = 0, transparent = 1,
	min_set = 0, max_set = 0, table_read = 0, nc = 0, nice_legend = 0,
	nc_max = 0, n_table = 0, nclass = 16, draw_map_box = 0, ticks = 0;
float min, max, interval, *table; 
char colour_table[256] = "", *title = NULL, tmp_str[256] = "";
PALETTE pal = GYR0;

#ifdef STANDALONE
# define map2gd main
#endif

#ifndef HAVE_LIBGD
int map2gd(int argc, char *argv[]) {
	printlog("for %s, install the gd library (http://www.boutell.com/gd/)\n", argv[0]);
	printlog("(note that gd versions below 1.6 support GIF, others only PNG)\n");
	printlog("then, use  ./configure --with-gd; make; make install\n");
	ErrMsg(ER_ARGOPT, argv[0]);
	return 0;
}

int one_map2gd(GRIDMAP *m, const char *f_name, TICKS *xt, TICKS *yt, 
		int legend) {
	char *argv = "ossfim gif/png output";
	return map2gd(1, &argv);
}
#else

#include "gd.h"

static int do_one_map2gd(char *map, const char *out_name, int legend);
static void draw_one_cell(gdImagePtr g, int x, int y, int pixx, int pixy, 
		int color);
static double find_reasonable_hash_start(double min, double max, 
		double hash_interval);
static double find_reasonable_hash_interval(double range);
static int read_table(char *s);
static int find_entry(float *table, int n_table, float value);

/* 
 * convert grid map to gif/png image 
 */
int map2gd(int argc, char *argv[]) {
	int c, legend = 0, ipal = 0;
	char *in = NULL, *out = "-";

#if (!HAVE_GDIMAGEPNG) && (!HAVE_GDIMAGEGIF) /* with gd versions before libgd 1.6.1+ */
# error gdImagePng/Gif not available -- get gd from http://www.boutell.com/gd/
#endif

#ifdef HAVE_GDIMAGEPNG
	if (strcmp(argv[0], "map2gif") == 0) {
		printlog("map2gif not supported -- only map2png will work\n");
		ErrMsg(ER_ARGOPT, argv[0]);
	}
#else /* GIF: */
	if (strcmp(argv[0], "map2png") == 0) {
		printlog("map2png not supported -- only map2gif will work\n");
		ErrMsg(ER_ARGOPT, argv[0]);
	}
#endif

	opterr = 0;
	while ((c = getopt(argc, argv, "abc:C:hil:m:M:n:p:P:tT:x")) != EOF) {
		switch (c) {
			case 'a': nice_legend = 1; break;
			case 'b': draw_map_box = 1; break;
			case 'c': sprintf(colour_table, "%s", optarg); break;
			case 'h': 
			case '?': 
				printf("%s: Copyright 1997 Edzer J. Pebesma\n", argv[0]);
				printf("usage: %s [options] in_map out_file\n", argv[0]);
				printf("options:\n%s", OPTIONS);
				return 0;
				break;
			case 'i': interlace = 1; break;
			case 'l': legend = atoi(optarg); break;
			case 'm': min = atof(optarg); min_set = 1; break;
			case 'M': max = atof(optarg); max_set = 1; break;
			case 'n': nclass = atoi(optarg); break;
			case 'p': px = atoi(optarg); break;
			case 'P': 
				if (read_int(optarg, &(ipal)) || ipal < 0)
					ErrClo(c);
				pal = int2pal(ipal);
				break;
			case 'C': nc_max = read_table(optarg); break;
			case 't': transparent = 0; break;
			case 'T': title = optarg; break;
			case 'x': ticks = 1; break;
			default: ErrClo(optopt);
		}
	}

	if (argc - optind == 0) {
		printf("%s: Copyright 1997 Edzer J. Pebesma\n", argv[0]);
		printf("usage: %s [options] in_map out_file\n", argv[0]);
		printf("options:\n%s", OPTIONS);
		return 0;
	}
	in = argv[optind];
	if (argc - optind == 2)
		out = argv[optind+1];
	if (argc - optind <= 2)
		return do_one_map2gd(in, out, legend);
	ErrMsg(ER_IMPOSVAL, "multiple gif merger not longer supported");
	return 0;
}

static int do_one_map2gd(char *map, const char *f_name, int legend) {
	GRIDMAP *m = NULL;
		
	m = new_map(READ_ONLY);
	m->filename = map;
	if ((m = map_read(m)) == NULL)
		ErrMsg(ER_READ, map);
	return one_map2gd(m, f_name, NULL, NULL, legend);
}

int one_map2gd(GRIDMAP *m, const char *f_name, TICKS *xt, TICKS *yt, 
		int legend) {
	FILE *out = NULL;
#ifdef HAVE_LIBCSF
	CSF_LEGEND *l = NULL;
#else
	typedef unsigned int UINT2;
#endif
	gdImagePtr g;
	unsigned int i, j, k, class, sx, sy;
	int gcol[256], black, si, offsetx;
	float value, R, G, B, *rgb;
	extern gdFontPtr gdFontSmall; 
	gdFontPtr font;
	UINT2 *col_buf = NULL, 
		nominal_table[3 * COLOUR_TSIZE] = {
			0,128,128,
			255,128,0,
			128,128,0,
			0,128,0,
			0,128,255,
			128,0,0,
			128,0,255,
			0,0,128,
			128,0,128,
			255,0,128,
			0,188,188,
			255,188,0,
			188,188,0,
			0,188,0,
			0,188,255,
			188,0,0 
		};

	font = gdFontSmall;
	if (xt || yt)
		ticks = 1;

	if (px < 1)
		px = MAX(1, (500 / MAX(m->cols, m->rows))); /* default max to 500 px */

#ifdef HAVE_LIBCSF
	if (m->CSF_MAP && (i = MgetNrLegendEntries(m->CSF_MAP)) > 0) {
		l = emalloc(i * sizeof(CSF_LEGEND));
		MgetLegend(m->CSF_MAP, l);
		if (legend > 0)
			legend = 3;
		title = l[0].descr;
	}
#endif
	if (title == NULL || title[0] == '\0')
		title = (char *) m->filename; /* suppress this with -T "" */

#ifdef HAVE_LIBCSF 
	/* read colour table from map, if present: */
	if (m->type == MT_CSF && m->CSF_MAP != NULL) { /* read table from map */
		if ((nc = MgetNrColourPaletteEntries(m->CSF_MAP)) > 0) {
		/* colour table is present in map: */
			col_buf = emalloc(nc * 3 * sizeof(UINT2));
			MgetColourPalette(m->CSF_MAP, col_buf);
		} else if ((nc = MgetNrGreyPaletteEntries(m->CSF_MAP)) > 0) {
			col_buf = emalloc(nc * 3 * sizeof(UINT2));
			MgetGreyPalette(m->CSF_MAP, col_buf);
			for (si = nc - 1; si >= 0; si--) { /* fill colours backwards */
				col_buf[3 * si + 2] = col_buf[si];
				col_buf[3 * si + 1] = col_buf[si];
				col_buf[3 * si] = col_buf[si];
			}
		}
		for (i = 0; i < nc; i++) {
			col_buf[i*3] = col_buf[i*3] >> 8; /* [0,65535] -> [0,255] */
			col_buf[i*3+1] = col_buf[i*3+1] >> 8;
			col_buf[i*3+2] = col_buf[i*3+2] >> 8;
		}
	}
#endif
	if (nc == 0) { /* no table read... */
		if (colour_table[0] == '\0' && getenv("PCR_DIR") != NULL)
			sprintf(colour_table, "%s%s", getenv("PCR_DIR"),
				IS_CLASSIFIED(m) ? "class.pal" : "con.pal");
		if ((out = fopen(colour_table, "r")) != NULL) { 
			/* found colour table file -- read it */
			col_buf = emalloc(256 * 3 * sizeof(UINT2));
			i = 0;
			while (fscanf(out, "%g %g %g", &R, &G, &B) == 3) {
				if (i == 256)
					ErrMsg(ER_IMPOSVAL, "colour table too big (max. 256)");
				col_buf[i * 3]     = (UINT2) floor((R * 2.55) + 0.5);
				col_buf[i * 3 + 1] = (UINT2) floor((G * 2.55) + 0.5);
				col_buf[i * 3 + 2] = (UINT2) floor((B * 2.55) + 0.5);
				i++;
			}
			fclose(out);
			nc_max = nc = i;
		} else if (IS_CLASSIFIED(m)) {
			/* use default, built-in colour table */
			col_buf = nominal_table;
			nc_max = nc = COLOUR_TSIZE;
		}
	}

	if (n_table > 0) {
		min = table[0];
		max = table[n_table-1];
		min_set = max_set = 1;
		nc = n_table - 1;
	}

	if (m->cellmax == m->cellmin) { /* prevent infinite loop */
		min = max = m->cellmin;
		min_set = max_set = 1;
		interval = 1;
		nclass = 2;
	} else { 
		if (nice_legend)
			interval = find_reasonable_hash_interval(m->cellmax - m->cellmin);
	
		if (! min_set) {
			if (nice_legend)
				min = find_reasonable_hash_start(m->cellmin, 
						m->cellmax, interval);
			else
				min = m->cellmin;
		}

		if (! max_set) {
			if (nice_legend) {
				max = min;
				nc = 0;
				while (max < m->cellmax) {
					max += interval;
					nc++;
				}
				nc_max = nc;
			} else
				max = m->cellmax;
		}
	}

	if (legend == 3 && max - min < nc)
		nc = ceil(max - min) + 1;

	if (nc > 254) {
		printf("too many colours in %s (max 254)\n", m->filename);
		exit(1);
	}

	if (col_buf == NULL) {
		if (nc == 0)
			nc_max = nc = nclass;
		col_buf = (UINT2 *) emalloc(nc * 3 * sizeof(UINT2));
		for (i = 0; i < nc; i++) {
			rgb = GetRGB(pal, i/(1.0 * (nc - 1)), 255);
			col_buf[3*i] = (UINT2) floor(rgb[0]);
			col_buf[3*i+1] = (UINT2) floor(rgb[1]);
			col_buf[3*i+2] = (UINT2) floor(rgb[2]);
		}
	}

	sx = m->cols * px + (legend && legend != 2) * LEGSIZEX + ticks * TICKSWIDTH;
	sy = m->rows * px + (legend == 2) * LEGSIZEY + ticks * TICKSHEIGHT;

	g = gdImageCreate(sx, sy);

	gcol[nc] = gdImageColorAllocate(g, 255, 255, 255); /* back ground: white */
	black = gdImageColorAllocate(g, 0, 0, 0); /* text colour */
	if (legend) {
		for (i = 0; i < nc; i++) {
			if (legend == 3)
				class = i + min;
			else
				class = i;
			if (class > nc_max)
				pr_warning("class: %d, nc_max: %d", class, nc_max);
			assert(class <= nc_max);
			gcol[i] = gdImageColorAllocate(g, col_buf[class * 3],
					col_buf[class * 3 + 1], col_buf[class * 3 + 2]);
		}
	} else {
		for (i = 0; i < nc; i++)
			gcol[i] = gdImageColorAllocate(g, col_buf[i * 3],
					col_buf[i * 3 + 1], col_buf[i * 3 + 2]);
	}

	if (transparent)
		gdImageColorTransparent(g, gcol[nc]); /* white -> transparent */

	gdImageInterlace(g, interlace);

	draw_one_cell(g, 0, 0, sx, sy, gcol[nc]); 
	/* initialize png white/transparent */

	/* draw cells: */
	for (i = 0; i < m->rows; i++) {
		for (j = 0; j < m->cols; j++) {
			if (! map_cell_is_mv(m, i, j)) {
				value = map_get_cell(m, i, j) - min; /* [0,...,max-min] */
				if (IS_CLASSIFIED(m)) {
					/* map as integers to [0,...,nc>: */
					class = floor(value);
					if (nc > 1 && class >= nc) /* OVERFLOW------>>> */
						class = class % (nc - 1); /* repeat colours */
				} else {
					/* map to [0,...,nc-1]: */
					if (table == NULL) {
						class = (max == min) ? 0 : floor(nc * value/(max-min));
						if (class == nc) /* occurs only when value == max-min */
						class--;
					} else
						class = find_entry(table, n_table, value);
				}
				if (class != -1)
					draw_one_cell(g, j * px, i * px, px, px, gcol[class]);
			} 
		}
	}
	/* draw map box */
	if (draw_map_box) {
		gdImageLine(g, 0, 0, 0, m->rows * px - 1, black);
		gdImageLine(g, 0, m->rows * px - 1, m->cols * px - 1, 
			m->rows * px - 1, black);
		gdImageLine(g, m->cols * px - 1, m->rows * px - 1, 
			m->cols * px - 1, 0, black);
		gdImageLine(g, m->cols * px - 1, 0, 0, 0, black);
	}
	switch (legend) {
		case 0: /* no legend: */ break;
		case 1:
			offsetx = m->cols * px + 10 + ticks * TICKSWIDTH;
			for (i = 0; i <= nc; i++) {
				j = px * m->rows / 2 + (i - nc / 2) * 10; /* half way up/down */
				if (i < nc)
					draw_one_cell(g, offsetx, j, 20, 10, gcol[i]);
				if (table)
					sprintf(tmp_str, "%4g", table[i]);
				else
					sprintf(tmp_str, "%4g", min + i*(max-min)/nc);
				gdImageString(g, font, offsetx + 25, j - 8, 
						(unsigned char *) tmp_str, black);
			}

			/* max value: */
			if (table)
				sprintf(tmp_str, "%4g", table[n_table-1]);
			else
				sprintf(tmp_str, "%4g", max);

			/* title: */
			j = px * m->rows / 2 - (nc + 4) * 5;
			gdImageString(g, font, offsetx, j, (unsigned char *) title, black);
			break;
		case 2:
			for (i = 0; i < nc; i++) {
				j = px * m->cols / 2 + (i - nc / 2) * 10; 
				draw_one_cell(g, j, m->rows * px + 10, 10, 20, gcol[i]);
			}
			sprintf(tmp_str, "%6g", min);
			j = px * m->cols / 2 - (nc + 4) * 5;
			gdImageString(g, font, j, m->rows * px + 35, 
					(unsigned char *) tmp_str, black);
			sprintf(tmp_str, "%6g", max);
			j = px * m->cols / 2 + (nc - 4) * 5;
			gdImageString(g, font, j, m->rows * px + 35, 
					(unsigned char *) tmp_str, black);
			j = px * m->cols / 2 + (1 + nc / 2) * 10; /* right */
			gdImageString(g, font, j, m->rows * px + 20, 
					(unsigned char *) title, black);
			break;
		case 3:
			for (class = min; class <= max; class++) {
				i = class - min; /* [0...nc> */

				assert(i < nc);

				j = px * m->rows / 2 + (i - nc / 2) * 15; /* half way up/down */
				draw_one_cell(g, m->cols * px + 10, j, 20, 10, gcol[i]);
#ifdef HAVE_LIBCSF
				tmp_str[0] = '\0';
				for (k = 1; k < MgetNrLegendEntries(m->CSF_MAP); k++) {
					if (i == l[k].nr) {
						sprintf(tmp_str, "%s", l[k].descr);
						printf("%d:[%s]\n", k, tmp_str);
						break;
					}
				}
#endif
				if (tmp_str[0] == '\0')
					sprintf(tmp_str, "%3d", (int) floor(class));
				gdImageString(g, font, m->cols * px + 35, j, 
						(unsigned char *) tmp_str, black);
			}
			j = px * m->rows / 2 + (-1 - nc / 2) * 15; /* above */
			gdImageString(g, font, m->cols * px + 10, j, 
					(unsigned char *) title, black);
			break;
		default:
			printf("unknown legend option, %d\n", legend);
			exit(1);
	} /* switch legend */

	/* draw ticks */
	if (ticks) {
		if (yt) {
			assert(m->rows == yt->n);
			for (i = 0; i < m->rows; i += xt->every)  {
				gdImageLine(g, m->cols * px, (i+0.5) * px, 
					m->cols * px + 3, (i+0.5) * px, black);
				gdImageString(g, font, m->cols * px + 5, (i+0.5) * px - 7, 
					(unsigned char *) xt->entries[i], black);
			}
		} else {
			for (i = 0; i < m->rows; i += 30/px)  {
				gdImageLine(g, m->cols * px, (i+0.5) * px, 
					m->cols * px + 3, (i+0.5) * px, black);
				sprintf(tmp_str, "%3g", m->y_ul - (i+0.5) * m->cellsizey);
				gdImageString(g, font, m->cols * px + 5, (i+0.5) * px - 7, 
					(unsigned char *) tmp_str, black);
			}
		}
		if (xt) {
			assert(m->cols == xt->n);
			for (i = 0; i < yt->n; i += yt->every)  {
				gdImageLine(g, (i+0.5) * px, m->rows * px, 
					(i+0.5) * px, m->rows * px + 3, black);
				gdImageString(g, font, (i+0.5) * px - 10, m->rows * px + 2,
					(unsigned char *) yt->entries[i], black);
			}
		} else {
			for (i = 0; i < m->cols; i += 30/px)  {
				gdImageLine(g, (i+0.5) * px, m->rows * px, 
					(i+0.5) * px, m->rows * px + 3, black);
				sprintf(tmp_str, "%3g", m->x_ul + (i+0.5) * m->cellsizex);
				gdImageString(g, font, (i+0.5) * px - 10, m->rows * px + 2,
					(unsigned char *) tmp_str, black);
			}
		}
	}

	/* write output to f_name */
	if (strcmp(f_name, "-")) /* they differ */
		out = fopen(f_name, "wb");
	else
		out = stdout;
#ifdef HAVE_GDIMAGEPNG
	gdImagePng(g, out);
#else
	gdImageGif(g, out);
#endif
	if (strcmp(f_name, "-")) /* they differ */
		fclose(out);

	/* clean up: */
	gdImageDestroy(g);
	map_free(m);

	return 0;
}

static void draw_one_cell(gdImagePtr g, int x, int y, 
		int pixx, int pixy, int color) {
	int i, j;
	for (i = 0; i < pixx; i++)
		for (j = 0; j < pixy; j++)
			gdImageSetPixel(g, x+i, y+j, color);
}

static double find_reasonable_hash_interval(double range) {
  double d = 1.0;

  /* range = max - min; */
  if (range > 5.0) {
    while(1) {
      if (range / d < 6.0) return d;
      d *= 2.0;
      if (range / d < 6.0) return d;
      d *= 2.5;
      if (range / d < 6.0) return d;
      d *= 2.0;
    }
  } else {
    while(1) {
      if (range / d > 2.0) return d;
      d /= 2.0;
      if (range / d > 2.0) return d;
      d /= 2.5;
      if (range / d > 2.0) return d;
      d /= 2.0;
    }
  }
}

static double find_reasonable_hash_start(double min, double max, double hash_interval) {
  int i;
  double d = 0.0;
  
  if (max > 0.0 && min < 0.0) {
  	/* return 0.0; */
  	while (d > min)
  		d -= hash_interval;
  	return d;
  }
  i = ((int) (min / hash_interval));
  return ((double) i) * hash_interval;
}

static int read_table(char *s) {
	char *dup, *token = NULL;

	dup = string_dup(s); /* mess this one up */
	table = (float *) emalloc(TABLE_SIZE * sizeof(float));
	printf("table:\n");
	while ((token = strtok(n_table == 0 ? dup : NULL, ", ")) != NULL) {
		table[n_table] = atof(token);
		printf("%f ", table[n_table]);
		n_table++;
		if (n_table == TABLE_SIZE)
			ErrMsg(ER_IMPOSVAL, "table size too big");
	}
	printf("\n");
	return n_table;
}

static int find_entry(float *table, int n_table, float value) {
	int i = 0;

	if (value < table[0])
		return -1;
	while (value > table[i])
		i++;
	if (i == n_table)
		return -1;
	return i - 1;
}
#endif /* #ifdef HAVE_LIBGD */
