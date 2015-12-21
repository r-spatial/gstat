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
    (at your option) any later version.

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

#include "random.h"
#include "utils.h"
#include "userio.h"
#include "mapio.h"

#define MAX_SAMPLE 1000

typedef struct {
	int r, c;
	float x, y;  /* x and y coordinate in [0,1]x[0,1] area */
} PTS;

typedef struct {
	int max_size, size;
	float xy_ratio; /* dx/dy */
	PTS *pts;
} PTS_LIST;

static PTS_LIST *add_to_pts_list(int row_offset, int col_offset, 
		PTS_LIST *sample, PTS_LIST *add_to);
static PTS_LIST *get_sample(PTS_LIST *idx, int size, int nrows, int ncols);
static PTS_LIST *resize_pts_list(PTS_LIST *p, int new_max_size);

#define SAMPLE_USAGE "\
usage: gstat -e sample [options] input_map output_file\n\
options:\n\
 -t sampling_type (1 = simple random, 2 = stratified random, 3 = \n\
      systematic aligned, 4 = stratified systematic unaligned)\n\
 -s sample_size (within blocks; default = 1)\n\
 -k block_size (in number of cells; blocks are allways square)\n\
 -c nr_of_clusters (default: all blocks)\n\
 -r nr_of_replicates (default: 1)\n\
 -m write to maps instead of printing means"

int sample_main(int argc, char *argv[]) {

	GRIDMAP *m = NULL, *outmap = NULL;
	FILE *outfile = NULL;
	int i, j, n, r, row, col, c, size = 1, replicates = 1,
		write_maps = 0, k = -1, row_blocks, col_blocks, nclus = -1;
	float value;
	double mean;
	enum { SRS = 1, /* simple random sampling */
		STRS, /* stratified random sampling */
		SSA, /* systematic sampling aligned */
		SSU  /* stratified systematic unaligned */
	} type = SRS;
	PTS_LIST *index = NULL, *sample = NULL, *clusters = NULL;
	const char *ofilename = NULL;
	char filename[128];

	set_seed(0); /* time seed */
	index = resize_pts_list(NULL, 1);
	sample = resize_pts_list(NULL, 1);
	clusters = resize_pts_list(NULL, 1);

	opterr = 0;
	while ((c = getopt(argc, argv, "mt:s:k:r:c:")) != EOF) {
		switch (c) {
			case 'm': write_maps = 1; break;
			case 't': type = atoi(optarg); break;
			case 's': size = atoi(optarg); break;
			case 'k': k = atoi(optarg); break;
			case 'r': replicates = atoi(optarg); break;
			case 'c': nclus = atoi(optarg); break;
			case '?': default: ErrClo(optopt);
		}
	}

	if (optind == argc - 2) {
		m = new_map(READ_ONLY);
		m->filename = argv[optind];
		ofilename = argv[optind+1];
	} else /* if (optind != argc)  */ {
		printf("%s\n", SAMPLE_USAGE);
		exit(1);
	}

	if (type < 1 || type > 4)
		ErrMsg(ER_IMPOSVAL, "type values are 1, 2, 3 or 4\n");

	if (m->filename == NULL)
		ErrMsg(ER_ARGOPT, "missing map name\n");

	if (map_read(m) == NULL)
		ErrMsg(ER_READ, m->filename);

	if (type == 1) {
		row_blocks = col_blocks = 1;
		if (k > 0)
			ErrMsg(ER_IMPOSVAL, "cannot use positive k to SRS");
		k = MAX(m->rows, m->cols);
	} else {
		if (k <= 1)
			ErrMsg(ER_IMPOSVAL, "k should be larger than 1");
		row_blocks = m->rows / k + ((m->rows % k) ? 1 : 0);
		col_blocks = m->cols / k + ((m->cols % k) ? 1 : 0);
	}

	printf("sample parameters:\n===================\n");
	printf(    "sample size:%7d (within each block)\n", size);
	if (type != 1)
		printf("block size (k):%4d (%d row blocks, %d column blocks)\n", 
				k, row_blocks, col_blocks);
	if (replicates > 1) 
		printf("replicates: %7d\n", replicates);
	if (nclus > 0)
		printf("clusters:   %7d\n", nclus);

	if (size > k * k) {
		pr_warning("size value (%d) > k * k (%d)", size, k * k);
		ErrMsg(ER_IMPOSVAL, "size");
	}

	if ((m->rows % k) != 0 || (m->cols % k) != 0) {
		pr_warning("rows (%d) and/or columns (%d) are not a multiple of k (%d)",
			m->rows, m->cols, k);
		printf("Note that the final sample size may differ from the intended size\n");
		/* but proceed */
	}

	for (r = 0; r < replicates; r++) {
		/* draw sample locations */
		index->size = 0;
		switch (type) {
			case SRS:
				if (nclus > 0)
					ErrMsg(ER_IMPOSVAL, "cannot apply clustered SRS");
				index = add_to_pts_list(0, 0, get_sample(index, size, k, k), index);
				break;
			case STRS:
				if (nclus <= 0)
					nclus = row_blocks * col_blocks; /* all */
				clusters = get_sample(clusters, nclus, row_blocks, col_blocks);
				for (c = 0; c < nclus; c++) {
					/* printf("cluster %d: r %d c %d\n", c, clusters.pts[c].r, clusters.pts[c].c); */
					sample = get_sample(sample, size, k, k);
					index = add_to_pts_list(k * clusters->pts[c].r, 
							k * clusters->pts[c].c, sample, index);
				}
				break;
			case SSA:
				if (nclus <= 0)
					nclus = row_blocks * col_blocks; /* all */
				sample = get_sample(sample, size, k, k);
				clusters = get_sample(clusters, nclus, row_blocks, col_blocks);
				for (c = 0; c < nclus; c++)
					index = add_to_pts_list(k * clusters->pts[c].r, 
							k * clusters->pts[c].c, sample, index);
				break;
			case SSU:
				if (nclus > 0 && r == 0)
					ErrMsg(ER_IMPOSVAL, "-c: clustered SSU is not allowed");
				
				sample = get_sample(sample, 1, 1, 1); /* just for allocation */

				for (c = 0; c < size; c++) {
					clusters = get_sample(clusters, MAX(row_blocks, col_blocks),
							k, k); /* locations */
					for (i = 0; i < row_blocks; i++) {
						/* for row i, fix col: */
						sample->pts[0].c = clusters->pts[i].c;
						for (j = 0; j < col_blocks; j++) {
							sample->pts[0].r = clusters->pts[j].r;
							index = add_to_pts_list(i * k, j * k,
									sample, index);
						}
					}
				}
				break;
		}

		/* output */
		if (write_maps) { /* open map */
			if (replicates > 1) {
				map_name_nr(m, ofilename, filename, r, replicates);
				outmap = map_dup(filename, m);
				printf("writing map %s\n", filename);
			} else
				outmap = map_dup(ofilename, m);
		} else if (r == 0) /* or open file */
			outfile = efopen(ofilename, "w");

		/* sample from map */
		for (i = 0, mean = 0.0, n = 0; i < index->size; i++) {
			row = index->pts[i].r;
			col = index->pts[i].c;
			if (row < m->rows && col < m->cols && !map_cell_is_mv(m, row, col)) {
				value = map_get_cell(m, row, col);
				/* printf("%d %d %g\n", row, col, value); */
				if (write_maps) 
					map_put_cell(outmap, row, col, value);
				else {
					mean += value;
					n += 1;
				}
			}
		}

		/* process sample */
		if (write_maps) {
			outmap->write(outmap);
			map_free(outmap);
		} else
			fprintf(outfile, "%12.6f\n", mean / n);

	} /* for replicates */
	if (! write_maps)
		efclose(outfile);
	return 0;
}

static PTS_LIST *get_sample(PTS_LIST *idx, int size, int nrows, int ncols) {
	/* static PTS_LIST idx = { 0, 0, NULL }; */
	static int old_blocksize = 0;
	static unsigned char **field = NULL;
	int i, j, row, col, block_size;

	assert(size <= nrows * ncols);
	assert(size > 0);

	if (size > idx->max_size) 
		resize_pts_list(idx, size);

	block_size = MAX(nrows, ncols);

	if (old_blocksize < block_size) {
		field = (unsigned char **) emalloc(block_size * sizeof(unsigned char *));
		for (i = 0; i < block_size; i++)
			field[i] = (unsigned char *) emalloc(block_size * sizeof(unsigned char));
		old_blocksize = block_size;
	}

	for (i = 0; i < nrows; i++) {
		for (j = 0; j < ncols; j++) {
			field[i][j] = 0;
		}
	}

	i = 0;
	while (i < size) {
		row = nrows * r_uniform();
		col = ncols * r_uniform();
		assert(row >= 0 && row < nrows);
		assert(col >= 0 && col < ncols);
		if (field[row][col] == 0) {
			field[row][col] = 1;
			idx->pts[i].r = row;
			idx->pts[i].c = col;
			i++;
		}
	}
	idx->size = size;
	return idx;
}

static PTS_LIST *add_to_pts_list(int row_offset, int col_offset, 
		PTS_LIST *sample, PTS_LIST *add_to) {

	int i, newsize;

	newsize = add_to->size + sample->size;
	resize_pts_list(add_to, newsize);

	for (i = 0; i < sample->size; i++) {
		add_to->pts[add_to->size + i].r = sample->pts[i].r + row_offset;
		add_to->pts[add_to->size + i].c = sample->pts[i].c + col_offset;
	}

	add_to->size = newsize;
	add_to->max_size = MAX(add_to->max_size, newsize);
	return add_to;
}

static PTS_LIST *resize_pts_list(PTS_LIST *p, int new_max_size) {
	if (new_max_size <= 0)
		ErrMsg(ER_IMPOSVAL, "size should be > 0");
	if (p == NULL) {
		p = (PTS_LIST *) emalloc(sizeof(PTS_LIST));
		p->pts = (PTS *) emalloc(new_max_size * sizeof(PTS));
	} else
		p->pts = (PTS *) erealloc(p->pts, new_max_size * sizeof(PTS));
	p->max_size = new_max_size;
	return p;
}
