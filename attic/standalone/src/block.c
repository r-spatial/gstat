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
 * block.c: calculate block discretization at a prediction location
 */
#include <stdio.h>
#include <math.h>

#include "defs.h"
#include "userio.h"
#include "debug.h"
#include "data.h"
#include "utils.h"
#include "glvars.h"
#include "block.h"

/* Gaussian quadrature adapted from Carr & Palmer, MG 25(5), p. 514:
gauss[4]: -0.8611363116, -0.3399810436, 0.3399810436, 0.8611363116
gauss_w[4]: 0.3478548452, 0.6521451548, 0.6521451548, 0.3478548452
both divided by 2
*/

static double gauss[4] = { 
		-0.4305681558,
		-0.1699905218,
		0.1699905218,
		0.4305681558
	};
double gauss_w[4] = {
		0.1739274226,
		0.3260725774,
		0.3260725774,
		0.1739274226
	};

int reset_anyway = 0;

void reset_block_discr(void) {
	reset_anyway = 1;
}

DATA *block_discr(DATA *d, const DPOINT *block, const DPOINT *where) {
/*
 * Purpose       : get locations that discretize the estimation area
 * Created by    : Edzer J. Pebesma                    
 * Date          : april 6th, 1992 
 * Modified      : jan 13, 1993, june 3, 1993, june 23 1993
 * Prerequisites : where;
 *                 the contents of block (if not NULL) may not change
 *                 during the run of the program: only once the full
 *                 discretisation is calculated.
 * Returns       : DATA list with 1 point (x,y,z) or a list;
 * Side effects  : dynamic memory (re)allocation for d;
 * uses IS_BLOCK(where) to determine whether point or block is needed;
 * discretisations is rectangular using dimensions of *block;
 * 
 * only the x, y and z field are used for locations, 
 * weights are put in the weight union.
 */
	int i, k, l, m, ndim;
	int ndiscr = 0, restart = 0;
/*
 * if restart != 0, block descretizations are re-calculated.
 */

/* 
 * ndim: number of dimensions (x, y, z, ..?) in data 
 * ndiscr: total number of points discretizing the block in this run
 * lastndiscr: number of discretizing points in last call
 * gl_nblockdiscr: number of points in each dimension in discretizing grid
 */

	static double *dx = NULL, *dy = NULL, *dz = NULL;
	static double *weight = NULL;
	static int max = 0;
	DATA *area = NULL, **data;

/* 
 * d has the points of the discretization of the estimation area;
 * dx, dy, dz and weight are remembered in case blockdiscretization has not
 * changed 
 */
	if (IS_POINT(where) && d != NULL) { /* point shortcut */
		d->list[0]->x = where->x; 
		d->list[0]->y = where->y; 
		d->list[0]->z = where->z; 
		d->list[0]->u.weight = 1.0; /* weight */
		d->n_list = 1;
		if (DEBUG_BLOCK) {
			printlog("block discretization (dist is weight):\n");
			d->mode =  X_BIT_SET | Y_BIT_SET | Z_BIT_SET;
			print_data_list(d);
		}
		return d;
	}

/* 
 * calculate the neccesary number of block discrimination points: 
 */
	ndim = 0; 
	ndiscr = 1;
	if (IS_BLOCK(where)) {
		if ((area = get_data_area()) == NULL) {
			if (block->x > 0.0) ndim++;
			if (block->y > 0.0) ndim++;
			if (block->z > 0.0) ndim++;
			for (i = 0; i < ndim; i++)
				ndiscr *= gl_nblockdiscr;
		} else
			ndiscr = area->n_list;
	}

	if (area == NULL && ndiscr > max) {
		dx = (double *) erealloc(dx, ndiscr * sizeof(double));
		dy = (double *) erealloc(dy, ndiscr * sizeof(double));
		dz = (double *) erealloc(dz, ndiscr * sizeof(double));
		weight = (double *) erealloc(weight, ndiscr * sizeof(double));
		max = ndiscr;
	}
/*
 * (re)allocate the memory neccesary: 
 */
	if (d == NULL) { /* allocate ndiscr: */
		d = (DATA *) emalloc (sizeof(DATA));
		init_one_data(d);
		assert(get_n_vars() > 0);
		data = get_gstat_data();
		d->pp_norm2 = data[0]->pp_norm2; /* assign distance function */
		d->what_is_u = U_ISWEIGHT;
		d->n_X = 0;
		d->colnx = d->colny = d->colnz = 1;
		d->list = (DPOINT **) emalloc(ndiscr * sizeof(DPOINT *));
		for (i = 0; i < ndiscr; i++) 
			d->list[i] = (DPOINT *) emalloc (sizeof(DPOINT));
		d->n_max = d->n_list = ndiscr;
		restart = 1;
	} else if (ndiscr > d->n_max) { /* resize if n_max < ndiscr: */
		d->list = (DPOINT **) erealloc(d->list, ndiscr * sizeof(DPOINT *));
		for (i = d->n_max; i < ndiscr; i++) 
			d->list[i] = (DPOINT *) emalloc (sizeof(DPOINT));
		d->n_max = d->n_list = ndiscr;
		restart = 1;
	} else if (reset_anyway) {
		reset_anyway = 0;
		restart = 1;
	}

	if (restart && ndim > 0 && area == NULL) { 
		/* set up block regular or Gaussian block discretisation */
		i = 0;
		switch (ndim) {
			case 1:
				if (block->y > 0 || block->z > 0)
					ErrMsg(ER_IMPOSVAL, 
						"block_discr(): block y and z dimensions must be 0");
				if (gl_gauss) {
					for (k = 0; k < 4; k++) {
						dx[i] = block->x * gauss[k];
						dy[i] = dz[i] = 0.0;
						weight[i] = gauss_w[k];
						i++;
					}
				} else {
					for (k = 0; k < gl_nblockdiscr; k++) {
						dx[i] = block->x * 
							(-0.5 + (1.0/gl_nblockdiscr) * (0.5 + k));
						dy[i] = dz[i] = 0.0;
						weight[i] = 1.0 / (1.0 * ndiscr);
						i++;
					}
				}
				break;
			case 2: 	
				if (block->z > 0)
					ErrMsg(ER_IMPOSVAL, "block_discr(): 2D block->z must be zero");
				if (gl_gauss) {
					for (k = 0; k < 4; k++) {
						for (l = 0; l < 4; l++) {
							dx[i] = block->x * gauss[k];
							dy[i] = block->y * gauss[l];
							weight[i] = gauss_w[k] * gauss_w[l];
							dz[i] = 0.0;
							i++;
						} /* for l */
					} /* for k */
				} else {
					for (k = 0; k < gl_nblockdiscr; k++) {
						for (l = 0; l < gl_nblockdiscr; l++) {
							dx[i] = block->x * 
								(-0.5 + (1.0/gl_nblockdiscr) * (0.5 + k));
							dy[i] = block->y * 
								(-0.5 + (1.0/gl_nblockdiscr) * (0.5 + l));
							dz[i] = 0.0;
							weight[i] = 1.0 / (1.0 * ndiscr);
							i++;
						} /* for l */
					} /* for k */
				}
				break;
			case 3:
				if (gl_gauss) {
					for (k = 0; k < 4; k++) {
						for (l = 0; l < 4; l++) {
							for (m = 0; m < 4; m++) {
								dx[i] = block->x * gauss[k];
								dy[i] = block->y * gauss[l];
								dz[i] = block->z * gauss[m];
								weight[i] = gauss_w[k] *
									gauss_w[l] * gauss_w[m];
								i++;
							} /* for m */
						} /* for l */
					} /* for k */
				} else {
					for (k = 0; k < gl_nblockdiscr; k++) {
						for (l = 0; l < gl_nblockdiscr; l++) {
							for (m = 0; m < gl_nblockdiscr; m++) {
								dx[i] = block->x * 
									(-0.5+(1.0/gl_nblockdiscr) * (0.5 + k));
								dy[i] = block->y * 
									(-0.5+(1.0/gl_nblockdiscr) * (0.5 + l));
								dz[i] = block->z * 
									(-0.5+(1.0/gl_nblockdiscr) * (0.5 + m));
								weight[i] = 1.0 / (1.0 * ndiscr);
								i++;
							} /* for m */
						} /* for l */
					} /* for k */
				}
				break;
			default:
				ErrMsg(ER_IMPOSVAL, "block_discr()");
		} /* switch */
		if (i != ndiscr)  /* silly check ! */
			ErrMsg(ER_IMPOSVAL, "block_discr()");
	} /* end if restart */

	if (IS_BLOCK(where)) {
		if (area != NULL) {
			for (i = 0; i < ndiscr; i++) {
				d->list[i]->x = area->list[i]->x + where->x; 
				d->list[i]->y = area->list[i]->y + where->y; 
				d->list[i]->z = area->list[i]->z + where->z; 
				if (area->colnvariance)
					d->list[i]->u.weight = area->list[i]->variance; 
				else
					d->list[i]->u.weight = 1.0 / ndiscr; 
			}
		} else {
			for (i = 0; i < ndiscr; i++) {
				d->list[i]->x = where->x - dx[i]; 
				d->list[i]->y = where->y - dy[i]; 
				d->list[i]->z = where->z - dz[i]; 
				d->list[i]->u.weight = weight[i]; 
			}
		}
	} else { /* where is point */
		d->list[0]->x = where->x; 
		d->list[0]->y = where->y; 
		d->list[0]->z = where->z; 
		d->list[0]->u.weight = 1.0; /* weight */
	}
	d->n_list = ndiscr;
	if (DEBUG_BLOCK) {
		printlog("block discretization (dist is weight):\n");
		d->mode =  X_BIT_SET | Y_BIT_SET | Z_BIT_SET;
		print_data_list(d);
	}
	return d;
}
