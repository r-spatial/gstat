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
 * vario_io.c: functions for point-point, point-block (co/semi)variances
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "userio.h"
#include "debug.h"
#include "data.h"
#include "glvars.h" /* get_n_outfile() */
#include "vario.h"
#include "block.h"
#include "vario_io.h"

static double sem_cov_blocks(VARIOGRAM *v, DATA *a, DATA *b, int sem);

double sem_cov_ab(VARIOGRAM *v, DPOINT *a, DPOINT *b, int sem)
/*
 * return Cov(a,b) or Sem(a,b),
 * taking care of IS_BLOCK(a) and IS_BLOCK(b):
 */
{
	static DATA *Discr_a = NULL, *Discr_b = NULL;
	static DPOINT *block_p = NULL;
	DPOINT *tmp;

	if (block_p == NULL)
		block_p = get_block_p();
	if (a == b) {
		if (IS_POINT(a))
			return sem_cov_blocks(v, NULL, NULL, sem);
		Discr_a = block_discr(Discr_a, block_p, a);
		return sem_cov_blocks(v, Discr_a, Discr_a, sem);
	}
	/*
	 * if one of them IS_BLOCK, make sure it's a: 
	 * (because block_discr() will otherwise store block
	 * discretisations in both Discr_a and Discr_b)
	 */
	if (IS_POINT(a) && IS_BLOCK(b)) {  
		tmp = a; a = b; b = tmp; /* swap a and b */
	}
	Discr_a = block_discr(Discr_a, block_p, a);
	Discr_b = block_discr(Discr_b, block_p, b);
	return sem_cov_blocks(v, Discr_a, Discr_b, sem);
}

static double sem_cov_blocks(VARIOGRAM *v, DATA *a, DATA *b, int sem) {
/*
 * Purpose       : calculates (once) and returns Cov(a,b);
 *                 if a==b && a denotes a block discretisation, the value
 *                 is put in a->block_variance and a->block_xxx_set gets 1
 * Created by    : Edzer J. Pebesma
 * Date          : 25 jan 1992, 3 june 1993
 * Prerequisites : none 
 * Returns       : sem == 1 ? Sem(a,b) : Cov(a,b) 
 * Side effects  : none 
 */
	int i, j;
	double block_value, dx, dy = 0.0, dz = 0.0, dist, ret, weight, dzero2;
	DPOINT *dpa, *dpb;

	/*
	if (a->what_is_u != U_ISWEIGHT || b->what_is_u != U_ISWEIGHT)
		ErrMsg(ER_IMPOSVAL, "weights needed in SevCov_Blocks()");
	*/
	if (a == NULL)
		return sem ? get_semivariance(v, 0.0, 0.0, 0.0) : 
				get_covariance(v, 0.0, 0.0, 0.0);
	if (a->n_list == 1 && b->n_list == 1) { /* point--point */
		if (gl_longlat) {
			if (! v->isotropic)
				ErrMsg(ER_IMPOSVAL, "for long/lat data, anisotropy cannot be defined");
			dist = pp_norm_gc(a->list[0], b->list[0]);
			ret = sem ?  get_semivariance(v, dist, 0.0, 0.0):
				get_covariance(v, dist, 0.0, 0.0);
			/* printf("ll dist: %g, ret.val %g\n", dist, ret); */
			return ret;
		} else {
			return sem ?
				get_semivariance(v, a->list[0]->x - b->list[0]->x,
					a->list[0]->y - b->list[0]->y, a->list[0]->z - b->list[0]->z):
				get_covariance(v, a->list[0]->x - b->list[0]->x,
					a->list[0]->y - b->list[0]->y, a->list[0]->z - b->list[0]->z);
		}
	}
	/* now a->n_list > 1 or b->n_list > 1: block--block or point--block */
	if (gl_longlat)
			ErrMsg(ER_IMPOSVAL, "block kriging for long-lat data undefined");
	if (a == b) { /* block--block for a single block */
		if (sem && v->block_semivariance_set)
	    	return v->block_semivariance;
		if (!sem && v->block_covariance_set)
	    	return v->block_covariance;
	}
	/* else: continue */
	dzero2 = gl_zero * gl_zero;
	block_value = 0.0;
	for (i = 0; i < a->n_list; i++) { 
		for (j = 0; j < b->n_list; j++) { /* compare all points with all p */
			dpa = a->list[i];
			dpb = b->list[j];
			weight = dpa->u.weight * dpb->u.weight;
			/* avoid the "zero-effect" if a or b is block: */
			dx = dpa->x - dpb->x;
			dy = dpa->y - dpb->y;
			dz = dpa->z - dpb->z;
			/* avoid ``zero-effect'': */
			if (a->pp_norm2(dpa, dpb) < dzero2) {
				dx = (dx >= 0 ? gl_zero : -gl_zero);
				if (a->mode & Y_BIT_SET)
					dy = (dy >= 0 ? gl_zero : -gl_zero);
				if (a->mode & Z_BIT_SET)
					dz = (dz >= 0 ? gl_zero : -gl_zero);
			}
			if (sem) 
				block_value += weight * get_semivariance(v, dx, dy, dz);
			else
				block_value += weight * get_covariance(v, dx, dy, dz);
		} /* for j */
	} /* for i */
	if (a == b) { /* remember within block cov./sem.: */
		if (sem) {
			v->block_semivariance = block_value;
			v->block_semivariance_set = 1;
		} else {
			v->block_covariance = block_value;
			v->block_covariance_set = 1;
		}
	}
	return block_value;
}
