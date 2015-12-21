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
 * direction.c: evaluation of directions for variogram pair inclusion
 * Optimized by Konstantin Malakhanov (K.M.), march 1998
 */
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "defs.h"
#include "userio.h"
#include "data.h"
#include "utils.h"
#include "debug.h"
#include "defaults.h"
#include "glvars.h"
#include "direct.h"

#define MAX_ANG (PI/2.0)

/*
 * direct's little private data base:
 */
static double alpha = 0.0, beta = 0.0,
	tol_hor = PI, tol_ver = PI,
	cos_tol_hor = -1.0, cos_tol_ver = -1.0, /* Changed K.M. Fri Feb 27 11:13:47 1998 */
	dir_v[2] = { 1.0, 0.0 } , dir_h[2] = { 0.0, 1.0 };
static int all_directions = -1;

void set_direction_values(double a, double b, double t_h, double t_v) {

	/* do some checks: */
	if (a < 0.0 || a > 360.0)
		pr_warning( "alpha must be in [0..360]");
	if (b < 0.0 || b > 360.0)
		pr_warning( "beta must be in [0..360]");
	if (t_h < 0.0 || t_h > DEF_tol_hor)
		pr_warning( "horizontal tolerance must be in <0..180>");
	if (t_v < 0.0 || t_v > DEF_tol_ver)
		pr_warning( "vertical tolerance must be in <0..180>");
	if (t_h == DEF_tol_hor && t_v == DEF_tol_ver) {
		all_directions = 1;
		return;
	} else
		all_directions = 0;
	alpha = a * PI / 180.0;
	beta = b * PI / 180.0;
	tol_hor = t_h * PI / 180.0;
	tol_ver = t_v * PI / 180.0;
	cos_tol_hor = cos(tol_hor); /* Changed K.M. Fri Feb 27 11:14:15 1998 */
	cos_tol_ver = cos(tol_ver); /* Changed K.M. Fri Feb 27 11:14:15 1998 */
	dir_h[0] = sin(alpha); /* <x> */
	dir_h[1] = cos(alpha); /* <y> */
	dir_v[0] = cos(beta); /* <x,y> */
	dir_v[1] = sin(beta); /* <z> */
	return;
}

double valid_direction(DPOINT *a, DPOINT *b, int symmetric, const DATA *d) {
/*
 * return distance when vector is within the tolerance section;
 * return -1.0 when vector is outside tolerance section.
 */
	double norm, inprod, dist, px, py, pz;
	/* Changed K.M. Fri Feb 27 11:06:07 1998 */
	
	/* dist = d->point_norm(p); */
	dist = sqrt(d->pp_norm2(a, b));
	
	if (all_directions == 1)
		return dist;
	
	px = a->x - b->x;
	py = a->y - b->y;
	pz = a->z - b->z;

	if (tol_hor >= MAX_ANG && tol_ver >= MAX_ANG)
		return dist;
	if (tol_hor >= MAX_ANG && pz == 0.0)
		return dist;
	if (tol_ver >= MAX_ANG && px == 0.0 && py == 0.0)
		return dist;
	if (tol_hor < MAX_ANG && (px != 0.0 || py != 0.0)) {
		/* 
		 * check in <x,y> plane:
		 */
		norm = sqrt(px * px + py * py);
		inprod = (px * dir_h[0] + py * dir_h[1])/norm;

		if (symmetric) { /* the most often case */
			if ( fabs(inprod) < cos_tol_hor) /* if cos(alpha) < cos(tol) then alpha > tol! */
				return -1.0;
		} else if (inprod < cos_tol_hor)
			return -1.0;

		/* Changed K.M. Fri Feb 27 11:18:24 1998: I have
		eliminated most of checks for the following reasons: if
		(because of numerical reasons) we get fabs(inprod)>1.0,
		then it really means that inprod is either +1.0 or
		-1.0, i.e. the angle between point vector and variogram
		direction is either 0 or pi. In case of SYMMETRIC
		it doesn't matter - fabs(inprod) is then bigger
		then cos_tol_hor, so we accept this point pair as it
		should be. In the case of not SYMMETRIC - if inprod>1
		then inprod is bigger then cos_tol_hor and we accept
		this pair, if inprod < -1.0 then it is smaller then
		cos_tol_hor and we reject this pair. This way we also
		do not need to calculate acos and to check for domain
		errors afterwards. */

	}
	if (tol_ver < MAX_ANG && (px != 0.0 || py != 0.0 || pz != 0.0)) { 
		/* 
		 * inproduct in <xy, z> 
		 */
		/* Changed K.M. Fri Feb 27 15:47:13 1998 */
		inprod = (sqrt(px * px + py * py) * dir_v[0] + pz * dir_v[1])/dist;
		if (symmetric) { /* the most often case */
			if ( fabs(inprod) < cos_tol_ver) /* if cos(alpha) < cos(tol) then alpha > tol! */
			  return -1.0;
		} else if (inprod < cos_tol_ver)
			return -1.0;
	}
	return dist;
}
#undef MAX_ANG
