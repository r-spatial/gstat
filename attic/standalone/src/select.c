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
 * select.c: neighborhood selection for local prediction
 */
#include <stdio.h>
#include <stdlib.h> /* qsort() */
#include <string.h> /* memcpy() */
#include <math.h> /* sqrt(), fabs() */

#include "defs.h"

#ifdef USING_R
void Rprintf(const char *,...);
#endif

#include "userio.h"
#include "data.h"
#include "utils.h"
#include "debug.h"
#include "vario.h"
#include "defaults.h"
#include "glvars.h"
#include "nsearch.h"
#include "select.h"
#include "polygon.h"
#ifdef HAVE_EXT_DBASE
# include "ext_dbase.h"
#endif

static int octant_select(DATA *d, DPOINT *where);
static int which_octant(DPOINT *where, DPOINT *p, int mode);

int CDECL dist_cmp(const DPOINT **ap, const DPOINT **bp);
static void zero_sel_dist2(DATA *d);
static void print_selection(DATA *d, DPOINT *where);

#define store_radius(d)    keep_radius(d, 0)
#define set_back_radius(d) keep_radius(d, 1)
/* beware-of-side-effect macro: don't call with ++/--'s */
#define DPSWAP(a,b) { if (a != b) { tmp = a; a = b; b = tmp; }}

static int select_qtree(DATA *d, DPOINT *where)
{
	if (d->n_list <= 0 || d->id < 0 || d->sel_max == 0)
		return (d->n_sel = 0);
/* 
 * global neighbourhood?
 */
	if (IS_GLOBAL(d) || where == NULL) {
		d->sel = d->list;
		d->n_sel = d->n_sel_max = d->n_list;
        if (get_n_edges()) {
            /* memcpy(d->sel, d->list, d->n_list * sizeof(DPOINT *)); */
            zero_sel_dist2(d);
            if (where != NULL)
				check_edges(d, where);
        }
        if (DEBUG_SEL) 
        	print_selection(d, where);
        
		return d->n_sel;
	}
	/* set up or resize selection list, d->sel */
	if (d->sel == NULL || d->sel == d->list) {
 		/* first time or after forced global selection: */
		d->sel = (DPOINT **) emalloc(d->n_list * sizeof(DPOINT *));
		d->n_sel_max = d->n_list;
	} else if (d->n_list > d->n_sel_max) { /* no, something has changed... */
		d->n_sel_max += MAX(MAX_DATA, d->n_list - d->n_sel_max);
		d->sel = (DPOINT **) erealloc(d->sel, d->n_sel_max * sizeof(DPOINT *));
	}
 	if (d->id > 0) { /* 2nd, 3d,..., n_vars-th variable */
		if (gl_coincide == DEF_coincide)
			gl_coincide = decide_on_coincide(); /* establish first ... */
 		if (gl_coincide) {
	        int i;
	        DATA **data = get_gstat_data();
 			d->n_sel = data[0]->n_sel;
			for (i = 0; i < d->n_sel; i++) /* copy previous selection: */
				d->sel[i] = d->list[GET_INDEX(data[0]->sel[i])]; 
            if (DEBUG_SEL) 
            	print_selection(d, where);
 			return d->n_sel; /* and we're done!! */
 		}
 	} 

/*
 * So far, no selection has been done.
 * Let's see if today's job is easily done:
 */
	memcpy(d->sel, d->list, d->n_list * sizeof(DPOINT *));
	if (d->sel_rad >= DBL_MAX && d->sel_max >= d->n_list && d->oct_max == 0) {
		if (get_n_edges()) {
			zero_sel_dist2(d);
			check_edges(d, where);
		} 
		d->n_sel = d->n_list;
		if (DEBUG_SEL) 
			print_selection(d, where);
		return d->n_sel;
	}

/*
 * now we're sure that (a) sel_rad is set, or (b) sel_max < n_list,
 * or (c) oct_max is set, so let's do the smart thing:
 */
	qtree_select(where, d);
	return -1; /* more work to do */
}

int select_at(DATA *d, DPOINT *where) {
/*
 * fill the array d->sel appropriatly given x,y,z,d->sel_min,d->sel_max
 * and d->sel_rad: possibly by corresponding semivariance value
 * first select all points within a distance d->sel_rad from where
 * then select at max the d->sel_max nearest and return no selection
 * if there are less then d->sel_min
 * if "FORCE", then select ALWAYS the d->sel_min nearest points.
 * 
 * corrected variogram distance feb. 16th 1993
 * changed search to normalizing to (0,0,0) first, aug 1993
 */

	if (get_method() == POLY) 
		return d->n_sel;

	if (d->what_is_u == U_UNKNOWN)
		d->what_is_u = U_ISDIST; /* we're going to fill this right now */
	else if (d->what_is_u != U_ISDIST)
		ErrMsg(ER_IMPOSVAL, "select_at() needs distances");

	/* CW */
	if (d->type.type==DATA_EXT_DBASE) {
#ifdef HAVE_EXT_DBASE
	  if (select_ext_dbase(d,where)!=-1) {
        if (DEBUG_SEL) 
           print_selection(d, where);
		return d->n_sel;
      }
#endif
	} else {
	  if (select_qtree(d,where)!=-1)
		return d->n_sel;
	}

	if (get_n_edges()) 
		check_edges(d, where);

/*
 * so, now we're at the stage where one of the following conditions holds:
 * (a) we selected all points within d->sel_rad
 * (b) we selected (at least) the nearest d->sel_max
 * (c) we selected (forced) at least d->sel_min, possibly beyond d->sel_rad
 * Now, should we select further, sorting on distance?
 *
 */

	if (d->vdist) { /* use variogram distance as sort criterium */
		int i;
		for (i = 0; i < d->n_sel; i++)
			d->sel[i]->u.dist2 = get_semivariance(get_vgm(LTI(d->id,d->id)),
				where->x - d->sel[i]->x,
				where->y - d->sel[i]->y,
				where->z - d->sel[i]->z);
	}

	if (d->oct_max) { /* do octant selection */
		d->oct_filled = octant_select(d, where);
		/* sorts, adjusts n_sel and destroys distance order, so only */
		if (get_method() == SPREAD) /* then we've got to re-order them */
			qsort(d->sel, (size_t) d->n_sel, sizeof(DPOINT *), (int
				CDECL (*)(const void *,const void *)) dist_cmp);
	} 

	if (d->vdist) {
		qsort(d->sel, (size_t) d->n_sel, sizeof(DPOINT *), (int 
			CDECL (*)(const void *, const void *)) dist_cmp);
		/* pick d->sel_[max|min] nearest: */
		if (d->sel_min && d->n_sel == d->sel_min
				&& d->sel[d->n_sel]->u.dist2 > d->sel_rad) /* we forced: */
			d->n_sel = d->sel_min;
		if (d->n_sel > d->sel_max)	
			d->n_sel = d->sel_max;
	} 

	if (DEBUG_SEL) 
		print_selection(d, where);
	return d->n_sel;
}

static int octant_select(DATA *d, DPOINT *where) {
/*
 * selects the d->oct_max nearest points per octant/quadrant/secant,
 * using euclidian or variogram distance
 */
	int i, j, noct = 1, n, start = 0, end = 0, n_notempty = 0;
	DPOINT **sel, *tmp;

	if (d->mode & Z_BIT_SET)
		noct = 8;
	else if (d->mode & Y_BIT_SET)
		noct = 4;
	else
		noct = 2;

	sel = d->sel;
	for (i = 0; i < noct; i++) { /* for each octant: */
		/* start == end: */
		for (j = start; j < d->n_sel; j++) { /* for the remaining of sel: */
			if (which_octant(where, sel[j], d->mode) == i) {
				DPSWAP(sel[end], sel[j]);
				end++; /* j >= end */
			}
		}
		n = end - start; /* # of pts in octant i */
		if (n > 0)
			n_notempty++; /* Yahoo, another non-empty octant! */
		if (n > d->oct_max) {
			/* to get the closest n: sort sel from start to end: */
			qsort(sel + start, (size_t) n, sizeof(DPOINT *), 
				(int CDECL (*)(const void *, const void *)) dist_cmp);
			/* swap the remaining ones to the end of sel and forget about'm: */
			for (j = start + d->oct_max; j < end; j++) {
				d->n_sel--;
				DPSWAP(sel[j], sel[d->n_sel]);
			}
			/* accept the first d->oct_max: */
			start += d->oct_max;
			/* proceed with the next octant: */
			end = start;
		} else /* accept all n: */
			start = end;
	}
	if (end != d->n_sel) {
		Rprintf("end: %d, n_sel: %d\n", end, d->n_sel);
		ErrMsg(ER_IMPOSVAL, "octant_select(): remaining points");
	}
	return n_notempty; /* # non-empty octants */
}

static int which_octant(DPOINT *where, DPOINT *p, int mode) {
/*
 * it's pretty hard to get this right for perfectly alligned 3D data:
 * in case of omax=1 and only one point at exactly equal distances
 * in each of the 6 directions, not all 6 are taken.
 * For 2D (omax=1, 4 points) this works, however.
 */
	double dx, dy, dz;
	int x = 0, y = 0, z = 0;

	dx = p->x - where->x; /* < 0 : p in west half */
	dy = p->y - where->y; /* < 0 : p in south half */
	dz = p->z - where->z; /* < 0 : p in lower half */
	if (mode & Z_BIT_SET)
		z = (dz < 0);
	if (mode & Y_BIT_SET) {
		x = (dy < 0 ? dx > 0 : dx >= 0);
		y = (dx < 0 ? dy >= 0 : dy > 0);
	} else 
		x = (where->x < p->x);
	return (x | (y << 1) | (z << 2));	
}

int CDECL dist_cmp(const DPOINT **pa, const DPOINT **pb) {
/* ANSI qsort() conformant dist_cmp */

	if ( (*pa)->u.dist2 < (*pb)->u.dist2 )
		return -1;
	if ( (*pa)->u.dist2 > (*pb)->u.dist2 )
		return 1;
	return 0;
} 

/* set u.idts2 of all selected points in DATA *d to 0 to be further
   proceeded by edges selection */
static void zero_sel_dist2(DATA *d) {
	int i;
    
	for (i = 0; i< d->n_sel; i++)
		d->sel[i]->u.dist2 = 0.0;
    
	return;        
}

static void print_selection(DATA *d, DPOINT *where) {

/* Add this statement to filter out
 * empty selections
   if (!d->n_sel)
     return;
 */

	if (where) {
		printlog("selection at "); 
		logprint_point(where, d);
	} else
		printlog("(NULL selection location)");
	print_data_selection(d);
}
