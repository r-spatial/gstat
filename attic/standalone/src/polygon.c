/*
    polygon.c -- as incorporated in
    Gstat, a program for geostatistical modelling, prediction and simulation
    1992, 2003, by Edzer J. Pebesma

    Authors: Edzer J. Pebesma (E.Pebesma@geo.uu.nl)
             Steve Joyce  (steve.joyce@resgeom.slu.se)
             Konstantin Malakhanov (malakhanov@iwwl.rwth-aachen.de)
             
    This code (polygon handling part) was initially written by Edzer J. Pebesma. 
    Edges handling part is originally by Steve Joyce.

    Simplified and merged by Konstantin Malakhanov.

    InPoly function is Written by Joseph O'Rourke, contributions by 
	Min Xu, June 1997.
    Questions to orourke@cs.smith.edu.
    --------------------------------------------------------------------
    InPoly is Copyright 1998 by Joseph O'Rourke.  It may be freely 
    redistributed in its entirety provided that this copyright notice is 
    not removed.
    --------------------------------------------------------------------

    Edzer J. Pebesma (e.pebesma@geo.uu.nl)
    Department of physical geography, Utrecht University
    P.O. Box 80.115, 3508 TC Utrecht, The Netherlands

    Steve Joyce                            mailto:steve.joyce@resgeom.slu.se
    Remote Sensing Laboratory                 http://www.resgeom.slu.se
    Dept. of Forest Resource Mgmt. & Geomatics          Tel: +46 90 16 69 57
    Swedish University of Agricultural Sciences         Fax: +46 90 14 19 15 
    S-901 83 Umeå, Sweden
    
    Konstantin Malakhanov
    Institut for Hydraulic Engineering & Groundwater Mgmt., 
    RWTH Aachen, Germany,

    --------------------------------------------------------------------
    
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

    (read also the files COPYING and Copyright) */

/*
 * polygon.c: point-in-polygon routine, if method=polygon;
 also edges handling for interpolation/simulation.
 Variogram modelling doesn't use this yet.
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "defs.h"
#include "data.h"
#include "utils.h"
#include "debug.h"
#include "userio.h"
#include "read.h"
#include "polygon.h"
#include "glvars.h"

/* masks  segment crossing codes - based on SPJ */

#define CROSS   0x01
#define SPOINT1 0x02
#define SPOINT2 0x04
#define EPOINT1 0x08
#define EPOINT2 0x10

/* point-to-segment, point-to-polygon relations */
#define POLY_IN  0x01
#define POLY_OUT 0x02
#define POLY_EDGE (POLY_IN|POLY_OUT)
#define IWW 1

/* for linear combination of coordinates:*/
#define DELTA_AB(_p,_a,_b,_r) (_p).x = (_a).x + (_r) * ((_b).x - (_a).x);\
(_p).y = (_a).y + (_r) * ((_b).y - (_a).y);

static char InPoly(PLOT_POINT q, POLYGON *Poly);
static int CDECL dist_cmp(const DPOINT **ap, const DPOINT **bp);
unsigned char line_of_sight(const PLOT_POINT data, const PLOT_POINT target,
                            const POLYGON *p_obj);
static unsigned int segment_cross_check(PLOT_POINT p1, PLOT_POINT p2, PLOT_POINT p3, PLOT_POINT p4, double *r, double *s);
static unsigned int segment_parallel_check(PLOT_POINT a, PLOT_POINT b, PLOT_POINT c, PLOT_POINT d);
static int segment_between_check(PLOT_POINT a, PLOT_POINT b, PLOT_POINT c);
static unsigned int check_open_edges(const PLOT_POINT data, const int target_side,
                                     const PLOT_POINT target, const POLYGON *p_obj);
static void print_poly_log(POLYGON *edge);

static const char *fname = NULL;
/*------------------------------ POLYGONS ------------------------------*/
#ifndef USING_R
POLYGON *read_polygons(const char *filename, int *n_polys, double **iso_values) {
	FILE *f = NULL;
	int i = 0, j,n, line_size = 0, d, np, iww_nl=-1;
	char *line = NULL, *cp;
	POLYGON *pol = NULL;
    double val, *values = NULL;
    
	fname = filename;
	f = efopen(filename, "r");
/*	pol = (POLYGON *) emalloc(sizeof(POLYGON));*/
	get_line(&line, &line_size, f);
	if (strstr(line, "EXP") == line) { /* read arc-info e00 file */
		if ((cp = strchr(line + 7, ' ')) != NULL)
			*cp = '\0';
		cp = line + 7;
		printlog("reading e00 polygon coverage `%s'", cp);
		get_line(&line, &line_size, f); /* ARC  # */
		n = 0;
		while (n >= 0) {
			get_line(&line, &line_size, f);  /* ARC  # */
            /* I have no examples where ARC comes before EVERY line..*/
			if (sscanf(line, "%d%d%d%d%d%d%d", &n, &d, &d, &d, &d, &d, &np) != 7) 
				pr_warning("cannot read 7 integers from `%s'", line);
			if (n > 0) {
				i++;
				pol = erealloc(pol, i * sizeof(POLYGON));
				pol[i - 1] = read_n_points(f, np);
                /* we do not delete polygons of 2 points, as we are
                   not neccessary reading closed polygons*/
                   
				if (np <= 1) { /* remove the degenerate polygon */
					efree(pol[i - 1].p);
					i--;
				}
			}
		}
	} else { /* read `plot' file */
		printlog("reading plot file `%s'", filename);
		do {
			d=sscanf(line, "%d %lf", &np,&val);
            if (d==0 || d==EOF) {
#ifdef IWW            
                if (i==0 && d==0) { /* try IWW poly/isolines */
                    printlog(" IWW poly/isolines? 7 header lines + number of isolines:\n");
/*                    printlog("%d) %s",1,line) ;*/
                    for (j=1;j<=7;j++)
                    {
                        printlog("%d) %s",j,line) ;
                        if (get_line(&line, &line_size, f) ==NULL) 
                            ErrMsg(ER_READ,filename);
                    }
                    if (sscanf(line, "%d", &iww_nl) !=1)
                        ErrMsg(ER_RDINT,line);
                    else
                        continue; /* get next line with np etc.*/
                    
                } else
#endif                
                {
                    pr_warning("file: %s, line without data encountered \n", fname);
                    continue;
                }
            }
            
            if (np<=1) {
                printf(" (polyline size %d encountered)",np);
                continue;
            }
            
			i++;
			pol = erealloc(pol, i * sizeof(POLYGON));

            if (d==2) {
                values = erealloc(values, i * sizeof(double));
                values[i - 1]=val;
            }
            
            
			pol[i - 1] = read_n_points(f, np);
			if (np <= 1) { /* remove the degenerate polygon */
				efree(pol[i - 1].p);
				i--;
			}
		} while (get_line(&line, &line_size, f));
	}
	printlog(", %d polygons read\n", i);
#ifdef    IWW
    if (iww_nl != -1 && iww_nl != i)
        pr_warning("the given number of isolines differs from number of actual loaded.");
#endif
	*n_polys = i;
    if (iso_values)
		*iso_values=values;
    
	efclose(f);
	efree(line);
	return pol;
}
#endif

void report_edges(void) {
    int i, j, n, *ip;
    char delim;
    POLYGON **edges;

    delim = (DEBUG_NORMAL || DEBUG_DUMP ? ' ': '\n');
    edges = get_edges();
    if (DEBUG_NORMAL || DEBUG_DUMP || DEBUG_DATA) {
        ip=get_n_edges_polys();
        printlog("edges:");
        for (i = 0; i < get_n_edges(); i++) {
            n=ip[i];
            /* report edge files and number of edges per file */
            printlog("%c %s (%d edge(s))",delim,get_edges_name(i),n); 
            if (! DEBUG_DATA)
            	delim=',';
            
            if (DEBUG_DATA) { /* report every line */
                printlog("\n");
                for (j = 0; j<n; j++) { /*over edges in a file */
                    print_poly_log(&(edges[i][j]));
                }
            }
        }
        if (!DEBUG_DATA) 
        	/* fputc('\n',logfile); */
        	printlog("\n");
    }
}

static void print_poly_log(POLYGON *edge) {
    int i,n=edge->lines;
    printlog("%+d",n);

    /* fputc('\n',a_stream); */
    printlog("\n");

    for (i=0;i<n;i++)
        printlog("%g %g\n",edge->p[i].x,edge->p[i].y);

    /* fflush(a_stream); */
    return;
}

#ifndef USING_R
POLYGON read_n_points(FILE *f, int np) {
	static char *line = NULL;
	static int line_size;
	POLYGON pl;
	int n = 0;
	char *cp;

	pl.lines = np; /* check later that first == last */
	pl.p = (PLOT_POINT *) emalloc(np * sizeof(PLOT_POINT));
	while (n < np) {
		get_line(&line, &line_size, f);
        /* pr_warning("file: %s, line '%s' %d\n",fname,line,n); */
		cp = strtok(line, DELIMITERS);
		if (cp == NULL) {
            if (feof(f)) {
                pr_warning("file: %s, not enough data points (%d requested, %d found '%s' )\n", fname, np, n, line);
                ErrMsg(ER_READ, fname);
            } else {
                pr_warning("file: %s, line without data encountered (%d)\n", fname, n+1);
                continue;
            }
		}
		if (read_double(cp, &(pl.p[n].x)))
			ErrMsg(ER_RDFLT, cp);
		cp = strtok(NULL, DELIMITERS);
		if (read_double(cp, &(pl.p[n].y)))
			ErrMsg(ER_RDFLT, cp);
		n++;
		cp = strtok(NULL, DELIMITERS);
		if (cp != NULL) { /* four on a line */
			if (read_double(cp, &(pl.p[n].x)))
				ErrMsg(ER_RDFLT, cp);
			cp = strtok(NULL, DELIMITERS);
			if (read_double(cp, &(pl.p[n].y)))
				ErrMsg(ER_RDFLT, cp);
			n++;
		}
	}
    pl.close = (pl.p[0].x == pl.p[np-1].x && pl.p[0].y == pl.p[np-1].y);
    setup_poly_minmax(&pl);
    
	return pl;
}
#endif

int point_in_polygon(PLOT_POINT point, POLYGON *pl)
{
    if ((point.x < pl->mbr.min.x) || (point.x > pl->mbr.max.x) ||
        (point.y < pl->mbr.min.y) || (point.y > pl->mbr.max.y))
        return 0;
    else
        return InPoly(point, pl) != 'o';
    
}

/*
This code is described in "Computational Geometry in C" (Second Edition),
Chapter 7.  It is not written to be comprehensible without the 
explanation in that book.

For each query point q, InPoly returns one of four char's:
	i : q is strictly interior to P
	o : q is strictly exterior to P
	v : q is a vertex of P
	e : q lies on the relative interior of an edge of P
These represent mutually exclusive categories.
For an explanation of the code, see Chapter 7 of 
"Computational Geometry in C (Second Edition)."

Written by Joseph O'Rourke, contributions by Min Xu, June 1997.
Questions to orourke@cs.smith.edu.
--------------------------------------------------------------------
This code is Copyright 1998 by Joseph O'Rourke.  It may be freely 
redistributed in its entirety provided that this copyright notice is 
not removed.
--------------------------------------------------------------------
*/

/*
InPoly returns a char in {i,o,v,e}.  See above for definitions.
*/

static char InPoly(PLOT_POINT q, POLYGON *Poly)
{
    int n = Poly->lines;
    PLOT_POINT *P=Poly->p;
    
    int	 i, i1;      /* point index; i1 = i-1 mod n */
    double x;          /* x intersection of e with ray */
    double xx=q.x, yy=q.y;
    int	 Rcross = 0; /* number of right edge/ray crossings */
    int    Lcross = 0; /* number of left edge/ray crossings */

    /* For each edge e=(i-1,i), see if crosses ray. */
    for( i = 0; i < n; i++ ) {
        /* First see if q=(0,0) is a vertex. */
        if (( P[i].x - xx )==0 &&( P[i].y - yy )==0 ) return 'v';
        i1 = ( i + n - 1 ) % n;
        /* printf("e=(%d,%d)\t", i1, i); */
    
        /* if e "straddles" the x-axis... */
        /* The commented-out statement is logically equivalent to the one 
           following. */
        /* if( ( ( P[i].y > 0 ) && ( P[i1].y <= 0 ) ) ||
           ( ( P[i1].y > 0 ) && ( P[i] .y <= 0 ) ) ) { */
    
        if( (( P[i].y - yy ) > 0 ) != (( P[i1].y - yy ) > 0 ) ) {
      
            /* e straddles ray, so compute intersection with ray. */
            x = (( P[i].x - xx) *( P[i1].y - yy ) -( P[i1].x - xx ) *( P[i].y - yy )) /
                (P[i1].y - P[i].y );
            /* printf("straddles: x = %g\t", x); */
      
            /* crosses ray if strictly positive intersection. */
            if (x > 0) Rcross++;
        }
        /* printf("Right cross=%d\t", Rcross); */
    
        /* if e straddles the x-axis when reversed... */
        /* if( ( ( P[i] .y < 0 ) && ( P[i1].y >= 0 ) ) ||
           ( ( P[i1].y < 0 ) && ( P[i] .y >= 0 ) ) )  { */
    
        if ( (( P[i].y - yy ) < 0 ) != (( P[i1].y - yy ) < 0 ) ) { 
      
            /* e straddles ray, so compute intersection with ray. */
            x = (( P[i].x - xx) *( P[i1].y - yy ) -( P[i1].x - xx ) *( P[i].y - yy )) /
                (P[i1].y - P[i].y);
            /* printf("straddles: x = %g\t", x); */

            /* crosses ray if strictly positive intersection. */
            if (x < 0) Lcross++;
        }
        /* printf("Left cross=%d\n", Lcross); */
    }	
  
    /* q on the edge if left and right cross are not the same parity. */
    if( ( Rcross % 2 ) != (Lcross % 2 ) )
        return 'e';
  
    /* q inside iff an odd number of crossings. */
    if( (Rcross % 2) == 1 )
        return 'i';
    else	return 'o';
}



#ifndef USING_R
/*------------------------------ EDGES ------------------------------*/
void read_edges(void) 
{
    int i,n;
    int n_edges /* ,n_edges_total=0 */ ;
    POLYGON **edges;
    int *n_edges_polys;
    
    n = get_n_edges();
    if (n <= 0)
    	return;
    /* get at local */
    /* edges = get_edges(); */
    edges = emalloc(n*sizeof(POLYGON*));
    /* n_edges_polys = get_n_edges_polys(); */
    n_edges_polys = emalloc(n*sizeof(int));

    /* ... read data in .. */
    for (i=0; i<n;  i++) { /* over all files */
        edges[i]=read_polygons(get_edges_name(i),&n_edges, NULL); 
/*         n_edges_total+=n_edges; */
        n_edges_polys[i]=n_edges;
    }
    /* ... and set back */
    set_edges(edges);
    set_n_edges_polys(n_edges_polys);
/*     set_n_edges(n_edges_total); */
    
}
#endif

void check_edges(DATA *d, const DPOINT *where) {
#define BAD FLT_MAX    
    /* Goes as follows : in the list d->sel set u.dist2 to
       BAD==FLT_MAX for the points failed through line-of-sight
       test. Then sort d->sel by u.dist2 and decrease d->n_sel by the
       number of points with u.dist2==BAD */
    
    int i,j,n_sel_edges,n;
    int p_i_p_result;
    char InPolyRes;
    int culled;
    unsigned int bad;
    int *n_edges_polys;
    POLYGON **edges;
    POLYGON *current_edge;
    PLOT_POINT p, p_where;
    static POLYGON **sel_edges = NULL;
    static int *point_edge_rel = NULL; 
    static int n_alloc_edges = 0;
    
	n = get_n_edges();
    if (n <= 0) 
		return;
    /* get at local */
    edges = get_edges();
    n_edges_polys = get_n_edges_polys();
    p_where.x=where->x;
    p_where.y=where->y;
    
    /* make a list of all relevant edges */
    n_sel_edges = 0;
    /* find relevant edges and also check if *where is inside/outside
       of edge for closed polynoms or which subdivision *where in for open edges */
    
    for (i=0; i<n; i++) { /* over files */
        for (j=0; j<n_edges_polys[i]; j++) {
            current_edge=&(edges[i][j]);
            
            /* check if edge extents are in vicinity of estimation point */
            /* basically we'll check it if point is in box+sel_rad (on both sides) */ 
            /* from check_edges of SPJ */
            
            if ((where->x > current_edge->mbr.min.x - d->sel_rad)       /*ltx */
                && (where->x  < current_edge->mbr.max.x + d->sel_rad)   /*rbx */
                && (where->y > current_edge->mbr.min.y - d->sel_rad)    /*rby */
                && (where->y  < current_edge->mbr.max.y + d->sel_rad))  /*lty*/
            {
                
#define IN_OUT(_a) ((_a) == 'i' ? POLY_IN : ((_a) == 'o' ? POLY_OUT : POLY_EDGE))

                if (current_edge->close) {
                    InPolyRes=InPoly(p_where, current_edge);
                    p_i_p_result=IN_OUT(InPolyRes);
                    /* don't use, if target on edge */
                    if (p_i_p_result == POLY_EDGE) continue;
                } else { /* do nothing at the moment, any ideas are welcome:
                            malakhanov@iwwl.rwth-aachen.de
                         */
                    p_i_p_result = POLY_IN;
                    
                }
                if (n_sel_edges+1>n_alloc_edges) {
                   n_alloc_edges = n_sel_edges+1;
                   sel_edges      = erealloc(sel_edges,n_alloc_edges*sizeof(POLYGON*));
                   point_edge_rel = erealloc(point_edge_rel,n_alloc_edges*sizeof(int)); 
                }
                
                sel_edges[n_sel_edges]=current_edge;
                point_edge_rel[n_sel_edges]=p_i_p_result;
                n_sel_edges++;       
            }
        }
    }
    if (n_sel_edges==0) 
		return; /* no relevant edges */

    culled=0;

    /* loop over all so far selected points */
    for (i=0; i<d->n_sel; i++) {
        p.x = d->sel[i]->x;
        p.y = d->sel[i]->y;
        
        for (j=0; j<n_sel_edges; j++) { /* over all selected edges */
            if ( sel_edges[j]->close) {
                /* InPolyRes=point_in_polygon(p,  sel_edges[j]); */
                InPolyRes = InPoly(p, sel_edges[j]);
                p_i_p_result = IN_OUT(InPolyRes);
                if (p_i_p_result == POLY_EDGE) continue; 
                bad = (p_i_p_result != point_edge_rel[j]);
            } else {
                /* check line-of-sight between data and target point for this edge */
                bad = check_open_edges(p, point_edge_rel[j], p_where, sel_edges[j]);
            }
            if (bad) {
                d->sel[i]->u.dist2 = BAD;
                culled++;
                break;
            }
        }
    } /* end of loop over all selected points */

    /* now look how bad (good) it was: */
    if (culled == 0) 
		return;
    /* otherwise sort d->sel by dist2 (u.dist2 is smaller then BAD, so it comes out first:*/
    qsort(d->sel, (size_t) d->n_sel, sizeof(DPOINT *),
          (int CDECL (*)(const void *,const void *)) dist_cmp);
    /* and set d->n_sel to the number of GOOD ones */
    d->n_sel -= culled;
    return;
}
#undef IN_OUT
#undef BAD

/* the rest is modified from SPJ: */
unsigned char line_of_sight(const PLOT_POINT data, const PLOT_POINT target,
                            const POLYGON *p_obj) {
    int i;
    unsigned int code;
    unsigned int crossings = 0;
    
	/* go through each segment to find out if it crosses line between data-target */
	for(i=0; i<p_obj->lines-1; i++) { /*we have really lines points */
        code = segment_cross_check(data, target, p_obj->p[i], p_obj->p[i+1], NULL, NULL);

        if (code==CROSS) { /* a simple intersection between data-target and the segment */
            crossings++;
            continue;
        }
        if ((code&SPOINT1) || (code&EPOINT1)) return 0; /* data or target at the edge
                                                           skip this edge completly */
        
	}
    return(crossings);
}

/*
  ----------------------------------------------------------------------
  http://www.cis.ohio-state.edu/hypertext/faq/usenet/graphics/algorithms-faq/faq.html
  Subject 1.03: How do I find intersections of 2 2D line segments?

    This problem can be extremely easy or extremely difficult depends
    on your applications.  If all you want is the intersection point,
    the following should work:

    Let A,B,C,D be 2-space position vectors.  Then the directed line
    segments AB & CD are given by:

        AB=A+r(B-A), r in [0,1]
        CD=C+s(D-C), s in [0,1]

    If AB & CD intersect, then

        A+r(B-A)=C+s(D-C), or

        Ax+r(Bx-Ax)=Cx+s(Dx-Cx)
        Ay+r(By-Ay)=Cy+s(Dy-Cy)  for some r,s in [0,1]

    Solving the above for r and s yields

            (Ay-Cy)(Dx-Cx)-(Ax-Cx)(Dy-Cy)
        r = -----------------------------  (eqn 1)
            (Bx-Ax)(Dy-Cy)-(By-Ay)(Dx-Cx)

            (Ay-Cy)(Bx-Ax)-(Ax-Cx)(By-Ay)
        s = -----------------------------  (eqn 2)
            (Bx-Ax)(Dy-Cy)-(By-Ay)(Dx-Cx)

    Let P be the position vector of the intersection point, then

        P=A+r(B-A) or

        Px=Ax+r(Bx-Ax)
        Py=Ay+r(By-Ay)

    By examining the values of r & s, you can also determine some
    other limiting conditions:

        If 0<=r<=1 & 0<=s<=1, intersection exists
            r<0 or r>1 or s<0 or s>1 line segments do not intersect

        If the denominator in eqn 1 is zero, AB & CD are parallel
        If the numerator in eqn 1 is also zero, AB & CD are coincident

    If the intersection point of the 2 lines are needed (lines in this
    context mean infinite lines) regardless whether the two line
    segments intersect, then

        If r>1, P is located on extension of AB
        If r<0, P is located on extension of BA
        If s>1, P is located on extension of CD
        If s<0, P is located on extension of DC

    Also note that the denominators of eqn 1 & 2 are identical.

    References:

    [O'Rourke (C)] pp. 249-51
[Gems III] pp. 199-202 "Faster Line Segment Intersection,"
*/
static unsigned int segment_cross_check(PLOT_POINT A, PLOT_POINT B, 
		PLOT_POINT C, PLOT_POINT D, double *R, double *S) {
    double uan,ubn,d,r,s;
    unsigned int code = 0;
/* calculates if two line segments cross  and checks if the crossing point is at 
   the start or endpoint of the segments */
    if (B.x == A.x && B.y == A.y) return code;
    
    d = (D.y-C.y)*(B.x-A.x) - (D.x-C.x)*(B.y-A.y);
  

    if (fabs(d)<FLT_MIN)  /* parallel OR coincident*/
        return segment_parallel_check(A,  B,  C,  D);
        
    uan = ((D.x-C.x)*(A.y-C.y) - (D.y-C.y)*(A.x-C.x)); /* for r */
    ubn = ((B.x-A.x)*(A.y-C.y) - (B.y-A.y)*(A.x-C.x)); /* for s */

    r=uan/d;
    s=ubn/d;
    if (R) *R=r;
    if (S) *S=s;
    
    /*  line segments do not intersect between (A,B) & (C,D)*/
    if ((r < 0.0) || (r > 1.0) || (s < 0.0) || (s > 1.0)) return code;
    
    code |= CROSS; /* lines cross somewhere between endpoints */
    
    if (r == 0.0 ) {
        /* intersection at start of segment 1 */
        code |= SPOINT1;
    } else if (r == 1.0) {
        /* cross at endpoint of line 1 */
        code |= EPOINT1;
    }


    if (s == 0.0) {
        /* intersection at start of segment 2 */
        code |= SPOINT2;
    } else if (s == 1.0) {
        /* cross at endpoint of line 2 */
        code |= EPOINT2;
    }
    
    return code;
    
}

/* check special case of lines through (A,B) and (C,D) being parallel */
static unsigned int segment_parallel_check(PLOT_POINT A, PLOT_POINT B, 
		PLOT_POINT C, PLOT_POINT D) {

    unsigned int code = 0;
    /* remember, A- data point, B- target, C - segment start, D -segment end */
    /* first check A,B,C on one line or parallel */
    if ((B.x - A.x)*(C.y - A.y) != (C.x - A.x)*(B.y - A.y))
        return code; /* not on one line, i.e. are parallel */

    if ( segment_between_check( C, D, A ) )  /* data  at polyline */
        code |= SPOINT1;
    
    if ( segment_between_check( C, D, B ) )  /* target  at polyline */
        code |= EPOINT1;
    
    if ( segment_between_check( A, B, C) )  /* data  at polyline */
        code |= SPOINT2;
    
    if ( segment_between_check( A, B, D ) )  /* target  at polyline */
        code |= EPOINT2;
    
    if (code) code|=CROSS;
    
    return code;
    
    
}

/* check the c is between a and b */
static int segment_between_check(PLOT_POINT a, PLOT_POINT b, PLOT_POINT c) {
    if ( a.x != b.x )
        return ((a.x <= c.x) && (c.x <= b.x)) ||
            ((a.x >= c.x) && (c.x >= b.x));
    else
        return ((a.y <= c.y) && (c.y <= b.y)) ||
            ((a.y >= c.y) && (c.y >= b.y));
}

/*
  ----------------------------------------------------------------------
  http://www.cis.ohio-state.edu/hypertext/faq/usenet/graphics/algorithms-faq/faq.html
  Subject 1.02: How do I find the distance from a point to a line?


    Let the point be C (Cx,Cy) and the line be AB (Ax,Ay) to (Bx,By).
    The length of the line segment AB is L:

        L= sqrt( (Bx-Ax)^2 + (By-Ay)^2 ) .

    Let P be the point of perpendicular projection of C onto AB.
    Let r be a parameter to indicate P's location along the line
    containing AB, with the following meaning:

          r=0      P = A
          r=1      P = B
          r<0      P is on the backward extension of AB
          r>1      P is on the forward extension of AB
          0<r<1    P is interior to AB

    Compute r with this:

            (Ay-Cy)(Ay-By)-(Ax-Cx)(Bx-Ax)
        r = -----------------------------
                        L^2

    The point P can then be found:

        Px = Ax + r(Bx-Ax)
        Py = Ay + r(By-Ay)

    And the distance from A to P = r*L.

    Use another parameter s to indicate the location along PC, with the 
    following meaning:
           s<0      C is left of AB
           s>0      C is right of AB
           s=0      C is on AB

    Compute s as follows:

            (Ay-Cy)(Bx-Ax)-(Ax-Cx)(By-Ay)
        s = -----------------------------
                        L^2


    Then the distance from C to P = s*L.


*/

static int CDECL dist_cmp(const DPOINT **pa, const DPOINT **pb) {
/* ANSI qsort() conformant dist_cmp */

	if ( (*pa)->u.dist2 < (*pb)->u.dist2 )
		return -1;
	if  ( (*pa)->u.dist2 > (*pb)->u.dist2 )
		return 1;
    return 0;
    
} 

static unsigned int check_open_edges(const PLOT_POINT data, const int target_side,
                                     const PLOT_POINT target, const POLYGON *p_obj) {
    
    int los=line_of_sight(data,target,p_obj);
    
    return (los % 2);
}


void setup_poly_minmax(POLYGON *pl) {
    int i, n=pl->lines;
    double minx,maxx,miny,maxy;
    
    minx=miny=DBL_MAX;
    maxx=maxy=-DBL_MAX;
    
    for (i=0;i<n;i++) {
        minx = MIN(minx, pl->p[i].x);
        miny = MIN(miny, pl->p[i].y);
        maxx = MAX(maxx, pl->p[i].x);
        maxy = MAX(maxy, pl->p[i].y);
    }
    pl->mbr.min.x = minx;
    pl->mbr.min.y = miny;
    pl->mbr.max.x = maxx;
    pl->mbr.max.y = maxy;
}


#undef CROSS  
#undef SPOINT1
#undef SPOINT2
#undef EPOINT1
#undef EPOINT2
#undef POLY_IN 
#undef POLY_OUT
#undef POLY_EDGE
#undef IWW
