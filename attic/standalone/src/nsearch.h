#ifndef SEARCH_H
# define SEARCH_H /* avoid multiple inclusion */

void qtree_free(QTREE_NODE *node);
void qtree_pop_point(DPOINT *p, DATA *d);
void qtree_push_point(DATA *d, DPOINT *p);
void qtree_rebuild(DATA *d);
int qtree_select(DPOINT *where, DATA *d);
/* 2-norm distances from point to block: */
double pb_norm_3D(const DPOINT *where, BBOX bbox);
double pb_norm_2D(const DPOINT *where, BBOX bbox);
double pb_norm_1D(const DPOINT *where, BBOX bbox);

/* define the maximum depth of the quadtree; 
 * Fri Jul  4 12:05:47 CEST 2003
 * if this is not defined, more than gl_split points at
 * a single spatial location cause infinite recursion
 * 10 seems a reasonable value: 1/2048 of the bbox dim
 * */
#define MAX_RECURSION_DEPTH 11

#endif /* SEARCH_H */
