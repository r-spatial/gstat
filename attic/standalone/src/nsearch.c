/*
   This software module is copyright 1997 (c):

   Steve Joyce                            mailto:steve.joyce@resgeom.slu.se
   Remote Sensing Laboratory                 http://www-umea.slu.se/~fjasj/
   Dept. of Forest Resource Mgmt. & Geomatics          Tel: +46 90 16 69 57
   Swedish University of Agricultural Sciences         Fax: +46 90 14 19 15 
   S-901 83 Umeå, Sweden

   Distributed freely under the terms of the GNU General Public License
   as a component of:

   Gstat, a program for geostatistical modelling, prediction and simulation
   Copyright 1992-2009 (C) Edzer J. Pebesma

   Edzer J. Pebesma (E.Pebesma@geo.uu.nl)
   Landscape and environmental research group
   Faculty of geographical sciences
   University of Amsterdam
   Nieuwe Prinsengracht 130
   1018 VZ  Amsterdam -- The Netherlands
*/

/*
 *    qtree.c: quick neighbourhood selection routines for gstat
 */

/*
 * Edzer's CHANGE LOG from Steve's original contribution:
 * converted most float to double;
 * bbox size initialization -1.0
 * added some ErrMsg() error message checks
 * changed qtree_quick_select() call (DPOINTX disappeared)
 * made is_leaf(node) a macro; removed max_ppnode -> Q_SPLIT_AT
 * added flexible 1/2/3d tree support: detect from data->mode;
 * added bbox.mode field;
 * modified routines like in_bbox, sub_bbox to detect tree dimension
 * qtree_free(): node->n_node should be N_NODES(node)
 * init_qtree is now called from qtree_quick_select()
 * added bbox_from_* routines
 * added get_index macro: searching is not required
 * added qtree_print() functions
 * might have missed some.
 * (after 2.0:)
 * Q_SPLIT_AT -> gl_split
 * implemented the priority queue mechanism, see PR bucket quad trees
 * demonstrated at http://www.cs.umd.edu/~brabec/quadtree/index.html.
 * (this enables efficient search when only max is defined, i.e. without
 *  a radius -- especially in case of simulation this seems really 
 *  efficient)
 * removed min_bbox -->> this became obsolete with the priority queue
 * 
 * I guess many of my modifications assume you have sufficient memory,
 * and that you want to use it to speed things up.
 */

#include <stdio.h>
#include <stdlib.h> /* qsort() */
#include <math.h>  /* sqrt() */
#include <float.h> /* DBL_MAX */
#include <limits.h> /* INT_MAX */

#include "defs.h"
#include "userio.h"
#include "debug.h"
#include "data.h"
#include "utils.h"
#include "glvars.h" /* get_method(), name_identifier() */
#include "mapio.h"
#include "predict.h" /* get_mask0() */
#include "nsearch.h"
#include "pqueue.h" /* include _after_ search.h! */

#define N_NODES(x) (x==NULL?0:(-((x)->n_node)))
	/* a negative n_node means it's not a leaf */

/*
 * this is a 'tuning' parameter defining the number of points in 
 *  a quad before it is split into 4
 *           #define Q_SPLIT_AT 4
 */

/*
 * another tuning parameter defining the minimum quad size allowed.
 * (multiplied by d->sel_rad)
   #define Q_STOP_AT 0.5 (now obsolete)
 */

#define is_leaf(node) (((node) == NULL) || (((node)->n_node) >= 0))
#define is_qtree_search(d) (d->qtree_root != NULL)
#define get_index(pt, bb) \
	((pt->x >= bb.x + 0.5 * bb.size) | \
	(((bb.mode & Y_BIT_SET) && (pt->y >= bb.y + 0.5 * bb.size)) << 1) | \
	(((bb.mode & Z_BIT_SET) && (pt->z >= bb.z + 0.5 * bb.size)) << 2)) 

static void init_qtree(DATA *d);
static void init_qnode(QTREE_NODE **p_node, int is_leaf, BBOX bb);
static void qtree_push(DPOINT *p, QTREE_NODE **p_node, int recursion_depth);
static BBOX sub_bbox(const BBOX bb, int i);
static int in_bbox(const DPOINT *p, BBOX bb);
static void qtree_split_node(QTREE_NODE *node, BBOX bb, int rec_level);
static QTREE_NODE *qtree_expand(const DPOINT *p, QTREE_NODE *root);
static QTREE_NODE **qtree_find_node(const DPOINT *p, QTREE_NODE **p_node,
	BBOX *bbox);
void qtree_print(DATA *d);
static void logprint_qtree(QTREE_NODE *node, int depth);
static BBOX bbox_from_grid(const GRIDMAP *gt, const DATA_GRIDMAP *dg);
static BBOX bbox_from_data(DATA *d);
static void qtree_zero_all_leaves(QTREE_NODE *node);
static int CDECL node_cmp(const QUEUE_NODE *a, const QUEUE_NODE *b);
void logprint_queue(QUEUE *queue);

static void init_qtree(DATA *d) {
/*
 * Initialize the root of the search tree:
 * Since this is called from qtree_quick_select(), we know allready
 * quite a lot. This helps choosing sensible values for the 
 * top level bbox origin and size.
 */
 	const GRIDMAP *gt = NULL;
 	DATA *simlocs = NULL;
 	int i, mode;
 	BBOX bbox;

	if (is_simulation(get_method())) {
		/*
		 * sequential simulation: the simulation path (through DATA or GRIDMAP)
		 * will make up for (most of) the search locations:
		 */
		gt = (const GRIDMAP *) get_mask0();
		simlocs = get_dataval();
		/* in case of simulation one of them will be non-NULL */
	} 

	/* get initial estimate for top level bbox: */
	if (gt != NULL)
		bbox = bbox_from_grid(gt, NULL);
	else if (simlocs != NULL) {
		bbox = bbox_from_data(simlocs);
		if (bbox.size <= 0.0)
			bbox = bbox_from_data(d);
	} else
		bbox = bbox_from_data(d);
	init_qnode(&(d->qtree_root), d->n_list < gl_split, bbox); /* ML1 */
	
	mode = bbox.mode;

	for (i = 0; i < d->n_list; i++)
		qtree_push_point(d, d->list[i]); /* now they won't be rejected */ /* ML2 */

	if (DEBUG_DUMP) {
		printlog("top level search tree statistics for data(%s):\n",
			name_identifier(d->id));
		printlog("bounding box origin [");
		if (mode & X_BIT_SET)
			printlog("%g", d->qtree_root->bb.x);
		if (mode & Y_BIT_SET)
			printlog(",%g", d->qtree_root->bb.y);
		if (mode & Z_BIT_SET)
			printlog(",%g", d->qtree_root->bb.z);
		printlog("]; dimension %g\n", d->qtree_root->bb.size);
	}
	/* qtree_print(d); */
	return;
}

static void init_qnode(QTREE_NODE **p_node, int isleaf, BBOX bb) {
/*
 * initialize a node in the search tree.
 * bb.mode tells the dimension of the tree (1D/2D/3D)
 */
	int i;

	if (*p_node == NULL) {
		*p_node = (QTREE_NODE *) emalloc(sizeof(QTREE_NODE)); /* ML1 */
		(*p_node)->bb = bb;
	}

	if (isleaf)
		(*p_node)->n_node = 0;
	else {
		if (bb.mode & Z_BIT_SET) /* 3D: oct-tree */
			(*p_node)->n_node = -8;
		else if (bb.mode & Y_BIT_SET) /* 2D: quad-tree */
			(*p_node)->n_node = -4;
		else if (bb.mode & X_BIT_SET) /* 1D: binary tree */
			(*p_node)->n_node = -2;
		else /* no x/y/z bit set...??? */
			ErrMsg(ER_IMPOSVAL, "init_qnode: invalid mode");
		(*p_node)->u.node = (QTREE_NODE **)
				emalloc(N_NODES(*p_node) * sizeof(QTREE_NODE *));
		for (i = 0; i < N_NODES(*p_node); i++)
			(*p_node)->u.node[i] = NULL;
	}
	return;
}

void qtree_push_point(DATA *d, DPOINT *where) {
/* add a single point to the search tree structure in d->qtree_root */

	/*
	 * do not do this while reading the data: suppose we'll never need
	 * a neighbourhood selection after all! 
	 */
	if (! is_qtree_search(d))
		return;

	/* min_bbox = d->sel_rad * Q_STOP_AT; */

	/*
	 * if this point is outside the current search tree,
	 * we need to add another level to the top and try again:
	 */
	while (! in_bbox(where, d->qtree_root->bb))
		d->qtree_root = qtree_expand(where, d->qtree_root);

	/*
	 * finally push the point onto the tree:
	 */
	qtree_push(where, &(d->qtree_root), 0);
	return;
}

static void qtree_push(DPOINT *where, QTREE_NODE **p_node, 
				int recursion_depth) {
/* add a data point to the quad tree starting at the node specified. */

	QTREE_NODE **p_leaf, *node;
	BBOX bb;

	bb = (*p_node)->bb;
	recursion_depth += 1;
	/* printf("recursion_depth: %d, max %d\n", recursion_depth, MAX_RECURSION_DEPTH); */
	/* find the leaf node where this point belongs */
	p_leaf = qtree_find_node(where, p_node, &bb);

	if (*p_leaf == NULL)
		init_qnode(p_leaf, 1, bb); /* leaf == 1: sets ->n_node to 0 */

	node = *p_leaf;

	/* If it is already full, split it into another level and try again: */
	if (node->n_node == gl_split && recursion_depth < MAX_RECURSION_DEPTH) {
		qtree_split_node(node, (*p_node)->bb, recursion_depth); 
		qtree_push(where, &node, recursion_depth);
		return;
	}

	/* XXX */
	if (node->n_node == 0)
		node->u.list = (DPOINT **) emalloc(sizeof(DPOINT *));
	else
		node->u.list = (DPOINT **) erealloc(node->u.list,
			(node->n_node + 1) * sizeof(DPOINT *)); /* ML2 */
	node->u.list[node->n_node] = where;
	node->n_node++;
	return;
}

void qtree_pop_point(DPOINT *where, DATA *d) {
	int i;
	QTREE_NODE *node, **p_node;
/* delete a point from the search tree */

	if (! is_qtree_search(d)) /* don't bother */
		return;

	p_node = qtree_find_node(where, &(d->qtree_root), NULL);

	if (*p_node == NULL)
		ErrMsg(ER_IMPOSVAL, "qtree_pop_point(): could not find node");

	node = *p_node;

	for (i = 0; i < node->n_node; i++) {
		if (where == node->u.list[i]) {
			/* don't preserve order: copy last to this one: */
			node->u.list[i] = node->u.list[node->n_node - 1];
			break; /* from for loop */
		}
	}
	node->n_node--;

	/* free memory if empty list: */
	if (node->n_node == 0) {
		efree(node->u.list);
		efree(node);
		*p_node = NULL;
	} 

	return;
}

void qtree_free(QTREE_NODE *node) {
/*
 * If a push or search fails, you might want to consider getting rid of 
 * whole tree and default to exhaustive search. (SJ)
 * [If you ever get so far, exhaustive search will take
 * a nearly infinite amount of time. Instead, tweek gl_split. --EJP]
 */
	int i;

	if (node == NULL)
		return;

	if (!is_leaf(node)) {
		for (i = 0; i < N_NODES(node); i++)
			qtree_free(node->u.node[i]);
		efree(node->u.node);
	} else
		efree(node->u.list);
	efree(node);
	return;
}

static void qtree_split_node(QTREE_NODE *node, BBOX bbox, int rec_level) {
/*
 * split the quadtree at 'node' and redistribute its points
 */
	DPOINT **list;
	int i, n;

	/* first copy the points to a temporary location and free the pointers */
	n = node->n_node;
	list = node->u.list; /* save temporary copy */

	/* make a node from this leaf, overwrite u: */
	init_qnode(&node, 0, bbox);

	/* redistribute the points into the child nodes where they belong */
	for (i = 0; i < n; i++)
		qtree_push(list[i], &node, rec_level);
	efree(list); 
	return;
}


static QTREE_NODE *qtree_expand(const DPOINT *where, QTREE_NODE *root) {
/*
 * expand the top level of the search tree
 */
	QTREE_NODE *new_top = NULL;
	DPOINT old_centre;
	BBOX old_bb, new_bb;
	int i;

	old_bb = root->bb;
	old_centre.x = old_centre.y = old_centre.z = 0.0;
	
	if (old_bb.mode & X_BIT_SET)
		old_centre.x = old_bb.x + old_bb.size / 2.0;
	if (old_bb.mode & Y_BIT_SET)
		old_centre.y = old_bb.y + old_bb.size / 2.0;
	if (old_bb.mode & Z_BIT_SET)
		old_centre.z = old_bb.z + old_bb.size / 2.0;

	new_bb = old_bb;
/* 
 * set the new bounding box: Steve, could you check this?
 * (I didn't grasp your original bbox setting here:)
 *
 * set the root bbox to the new_top bbox:
 */
	if ((old_bb.mode & X_BIT_SET) && (where->x < old_bb.x))
		new_bb.x -= old_bb.size;

	if ((old_bb.mode & Y_BIT_SET) && (where->y < old_bb.y))
		new_bb.y -= old_bb.size;

	if ((old_bb.mode & Z_BIT_SET) && (where->z < old_bb.z))
		new_bb.z -= old_bb.size;

	new_bb.size *= 2.0; 

	/* link the old root node to the proper spot in new_top: */
	init_qnode(&new_top, 0, old_bb);
	i = get_index((&old_centre), new_bb);
	new_top->u.node[i] = root;
	new_top->bb = new_bb;

	/* make this one the new root */
	return new_top;
}

static QTREE_NODE **qtree_find_node(const DPOINT *where, QTREE_NODE **p_node,
	BBOX *bb) {
/*
 * find the deepest leaf (end node) in the tree that bounds this point's
 * coordinates.
 * It's bounding box is saved in the location pointed to by p_bbox
 */
	int i;

	if (is_leaf(*p_node))
		return p_node;

	/* find in which node we are: */
	i = get_index(where, (*p_node)->bb);
	if (bb != NULL)
		*bb = sub_bbox(*bb, i);
	/* recurse into this node: */
	return qtree_find_node(where, &((*p_node)->u.node[i]), bb);
}


static BBOX sub_bbox(const BBOX bbox, int index) {
/*
 * return the bounding box of a quad-tree child node based on the index
 * layout of octree index:
 * 
 *          | dz <= 0    dz > 0 
 *  --------|--------------------
 *  dy  > 0 |  2  3       6  7
 *  dy <= 0 |  0  1       4  5
 *          |
 *  dx ? 0  |  <= >       <= >
 */
 	double size;
 	BBOX b;

	b = bbox;
	b.size = size = bbox.size / 2.0;

	if (index & X_BIT_SET) /* 1, 3, 5, 7 */
		b.x += size;
	if (index & Y_BIT_SET) /* 2, 3, 6, 7 */
		b.y += size;
	if (index & Z_BIT_SET) /* 4, 5, 6, 7 */
		b.z += size;
	return b;
}

static int in_bbox(const DPOINT *where, BBOX bbox) {
/*
 * check if where is inside the bounding box:
 * _on_ the left/lower/downside boundary or is inside the box
 */

	if ((bbox.mode & X_BIT_SET)  &&
		((where->x < bbox.x) || (where->x >= bbox.x + bbox.size)))
			return 0;
	if ((bbox.mode & Y_BIT_SET) &&
		((where->y < bbox.y) || (where->y >= bbox.y + bbox.size)))
			return 0;
	if ((bbox.mode & Z_BIT_SET) &&
		((where->z < bbox.z) || (where->z >= bbox.z + bbox.size)))
			return 0;
	/* so, where apparently in ... */
	return 1;
}

/*
 * pb_norm2_?D() functions: 
 * calculate shortest (squared) euclidian distance from a point to a BBOX,
 * for ? being 1, 2 or 3 dimensions
 */
double pb_norm_1D(const DPOINT *where, BBOX bbox) {
	double x, dx;

	x = where->x;
	if (x < bbox.x) {
		dx = bbox.x - x;
		return dx * dx;
	}
	bbox.x += bbox.size;
	if (x > bbox.x) {
		dx = x - bbox.x;
		return dx * dx;
	}
	return 0.0; /* inside box */
}

double pb_norm_2D(const DPOINT *where, BBOX bbox) {
	double x, y, dx = 0.0, dy = 0.0;

	x = where->x;
	y = where->y;
	if (x < bbox.x)
		dx = bbox.x - x;
	else {
		bbox.x += bbox.size;
		if (x > bbox.x)
			dx = x - bbox.x;
	}
	if (y < bbox.y)
		dy = bbox.y - y;
	else {
		bbox.y += bbox.size;
		if (y > bbox.y)
			dy = y - bbox.y;
	}
	return dx * dx + dy * dy;
}

double pb_norm_3D(const DPOINT *where, BBOX bbox) {
	double x, y, z, dx = 0.0, dy = 0.0, dz = 0;

	x = where->x;
	y = where->y;
	z = where->z;
	if (x < bbox.x)
		dx = bbox.x - x;
	else {
		bbox.x += bbox.size;
		if (x > bbox.x)
			dx = x - bbox.x;
	}
	if (y < bbox.y)
		dy = bbox.y - y;
	else {
		bbox.y += bbox.size;
		if (y > bbox.y)
			dy = y - bbox.y;
	}
	if (z < bbox.z)
		dz = bbox.z - z;
	else {
		bbox.z += bbox.size;
		if (z > bbox.z)
			dz = z - bbox.z;
	}
	return dx * dx + dy * dy + dz * dz;
}

void qtree_print(DATA *d) {
/*
 * plot the full tree (2D), in a format that can be read by jgraph, found
 * at netlib or at http://kenner.cs.utk.edu/~plank/plank/jgraph/jgraph.html
 */
	printlog("newgraph\nxaxis size 3\nyaxis size 3\n");
	printlog("title : %s [n = %d]\n",
		name_identifier(d->id), d->n_list);
	logprint_qtree(d->qtree_root, 0);
	return;
}

static void logprint_qtree(QTREE_NODE *node, int depth) {
	BBOX b;
	int i;

	if (node == NULL)
		return;
	b = node->bb;
	if (!is_leaf(node)) {
		printlog("newline linethickness 0.3 pts %g %g %g %g %g %g %g %g %g %g\n",
			b.x, b.y, b.x+b.size, b.y, b.x+b.size,
			b.y+b.size, b.x, b.y+b.size, b.x, b.y);
		for (i = 0; i < N_NODES(node); i++)
			logprint_qtree(node->u.node[i], depth+1);
	} else {
		printlog("newline pts %g %g %g %g %g %g %g %g %g %g\n",
			b.x, b.y, b.x+b.size, b.y, b.x+b.size,
			b.y+b.size, b.x, b.y+b.size, b.x, b.y);
		/*
		if (node == NULL)
			printlog("newcurve marktype circle fill 1 pts %g %g\n",
				b.x+0.5*b.size, b.y+0.5*b.size);
		*/
		if (node->n_node > 0) {
			printlog("newcurve marktype cross pts");
			for (i = 0; i < node->n_node; i++)
				printlog(" %g %g",	node->u.list[i]->x, node->u.list[i]->y);
			printlog("\n");
		}
	}
}

static BBOX bbox_from_grid(const GRIDMAP *gt, const DATA_GRIDMAP *dg) {
/* derive a sensible top level bounding box from grid map topology */
	double sizex, sizey;
	BBOX bbox;

	if (gt) {
		sizex = gt->cols * gt->cellsizex;
		sizey = gt->rows * gt->cellsizey;
		bbox.x = gt->x_ul;
		bbox.y = gt->y_ul - sizey;
		/*
	 	 * bbox.size should be set to such a value that the smallest
	 	 * leaf fits exactly over 4 grid map cells 
	 	 */
		bbox.size = MIN(gt->cellsizex, gt->cellsizey);
	} else {
		sizex = dg->cols * dg->cellsizex;
		sizey = dg->rows * dg->cellsizey;
		bbox.x = dg->x_ul;
		bbox.y = dg->y_ul - sizey;
		bbox.size = MIN(dg->cellsizex, dg->cellsizey);
	}
	bbox.z = DBL_MAX;
	while (bbox.size < MAX(sizex, sizey))
		bbox.size *= 2;
	bbox.mode = (X_BIT_SET | Y_BIT_SET); /* i.e. 3 */
	return bbox;
}

static BBOX bbox_from_data(DATA *d) {
/* derive a sensible top level bounding box from a data var */
	double maxspan, dy, dz;
	BBOX bbox;

	if (d->grid)
		return bbox_from_grid(NULL, d->grid);
	bbox.x = d->minX;
	bbox.y = d->minY;
	bbox.z = d->minZ;
	bbox.mode = d->mode; /* ??? */
	/*
	bbox.mode = d->mode & (X_BIT_SET|Y_BIT_SET|Z_BIT_SET);
	maxspan = MAX((d->maxX-d->minX), MAX((d->maxY-d->minY),(d->maxZ-d->minZ)));
	*/
	maxspan = fabs(d->maxX - d->minX);
	dy = fabs(d->maxY - d->minY);
	if (dy > maxspan)
		maxspan = dy;
	dz = fabs(d->maxZ - d->minZ);
	if (dz > maxspan)
		maxspan = dz;
		
	/* with d->grid_size entered by user:
	if (d->grid_size > 0.0) {
		bbox.x -= 0.5 * d->grid_size;
		bbox.y -= 0.5 * d->grid_size;
		bbox.z -= 0.5 * d->grid_size;
		bbox.size = d->grid_size;
		do {
			bbox.size *= 2;
		} while (bbox.size < (maxspan + d->grid_size));
	} 
	*/
	bbox.size = maxspan * 1.01;
	return bbox;
}

static void qtree_zero_all_leaves(QTREE_NODE *node) {
	int i;

	if (!is_leaf(node)) {
		for (i = 0; i < N_NODES(node); i++)
			qtree_zero_all_leaves(node->u.node[i]);
	} else if (node != NULL)
		node->n_node = 0;
	return;
}

void qtree_rebuild(DATA *d) {
/* rebuild tree */
	int i;
	QTREE_NODE **p_leaf, *leaf;

	if (d->n_list <= 0 || d->qtree_root == NULL)
		return;

	qtree_zero_all_leaves(d->qtree_root);
	for (i = 0; i < d->n_list; i++) {
		p_leaf = qtree_find_node(d->list[i], &(d->qtree_root), NULL);
		leaf = *p_leaf;
		leaf->u.list[leaf->n_node] = d->list[i];
		leaf->n_node++;
	}
	return;
}

static DPOINT *get_nearest_point(QUEUE *q, DPOINT *where, DATA *d) {
/*
 * returns the first (closest) DPOINT in the priority queue q, after all
 * unwinding necessary (which is effectively recursion into the tree).
 *
 * this and the following functions: Copyright (GPL) 1998 Edzer J. Pebesma
 */
	QUEUE_NODE head, *el = NULL /* temporary storage */ ;
	QTREE_NODE *node;
	int i, n;

	while (q->length > 0) {  /* try: */
		/* logprint_queue(q); */
		head = dequeue(q);
		if (! head.is_node) { /* nearest element is a point: */
			if (el != NULL)
				efree(el);
			return head.u.p;
		}
		node = head.u.n;
		if (is_leaf(node)) { /* ah, the node dequeued is a leaf: */
			/* printf("node->n_node: %d\n", node->n_node); */
			if (node->n_node > 0)
				el = (QUEUE_NODE *) erealloc(el, node->n_node * sizeof(QUEUE_NODE));
			for (i = 0; i < node->n_node; i++) { /* enqueue it's DPOINT's: */
				el[i].is_node = 0;
				el[i].u.p = node->u.list[i];
				el[i].dist2 = node->u.list[i]->u.dist2 =
						d->pp_norm2(where, node->u.list[i]);
			}
			n = node->n_node;
		} else { /* nope, but enqueue its sub-nodes: */
			if (N_NODES(node) > 0)
				el = (QUEUE_NODE *) erealloc(el, N_NODES(node) * sizeof(QUEUE_NODE));
			for (i = n = 0; i < N_NODES(node); i++) {
				if (node->u.node[i] != NULL) {
					el[n].is_node = 1;
					el[n].u.n = node->u.node[i];
					el[n].dist2 = d->pb_norm2(where, node->u.node[i]->bb);
					n++;
				}
			}
		}
		if (n > 0)
			enqueue(q, el, n);
	}
	/* the while-loop terminates when the queue is empty */
	if (el != NULL)
		efree(el);
	return NULL;
}

void logprint_queue(QUEUE *queue) {
	Q_ELEMENT *q;
	QUEUE_NODE *e;

	printlog("current priority queue size: %d\n", queue->length);
	for (q = queue->head; q != NULL; q = q->next) {
		e = &(q->el);
		printlog("%s %12.6g", e->is_node ?
			"Node at " : "Point at", sqrt(e->dist2));
		if (e->is_node)
			printlog(" [xll=%g,yll=%g,size=%g] (with %d %s)\n",
				e->u.n->bb.x, e->u.n->bb.y,
				e->u.n->bb.size, ABS(e->u.n->n_node),
				e->u.n->n_node < 0 ? "nodes" : "points");
		else
			printlog(" [index %d, value %g]\n",
				GET_INDEX(e->u.p), e->u.p->attr);
	}
}

static int CDECL node_cmp(const QUEUE_NODE *a, const QUEUE_NODE *b) {
/* ANSI qsort() conformant comparison function */

	if (a->dist2 < b->dist2)
		return -1;
	if (a->dist2 > b->dist2)
		return 1;
	/* equal distances: prefer DPOINT over a node to speed up things */
	if (a->is_node != b->is_node)
		return (a->is_node ? 1 : -1);
	return 0;
}

int qtree_select(DPOINT *where, DATA *d) {
	DPOINT *p;
	static QUEUE *q = NULL;
	static QUEUE_NODE root;
	int sel_max;
	double rad2;

	/* if this is the first time calling with this d: */
	if (d->qtree_root == NULL)
		init_qtree(d);

	root.is_node = 1;
	root.u.n = d->qtree_root;
	root.dist2 = 0.0;

	q = init_queue(q, node_cmp);
	enqueue(q, &root, 1);

	if (d->sel_rad >= DBL_MAX) {
		/*
		 * simply get the d->sel_max nearest:
		 */
		for (d->n_sel = 0; d->n_sel < d->sel_max; d->n_sel++)
			d->sel[d->n_sel] = get_nearest_point(q, where, d);
	} else {
		/*
		 * also consider a maximum distance to where
		 */
		if (d->vdist) /* select everything within sel_rad; cut later on */
			sel_max = INT_MAX;
		else 
			sel_max = d->sel_max;
		rad2 = d->sel_rad * d->sel_rad;
		d->n_sel = 0;
		while (d->n_sel < sel_max) {
			p = get_nearest_point(q, where, d);
			if (p != NULL && p->u.dist2 <= rad2) { /* accept this point */
				d->sel[d->n_sel] = p;
				d->n_sel++;
			} else
				break; /* reject, and break while loop */
		}
		if (d->n_sel < d->sel_min) {
			/*
			 * d->sel_min was set: consider points beyond radius
			 */
			if (d->force) /* proceed beyond d->sel_max */
				for ( ; d->n_sel < d->sel_min; d->n_sel++)
					d->sel[d->n_sel] = get_nearest_point(q, where, d);
			else /* stop: a zero d->n_sel will result in a missing value */
				d->n_sel = 0;
		}
	}
	return d->n_sel;
}
