%{
/*
    Gstat, a program for geostatistical modelling, prediction and simulation
    Copyright 1992, 2003 (C) Edzer J. Pebesma

    Edzer J. Pebesma, e.pebesma@geo.uu.nl
    Department of physical geography, Utrecht University
    P.O. Box 80.115, 3508 TC Utrecht, The Netherlands

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

/*
 * parse.y: LALR(1) grammar for the gstat command syntax.
 * to make parse.c, type ``make parse.c'', it will use bison or yacc.
 *
 * If you fail (or don't have bison or yacc), then copy the file parse.c_
 * to parse.c. All this can be prevented by running configure while NO_YACC
 * is defined.
 * 
 * The parser assumes that in the function yylex() each identifier is
 * duplicated to ylval.sval, not just a pointer-copy. (some memory loss
 * will occur as a result)
 *
 * hints to extend the parser: 
 * o add a command: copy all from the most similar available command
 *   (add a %token and %type declaration, add a rule, add a return value
 *   from yylex() -> see the IDENT actions in lex.l)
 * o add a ``set'' variable: modify is_set_expr(), glvars.[ch] and defaults.h
 * o add a data() command: modify data.[ch] and is_data_expr()
 * o add a variogram model: vario*.[ch]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "defs.h"

#ifdef HAVE_UNISTD_H
# include <unistd.h> /* isatty() */
#endif

#include "data.h"
#include "vario.h"
#include "debug.h"
#include "glvars.h"
#include "userio.h"
#include "utils.h"
#include "lex.h"

static DATA *d = NULL, **dpp = NULL;
static DPOINT *bp = NULL;
static VARIOGRAM *v = NULL;
static int id = -1, id1 = -1, id2 = -1, col1 = -1, col2 = -1,
	fit_sill = 0, fit_range = 0, nrangepars = 1,
	vector_only = 0, allow_vector_only = 0;
static double range[NRANGEPARS], anis[5];
static char **ofn = NULL, *boundary_file = NULL;
static VARIOGRAM *parse_variogram = NULL;
static D_VECTOR *sd_vector = NULL;

#ifdef YYBISON
# ifndef __STDC__
#  define __STDC__
/* or else all const's will be defined empty */
# endif
#endif

typedef struct {
	const char *name;
	void *ptr;
	enum { 
		UNKNOWN, 
		IS_INT, 
		IS_UINT, 
		IS_REAL, 
		IS_STRING, 
		IS_D_VECTOR, 
		NO_ARG 
	} what;
	enum { 
		NOLIMIT, 
		GEZERO, 
		GTZERO 
	} limit;
} GSTAT_EXPR;

GSTAT_EXPR expr = { NULL, NULL, UNKNOWN, NOLIMIT };

static void push_data_X(DATA *d, int id);
static int is_data_expr(DATA *d, GSTAT_EXPR *expr, const char *fld);
static int is_set_expr(GSTAT_EXPR *expr, const char *fld);
static int is_block_expr(GSTAT_EXPR *expr, const char *s);
static void push_marginal(char *name, double val);
static void check_assign_expr(GSTAT_EXPR *expr);
static void reset_parser(void);
static void verify_data(DATA *d);

#define gstat_yyerror(s) lex_error()

%}

%union {
	int ival;
	unsigned int uval;
	double dval;
	char *sval;
}

%token <ival> INT
%token <uval> UINT
%token <dval> REAL
%token <sval> QSTR  IDENT ID_DATA ID_X ID_VARIOGRAM ID_PREDICTIONS
%token <sval> ID_VARIANCES ID_COVARIANCES ID_OUTPUT ID_MASKS ID_EDGES ID_SET
%token <sval> ID_MERGE ID_AREA ID_BLOCK ID_METHOD ID_BOUNDS ID_MARGINALS

%type <sval> input assign any_id comcol command data_cmd data_decl 
%type <sval> data_cont data_exp data_what data_X data_X_what
%type <sval> vgm_cmd vgm_decl vgm_cont vgm_model vgm_model_type vgm_range
%type <sval> merge_cmd set_cmd set_exp set_lhs mask_cmd mask_cont
%type <sval> edges_cmd edges_cont block_cmd block_cont block_exp block_lhs
%type <sval> output_cmd method_cmd bounds_cmd bounds_exp 
%type <sval> marginals_cmd marginals_cont
%type <dval> val sill_val range_val d_vector d_list d_val

%%
input: { ; } 						    /* allow empty input */
	| input command { reset_parser(); }	/* or a list of commands */
	| d_list { vector_only = 1; }       /* or a list of numbers (stat.c) */
	;

command: ';' { ; }
	| data_cmd    ';' { ; }
	| vgm_cmd     ';' { update_variogram(v); }
	| merge_cmd   ';' { ; }
	| mask_cmd    ';' { ; }
	| edges_cmd  ';'  { ; }
    | block_cmd   ';' { ; }
	| area_cmd    ';' { ; }
	| output_cmd  ';' { ; }
	| method_cmd  ';' { ; }
	| bounds_cmd  ';' { ; }
	| marginals_cmd ';' { ; }
	| set_cmd     ';' { ; }
	;

val: INT { $$ = (double) $1; }
	| REAL ;

assign: '=' { ; }
	| ':' { ; }
	;

comcol: ':' { ; }
	| ',' { ; }
	;

d_vector: '[' d_list ']' { ; }
	;

d_list: d_val 
	| d_list ',' d_val 
	| d_list d_val
	;

d_val: val { 
			if (d == NULL)
				sd_vector = push_d_vector($1, sd_vector);
			else
				d->beta = push_d_vector($1, d->beta);
		}
	;

any_id: IDENT | ID_AREA | ID_BLOCK | ID_BOUNDS | ID_COVARIANCES |
	ID_DATA | ID_MARGINALS | ID_MASKS | ID_METHOD | ID_OUTPUT |
	ID_PREDICTIONS | ID_SET | ID_VARIANCES | ID_VARIOGRAM | ID_X |
	ID_EDGES { ; }
	; /* allows things like  data(data) : ... ; etc. */

data_cmd: data_decl ':' data_cont { verify_data(d); } /* declaration : contents */
	| data_decl { d->dummy = 1; }
	; 

data_decl: ID_DATA '(' any_id ')' {
			id = which_identifier($3);
			dpp = get_gstat_data();
			d = dpp[id];
			d->id = id;
		}
	| ID_DATA '(' ')' {
			d = get_dataval();
			d->id = ID_OF_VALDATA;
		}
	| ID_DATA '(' error { ErrMsg(ER_SYNTAX, "invalid identifier"); }
	;

data_cont: data_exp				/* one data expression */
	| data_cont ',' data_exp	/* a list of data expressions */
	;

data_exp: { ; } /* can be empty */
	| ID_X '=' data_X
	| QSTR { d->fname = $1; }
	| data_what '=' INT {
			switch (expr.what) { 
				case IS_INT: *((int *)expr.ptr) = $3; break;
				case IS_REAL: *((double *)expr.ptr) = (double) $3; break;
				default: lex_error(); YYERROR; break;
			}
			check_assign_expr(&expr);
		}
	| data_what '=' REAL {
			if (expr.what != IS_REAL) {
				lex_error();
				YYERROR;
			}
			*((double *)expr.ptr) = $3;
			check_assign_expr(&expr);
		}
	| data_what '=' QSTR {
			if (expr.what != IS_STRING) {
				lex_error();
				YYERROR;
			}
			*((char **)expr.ptr) = $3;
		}
	| data_what '=' d_vector {
			if (expr.what != IS_D_VECTOR) {
				lex_error();
				YYERROR;
			}
			/*
			*((D_VECTOR **)expr.ptr) = sd_vector; 
			printf("[[ %d ]]\n", sd_vector->size); 
			sd_vector = NULL;
			*/
		}
	| data_what {
			if (expr.what != NO_ARG) {
				lex_error();
				YYERROR;
			}
		}
	;

data_X: data_X_what
	| data_X '&' data_X_what
	;

data_X_what: IDENT {
			for (id = 0; id < N_POLY; id++) {
				if (almost_equals($1, polynomial[id].name)) {
					id += POLY_MIN;
					break; /* i-loop */
				}
			}
			if (id < 0)
				data_add_X(d, id);
			else {
				lex_error();
				YYERROR;
			}
		}
	| INT { push_data_X(d, $1); }
	;

data_what: IDENT {
			if (! is_data_expr(d, &expr, $1)) {
				lex_error();
				YYERROR;
			}
		}
	;

block_cmd: ID_BLOCK {
			bp = get_block_p();
			bp->x = -1.0; /* will be set to grid cell size in predict.c */
		} 
	| ID_BLOCK ':' block_cont
	;

block_cont: block_exp
	| block_cont ',' block_exp
	;

block_exp: block_lhs '=' val {
			*((double *)expr.ptr) = $3;
			check_assign_expr(&expr);
		}
	;

block_lhs: IDENT { if (! is_block_expr(&expr, $1)) { lex_error(); YYERROR; }}
	;

area_cmd: area_decl ':' data_cont { ; }
	;

area_decl: ID_AREA {
			d = create_data_area();
			d->id = ID_OF_AREA;
		}
	| ID_AREA '(' ')' {
			d = create_data_area();
			d->id = ID_OF_AREA;
		}
	;

vgm_cmd: vgm_decl ':' QSTR comcol vgm_cont { v->fname = $3; }
	| vgm_decl ':' QSTR comcol QSTR comcol vgm_cont { v->fname = $3; v->fname2 = $5; }
	| vgm_decl ':' QSTR {v->fname = $3; }
	| vgm_decl ':' QSTR comcol QSTR {v->fname = $3; v->fname2 = $5; }
	| vgm_decl ':' vgm_cont
	| vgm_decl ':' error ';' { 
			/* this will eat the ';' as well, but we're bailing out anyway: */
			YYERROR; 
		}
	;

vgm_decl: ID_VARIOGRAM '(' ')' { 
			/* only allow this when called through read_variogram(): */
			assert(parse_variogram != NULL);
			v = parse_variogram;
			v->n_models = v->n_fit = 0;
		}
	| ID_VARIOGRAM '(' any_id ')' {
			id = which_identifier($3);
			v = get_vgm(LTI(id,id));
			v->id = v->id1 = v->id2 = id;
			v->n_models = v->n_fit = 0;
		}
	| ID_VARIOGRAM '(' any_id ',' any_id ')' {
			id1 = which_identifier($3);
			id2 = which_identifier($5);
			id = LTI(id1,id2);
			v = get_vgm(id);
			v->id = id;
			v->id1 = id1;
			v->id2 = id2;
			v->n_models = v->n_fit = 0;
		}
	;

vgm_cont: vgm_model
	| vgm_cont vgm_model
	;

vgm_model: sill_val vgm_model_type '(' ')' {
			range[0] = 0.0;
			push_to_v(v, $2, $1, range, 1, NULL, fit_sill, fit_range);
		}
	| sill_val vgm_model_type '(' vgm_range ')' {
			push_to_v(v, (const char *) $2, $1, range, nrangepars, anis, fit_sill, fit_range);
		}
	| '+' sill_val vgm_model_type '(' vgm_range ')' {
			push_to_v(v, $3, $2, range, nrangepars, anis, fit_sill, fit_range);
		}
	| '-' sill_val vgm_model_type '(' vgm_range ')' {
			push_to_v(v, $3, -1.0 * $2, range, nrangepars, anis, 
				fit_sill, fit_range);
		}
	;

vgm_model_type:	IDENT { 
			if (which_variogram_model($1) == NOT_SP) {
				lex_error(); YYERROR;
			}
	}
	;

vgm_range: range_val { 
			range[0] = $1; 
			nrangepars = 1;
			anis[0] = -9999.0; 
		}
	| range_val ',' val {
			range[0] = $1;
			range[1] = $3;
			nrangepars = 2;
		}
	| range_val ',' val ',' val {
			range[0] = $1;
			nrangepars = 1;
			anis[0] = $3;
			anis[3] = $5;
			anis[1] = anis[2] = 0.0;
			anis[4] = 1.0;
		}
	| range_val ',' val ',' val ',' val {
			range[0] = $1;
			range[1] = $3;
			nrangepars = 2;
			anis[0] = $5;
			anis[3] = $7;
			anis[1] = anis[2] = 0.0;
			anis[4] = 1.0;
		}
	| range_val ',' val ',' val ',' val ',' val ',' val {
			range[0] = $1;
			nrangepars = 1;
			anis[0] = $3;
			anis[1] = $5;
			anis[2] = $7;
			anis[3] = $9;
			anis[4] = $11;
		}
	| range_val ',' val ',' val ',' val ',' val ',' val ',' val {
			range[0] = $1;
			range[1] = $3;
			nrangepars = 2;
			anis[0] = $5;
			anis[1] = $7;
			anis[2] = $9;
			anis[3] = $11;
			anis[4] = $13;
		}
	;

sill_val: val { fit_sill = 1; }
	| '@' val { fit_sill = 0; $$ = $2; }
	| val '@' { fit_sill = 0; $$ = $1; }
	;

range_val: val { fit_range = 1; }
	| '@' val { fit_range = 0; $$ = $2; }
	;

output_cmd: ID_OUTPUT '=' QSTR { o_filename = $3; }
	| ID_PREDICTIONS '(' any_id ')' ':' QSTR { 
			id = which_identifier($3);
			ofn = (char **) get_outfile_name();
			ofn[2 * id] = $6;
		}
	| ID_PREDICTIONS ':' QSTR { 
			if (get_n_vars() == 0) {
				lex_error();
				ErrMsg(ER_SYNTAX, "define data first");
			}
			ofn = (char **) get_outfile_name();
			ofn[0] = $3;
		}
	| ID_VARIANCES '(' any_id ')' ':' QSTR { 
			id = which_identifier($3);
			ofn = (char **) get_outfile_name();
			ofn[2 * id + 1] = $6;
		}
	| ID_VARIANCES ':' QSTR { 
			if (get_n_vars() == 0) {
				lex_error();
				ErrMsg(ER_SYNTAX, "define data first");
			}
			ofn = (char **) get_outfile_name();
			ofn[1] = $3;
		}
	| ID_COVARIANCES '(' any_id ',' any_id ')' ':'  QSTR { 
			id = get_n_vars();
			id1 = which_identifier($3);
			id2 = which_identifier($5);
			if (id != get_n_vars())	
				ErrMsg(ER_SYNTAX, "define all data(..) before covariances(..,..)");
			ofn = (char **) get_outfile_name();
			id = 2 * id + LTI2(id1, id2);
			ofn[id] = $8;
		}
	;

set_cmd: ID_SET set_exp /* e.g. set nsim = 1 _or_ set out: "out" */
	| set_exp           /* e.g. nsim = 1  or  plot : 'file.plt" */
	; 

set_exp: set_lhs assign INT {
			switch (expr.what) { 
				case IS_INT: *((int *)expr.ptr) = $3; break;
				case IS_REAL: *((double *)expr.ptr) = (double) $3; break;
				default: lex_error(); YYERROR;
			}
			check_assign_expr(&expr);
		}
	| set_lhs assign UINT {
			switch (expr.what) { 
				case IS_UINT: *((unsigned int *)expr.ptr) = $3; break;
				default: lex_error(); YYERROR;
			}
			check_assign_expr(&expr);
		}
	| set_lhs assign REAL {
			switch (expr.what) {
				case IS_REAL: *((double *)expr.ptr) = $3; break;
				default: lex_error(); YYERROR; break;
			}
			check_assign_expr(&expr);
		}
	| set_lhs assign QSTR {
			if (expr.what != IS_STRING) {
				lex_error();
				YYERROR;
			}
			*((char **) expr.ptr) = $3;
		}
	;

set_lhs: IDENT { if (! is_set_expr(&expr, $1)) { lex_error(); YYERROR; }}
	;

method_cmd: ID_METHOD ':' IDENT {
			for (id = 1; methods[id].name != NULL; id++) {
				if (almost_equals($3, methods[id].name)) {
					set_method(methods[id].m);
					break; /* id-loop */
				}
			}
			if (methods[id].m == NSP) {
				lex_error();
				YYERROR;
			}
		}
	;

mask_cmd: ID_MASKS ':' mask_cont
	;

mask_cont: QSTR { push_mask_name($1); }
	| mask_cont ',' QSTR { push_mask_name($3); }
	;

edges_cmd: ID_EDGES ':' edges_cont
	;

edges_cont: QSTR { push_edges_name($1); }
	| edges_cont ',' QSTR { push_edges_name($3); }
	;

merge_cmd: ID_MERGE any_id IDENT any_id {
			if (!almost_equals($3, "w$ith"))
				lex_error();
			id1 = which_identifier($2);
			id2 = which_identifier($4);
			dpp = get_gstat_data();
			if (dpp[id1]->id != id1 || dpp[id2]->id != id2) {
				lex_error();
				ErrMsg(ER_IMPOSVAL, "define data before attempting to merge");
			}
			col1 = col2 = 0;
			if (id1 < id2) { /* swap id's */
				id = id1; id1 = id2; id2 = id;
			}
			if (push_to_merge_table(dpp[id1], id2, col1, col2)) {
				lex_error();
				ErrMsg(ER_IMPOSVAL, "attempt to merge failed");
			}
		}
	| ID_MERGE any_id '(' INT ')' IDENT any_id '(' INT ')' {
			if (!almost_equals($6, "w$ith"))
				lex_error();
			id1 = which_identifier($2);
			id2 = which_identifier($7);
			dpp = get_gstat_data();
			if (dpp[id1]->id != id1 || dpp[id2]->id != id2) {
				lex_error();
				ErrMsg(ER_IMPOSVAL, "define data before attempting to merge");
			}
			col1 = $4;
			col2 = $9;
			if (id1 < id2) { /* swap id and col */
				id = id1; id1 = id2; id2 = id;
				id = col1; col1 = col2; col2 = id;
			}
			if (push_to_merge_table(dpp[id1], id2, col1, col2)) {
				lex_error();
				ErrMsg(ER_IMPOSVAL, "attempt to merge failed");
			}
		}
	| ID_MERGE INT '(' INT ')' IDENT INT '(' INT ')' {
			if (!almost_equals($6, "w$ith"))
				lex_error();
			id1 = $2;
			id2 = $7;
			if (id1 >= get_n_vars() || id2 >= get_n_vars() || id1 < 0 || id2 < 0) {
				lex_error();
				ErrMsg(ER_IMPOSVAL, "id values out of range");
			}
			col1 = $4;
			col2 = $9;
			if (id1 < id2) { /* swap id and col */
				id = id1; id1 = id2; id2 = id;
				id = col1; col1 = col2; col2 = id;
			}
			dpp = get_gstat_data();
			if (push_to_merge_table(dpp[id1], id2, col1, col2)) {
				lex_error();
				ErrMsg(ER_IMPOSVAL, "attempt to merge failed");
			}
		}
	;

bounds_cmd: ID_BOUNDS ':' QSTR { boundary_file = $3; }
	| ID_BOUNDS ':' bounds_exp
	;

bounds_exp: val { push_bound($1); }
	| bounds_exp val { push_bound($2); }
	| bounds_exp ',' val { push_bound($3); }
	;

marginals_cmd : ID_MARGINALS ':' marginals_cont
	;

marginals_cont: val             { push_marginal(NULL, $1); }
	| marginals_cont ',' val    { push_marginal(NULL, $3); }
	| QSTR                      { push_marginal($1, -1.0); }
	| marginals_cont ',' QSTR   { push_marginal($3, -1.0); }
	;

%%

static int is_data_expr(DATA *d, GSTAT_EXPR *expr, const char *fld) {
#define TABLE_SIZE 32
	GSTAT_EXPR data_options[TABLE_SIZE];
	int i = 0;
#define FILL_TABLE(n, p, w, l) \
 data_options[i].name = n; data_options[i].ptr = p; \
 data_options[i].what = w; data_options[i].limit = l; i++;

/* set up table: */
	FILL_TABLE("x",          &(d->colnx),        IS_INT, GTZERO  )
	FILL_TABLE("y",          &(d->colny),        IS_INT, GTZERO  )
	FILL_TABLE("z",          &(d->colnz),        IS_INT, GTZERO  )
	FILL_TABLE("v",          &(d->colnvalue),    IS_INT, GTZERO  )
	FILL_TABLE("V",          &(d->colnvariance), IS_INT, GTZERO  )
	FILL_TABLE("d",          &(d->polynomial_degree), IS_INT, GEZERO  )
	FILL_TABLE("max",        &(d->sel_max),      IS_INT, GEZERO  )
	FILL_TABLE("omax",       &(d->oct_max),      IS_INT, GEZERO  )
	FILL_TABLE("min",        &(d->sel_min),      IS_INT, GEZERO  )
	FILL_TABLE("n$max",      &(d->init_max),     IS_INT, GEZERO  )
	FILL_TABLE("togrid",     &(d->togrid),       IS_INT,  GEZERO )
	FILL_TABLE("I",          &(d->Icutoff),      IS_REAL, NOLIMIT )
	FILL_TABLE("mv",         &(d->mv),           IS_REAL, NOLIMIT )
	FILL_TABLE("rad$ius",    &(d->sel_rad),      IS_REAL, GTZERO  )
	FILL_TABLE("dX",         &(d->dX),           IS_REAL, GEZERO  )
	FILL_TABLE("b$eta",      &(d->beta),         IS_D_VECTOR, NOLIMIT )
	FILL_TABLE("stan$dard",  &(d->standard),     NO_ARG, NOLIMIT )
	FILL_TABLE("log",        &(d->log),          NO_ARG, NOLIMIT )
	FILL_TABLE("av$erage",   &(d->average),      IS_INT, GEZERO )
	FILL_TABLE("re$gion",    &(d->region),       NO_ARG, NOLIMIT )
	FILL_TABLE("du$mmy",     &(d->dummy),        NO_ARG, NOLIMIT )
	FILL_TABLE("res$idual",  &(d->calc_residuals), NO_ARG, NOLIMIT )
	FILL_TABLE("vdist",      &(d->vdist),        NO_ARG, NOLIMIT )
	FILL_TABLE("force",      &(d->force),        NO_ARG, NOLIMIT )
	FILL_TABLE("Cat$egory",   &(d->Category),    IS_STRING, NOLIMIT )
	FILL_TABLE("ID",         &(d->coln_id),      IS_INT, GTZERO  )
	FILL_TABLE("VarF$unction", &(d->var_fn_str), IS_STRING, NOLIMIT  )
	FILL_TABLE("nscore",     &(d->nscore_table), IS_STRING, NOLIMIT )
	FILL_TABLE("every",      &(d->every),        IS_INT, GTZERO )
	FILL_TABLE("offset",     &(d->offset),       IS_INT, NOLIMIT )
	FILL_TABLE("prob",       &(d->prob),         IS_REAL, GTZERO )
	FILL_TABLE(NULL, NULL, IS_INT, NOLIMIT )

/* check TABLE_SIZE was set correctly... */
	assert(i == TABLE_SIZE); 

	expr->ptr = NULL;
	expr->what = UNKNOWN;
	expr->limit = NOLIMIT;

	for (i = 0; data_options[i].name != NULL; i++) {
		if (almost_equals(fld, data_options[i].name)) {
			expr->name = fld;
			expr->ptr = data_options[i].ptr;
			expr->what = data_options[i].what;
			expr->limit = data_options[i].limit;
			if (expr->what == NO_ARG)
				*((int *) expr->ptr) = 1;
			return 1;
		}
	}

	/* non-standard cases not in data_options[] table: */
	if (almost_equals(fld, "s$tratum")) {
		if (d->id != ID_OF_VALDATA)
			return 0;
		expr->ptr = &(d->colns); 
		expr->what = IS_INT; 
		expr->limit = GTZERO;
		d->what_is_u = U_ISSTRATUM;
	} else if (almost_equals(fld, "av$erage")) {
		d->average = 1; expr->what = NO_ARG;
	} else if (almost_equals(fld, "noav$erage")) {
		d->average = 0; expr->what = NO_ARG;
	} else if (almost_equals(fld, "nores$idual")) {
		d->calc_residuals = 0; expr->what = NO_ARG;
	} else if (almost_equals(fld, "square")) {
		d->square = 1; expr->what = NO_ARG;
	} else if (almost_equals(fld, "c")) {
		pr_warning("use `v' instead of `c' in data definition");
	} else if (almost_equals(fld, "sk_mean")) { /* move it to beta: */
		d->beta = NULL;
		d->beta = push_d_vector(-9999.0, d->beta);
		expr->ptr = &(d->beta->val[0]);
		expr->what = IS_REAL; 
		expr->limit = NOLIMIT;
	} 

	return (expr->what != UNKNOWN);
}

static int is_set_expr(GSTAT_EXPR *expr, const char *name) {
/*
 * parse sequences like `set zmap = 50.0;' or `set zmap = 50, idp = 2.5;'
 * (int, float or string)
 */
	int i;

	const GSTAT_EXPR set_options[] = {
	{ "cn$_max",        &gl_cn_max,       IS_REAL, GTZERO  },
	{ "co$incide",      &gl_coincide,     IS_INT,  GEZERO  },
	{ "Cr$essie",       &gl_cressie,      IS_INT,  GEZERO  },
	{ "a$lpha",         &gl_alpha,        IS_REAL, GEZERO  },
	{ "b$eta",          &gl_beta,         IS_REAL, GEZERO  },
	{ "c$utoff",        &gl_cutoff,       IS_REAL, GTZERO  },
	{ "de$bug",         &debug_level,     IS_INT,  GEZERO  },
	{ "display",        &gl_display,      IS_STRING, NOLIMIT },
	{ "do$ts",          &gl_dots,         IS_INT,  GEZERO  },
	{ "fit",            &gl_fit,          IS_INT,  GEZERO  },
	{ "fit_l$imit",     &gl_fit_limit,    IS_REAL, GTZERO  },
	{ "fo$rmat",        &gl_format,       IS_STRING, NOLIMIT },
	{ "fr$action",      &gl_fraction,     IS_REAL, GTZERO  },
	{ "gcv",            &gl_gcv,          IS_REAL, GTZERO  },
	{ "gls$_residuals", &gl_gls_residuals, IS_INT, GEZERO  },
	{ "gnuplot",        &gl_gnuplot,      IS_STRING, NOLIMIT },
	{ "gnuplot35",      &gl_gnuplot35,    IS_STRING, NOLIMIT },
	{ "gpt$erm",        &gl_gpterm,       IS_STRING, NOLIMIT  },
	{ "id$p",           &gl_idp,          IS_REAL, GEZERO  },
	{ "in$tervals",     &gl_n_intervals,  IS_INT,  GTZERO  },
	{ "it$er",          &gl_iter,         IS_INT,  GEZERO  },
	{ "j$graph",        &gl_jgraph,       IS_INT,  GEZERO  },
	{ "lhs",            &gl_lhs,          IS_INT,  GEZERO  },
	{ "log$file",       &logfile_name,    IS_STRING, NOLIMIT },
	{ "longlat",        &gl_longlat,      IS_INT, GEZERO },
	{ "sim_beta",  		&gl_sim_beta,     IS_INT, GEZERO },
	{ "mv$string",		&gl_mv_string,    IS_STRING, NOLIMIT },
	{ "n_uk",           &gl_n_uk,         IS_INT,  GEZERO  },
	{ "numbers",        &gl_numbers,      IS_INT,  GEZERO  },
	{ "nb$lockdiscr",   &gl_nblockdiscr,  IS_INT,  GTZERO  },
	{ "no$check",       &gl_nocheck,      IS_INT,  GEZERO  },
	{ "ns$im",          &gl_nsim,         IS_INT,  GTZERO  },
	{ "o$utputfile",    &o_filename,      IS_STRING, NOLIMIT },
	{ "or$der",         &gl_order,        IS_INT,  GEZERO },
	{ "pag$er",         &gl_pager,        IS_STRING, NOLIMIT },
	{ "pl$otfile",      &gl_plotfile,     IS_STRING, NOLIMIT },
	{ "q$uantile",      &gl_quantile,     IS_REAL, GEZERO  },
	{ "rowwise",        &gl_rowwise,      IS_INT,  GEZERO  },
	{ "rp",             &gl_rp,           IS_INT,  GEZERO  },
	{ "sec$ure",        &gl_secure,       IS_INT,  GTZERO  },
	{ "see$d",          &gl_seed,         IS_INT,  GTZERO  },
	{ "useed",          &gl_seed,         IS_UINT,  GEZERO  },
	{ "spa$rse",        &gl_sparse,       IS_INT,  GEZERO  },
	{ "spi$ral",        &gl_spiral,       IS_INT,  GEZERO  },
	{ "spl$it",         &gl_split,        IS_INT,  GTZERO  },
	{ "sy$mmetric",     &gl_sym_ev,       IS_INT,  GEZERO  },
	{ "tol_h$or",       &gl_tol_hor,      IS_REAL, GEZERO  },
	{ "tol_v$er",       &gl_tol_ver,      IS_REAL, GEZERO  },
	{ "v$erbose",       &debug_level,     IS_INT,  GEZERO  },
	{ "w$idth",         &gl_iwidth,       IS_REAL, GEZERO  },
	{ "x$valid",        &gl_xvalid,       IS_INT,  GEZERO  },
	{ "zero_di$st",     &gl_zero_est,     IS_INT,  GEZERO  },
	{ "zero",           &gl_zero,         IS_REAL, GEZERO  },
	{ "zm$ap",          &gl_zmap,         IS_REAL, NOLIMIT },
	{ "plotw$eights",   &gl_plotweights,  IS_INT, GEZERO   },
	{ NULL, NULL, 0, 0 }
	};

	for (i = 0; set_options[i].name; i++)
		if (almost_equals(name, set_options[i].name))
			break; /* break out i-loop */
	if (set_options[i].name == NULL)
		return 0;

	if (almost_equals((const char *)name,"nb$lockdiscr"))
		gl_gauss = 0; /* side effect */

	expr->name = name;
	expr->ptr = set_options[i].ptr;
	expr->what = set_options[i].what;
	expr->limit = set_options[i].limit;

	return 1;
}

static void check_assign_expr(GSTAT_EXPR *expr) {
/* for INT and REAL expressions, check range */
	double val;

	switch(expr->what) {
		case IS_INT: 
			val = (double) (*((int *)(expr->ptr)));
			break;
		case IS_REAL: 
			val = (*((double *)(expr->ptr)));
			break;
		default:
			return;
	}
	if (expr->limit == GEZERO && val < 0.0) {
		lex_error();
		pr_warning("value should be non-negative");
		ErrMsg(ER_IMPOSVAL, expr->name);
	}
	if (expr->limit == GTZERO && val <= 0.0) {
		lex_error();
		pr_warning("value should be positive");
		ErrMsg(ER_IMPOSVAL, expr->name);
	}
}

static void push_data_X(DATA *d, int id) {
	if (id == -1) { /* remove default intercept */
		if (d->n_X > 1) {
			lex_error();
			ErrMsg(ER_SYNTAX, "-1 only as first argument following X="); 
		}
		d->n_X = 0;
	} else if (id == 0) {
		lex_error();
		ErrMsg(ER_SYNTAX, "intercept is default"); 
	} else /* id > 0 */
		data_add_X(d, id);
}

static int is_block_expr(GSTAT_EXPR *expr, const char *s) {
	DPOINT *bp;

	bp = get_block_p();
	expr->name = s;
	expr->limit = GEZERO;
	expr->what = IS_REAL;
	if (almost_equals(s, "dx"))
		expr->ptr = &(bp->x);
	else if (almost_equals(s, "dy"))
		expr->ptr = &(bp->y);
	else if (almost_equals(s, "dz"))
		expr->ptr = &(bp->z);
	else
		return 0;
	return 1;
}

static void push_marginal(char *name, double value) {
	static int names = -1;

	if (names == -1)
		names = (name != NULL);

	if (name) {
		if (!names) {
			lex_error();
			ErrMsg(ER_SYNTAX, "only real values allowed"); 
		}
		gl_marginal_names = (char **) erealloc(gl_marginal_names,
			++gl_n_marginals * sizeof(char *));
		gl_marginal_names[gl_n_marginals - 1] = name;
	} else {
		if (names) {
			lex_error();
			ErrMsg(ER_SYNTAX, "only quoted strings allowed"); 
		}
		gl_marginal_values = (double *) erealloc (gl_marginal_values,
			++gl_n_marginals * sizeof(double));
		gl_marginal_values[gl_n_marginals - 1] = value;
	}
	return;
}

static void reset_parser(void) {
/* savety first: reset all static globals (should be unnessesary) */
	v = NULL;
	d = NULL;
	bp = NULL;
	ofn = NULL;
	expr.ptr = NULL;
	expr.what =  UNKNOWN;
	expr.limit =  NOLIMIT;
	id = id1 = id2 = col1 = col2 = -1;
}

int parse_cmd(const char *cmd, const char *fname) {
	set_lex_source(cmd, fname);
	reset_parser();
	return yyparse();
}

#ifndef USING_R
int parse_file(const char *fname) {
/* 
 * parse commands in file fname
 */
	int stdin_isatty = 1;
	char *cp;

	if (fname == NULL || strcmp(fname, "-") == 0) {
#ifdef HAVE_UNISTD_H
		stdin_isatty = isatty(fileno(stdin));
#endif
		if (stdin_isatty)
			cp = string_prompt("gstat> ");
		else
			cp = string_file(NULL);
	} else /* read from file */
		cp = string_file(fname);

	if (parse_cmd(cp, fname))
		ErrMsg(ER_SYNTAX, fname);
	efree(cp);

	if (boundary_file != NULL) {
		cp = string_file(boundary_file);
		if (parse_cmd(cp, boundary_file))
			ErrMsg(ER_SYNTAX, boundary_file);
		efree(cp);
	}

	if (vector_only && !allow_vector_only)
		ErrMsg(ER_SYNTAX, fname);

	return 0;
}

void parse_gstatrc(void) {
	char *fname = NULL, *cp;

	if ((fname = getenv(GSTATRC)) != NULL) {
		if (! file_exists(fname)) {
			message("environment variable %s:\n", GSTATRC);
			ErrMsg(ER_READ, fname);
		}
		parse_file(fname);
	} else if ((cp = getenv("HOME")) != NULL) {
		fname = (char *) emalloc(strlen(cp) + strlen(HOMERCFILE) + 2);
		sprintf(fname, "%s/%s", cp, HOMERCFILE);
		if (file_exists(fname))
			parse_file(fname);
		efree(fname);
	}
	return;
}

int read_variogram(VARIOGRAM *v, const char *source) {
	char *cp;
	int rval;

	parse_variogram = v;
	cp = (char *) emalloc((strlen(source) + 20) * sizeof(char));
	sprintf(cp, "variogram(): %s;", source);
	rval = parse_cmd(cp, NULL);
	parse_variogram = NULL; /* for safety */
	efree(cp);
	return rval;
}

int read_vector(D_VECTOR *d, char *fname) {
	int rval;

	assert(d != NULL);
	sd_vector = d;

	allow_vector_only = 1;

	rval = parse_file(fname);

	if (! vector_only)  {
		message("stat: only numeric input allowed -- \n");
		ErrMsg(ER_IMPOSVAL, fname);
	}

	return rval;
}
#endif

static void verify_data(DATA *d) { /* declaration : contents */

	if (d->var_fn_str != NULL) {
		if (almost_equals(d->var_fn_str, "mu"))
			d->variance_fn = v_mu;
		else if (almost_equals(d->var_fn_str, "mu^2"))
			d->variance_fn = v_bin;
		else if (almost_equals(d->var_fn_str, "mu^3"))
			d->variance_fn = v_bin;
		else if (almost_equals(d->var_fn_str, "mu(1-mu)"))
			d->variance_fn = v_bin;
		else if (almost_equals(d->var_fn_str, "identity"))
			d->variance_fn = v_identity;
		else {
			lex_error();
			message("variance function %s not supported:\n", d->var_fn_str);
			ErrMsg(ER_SYNTAX, d->var_fn_str);
		}
	}
}
