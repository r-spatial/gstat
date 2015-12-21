
/* A Bison parser, made by GNU Bison 2.4.1.  */

/* Skeleton implementation for Bison's Yacc-like parsers in C
   
      Copyright (C) 1984, 1989, 1990, 2000, 2001, 2002, 2003, 2004, 2005, 2006
   Free Software Foundation, Inc.
   
   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.
   
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
   
   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.  */

/* As a special exception, you may create a larger work that contains
   part or all of the Bison parser skeleton and distribute that work
   under terms of your choice, so long as that work isn't itself a
   parser generator using the skeleton or a modified version thereof
   as a parser skeleton.  Alternatively, if you modify or redistribute
   the parser skeleton itself, you may (at your option) remove this
   special exception, which will cause the skeleton and the resulting
   Bison output files to be licensed under the GNU General Public
   License without this special exception.
   
   This special exception was added by the Free Software Foundation in
   version 2.2 of Bison.  */

/* C LALR(1) parser skeleton written by Richard Stallman, by
   simplifying the original so-called "semantic" parser.  */

/* All symbols defined below should begin with yy or YY, to avoid
   infringing on user name space.  This should be done even for local
   variables, as they might otherwise be expanded by user macros.
   There are some unavoidable exceptions within include files to
   define necessary library symbols; they are noted "INFRINGES ON
   USER NAME SPACE" below.  */

/* Identify Bison output.  */
#define YYBISON 1

/* Bison version.  */
#define YYBISON_VERSION "2.4.1"

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 0

/* Push parsers.  */
#define YYPUSH 0

/* Pull parsers.  */
#define YYPULL 1

/* Using locations.  */
#define YYLSP_NEEDED 0

/* Substitute the variable and function names.  */
#define yyparse         gstat_yyparse
#define yylex           gstat_yylex
#define yyerror         gstat_yyerror
#define yylval          gstat_yylval
#define yychar          gstat_yychar
#define yydebug         gstat_yydebug
#define yynerrs         gstat_yynerrs


/* Copy the first part of user declarations.  */

/* Line 189 of yacc.c  */
#line 1 "parse.y"

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



/* Line 189 of yacc.c  */
#line 199 "y.tab.c"

/* Enabling traces.  */
#ifndef YYDEBUG
# define YYDEBUG 0
#endif

/* Enabling verbose error messages.  */
#ifdef YYERROR_VERBOSE
# undef YYERROR_VERBOSE
# define YYERROR_VERBOSE 1
#else
# define YYERROR_VERBOSE 0
#endif

/* Enabling the token table.  */
#ifndef YYTOKEN_TABLE
# define YYTOKEN_TABLE 0
#endif


/* Tokens.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
   /* Put the tokens into the symbol table, so that GDB and other debuggers
      know about them.  */
   enum yytokentype {
     INT = 258,
     UINT = 259,
     REAL = 260,
     QSTR = 261,
     IDENT = 262,
     ID_DATA = 263,
     ID_X = 264,
     ID_VARIOGRAM = 265,
     ID_PREDICTIONS = 266,
     ID_VARIANCES = 267,
     ID_COVARIANCES = 268,
     ID_OUTPUT = 269,
     ID_MASKS = 270,
     ID_EDGES = 271,
     ID_SET = 272,
     ID_MERGE = 273,
     ID_AREA = 274,
     ID_BLOCK = 275,
     ID_METHOD = 276,
     ID_BOUNDS = 277,
     ID_MARGINALS = 278
   };
#endif
/* Tokens.  */
#define INT 258
#define UINT 259
#define REAL 260
#define QSTR 261
#define IDENT 262
#define ID_DATA 263
#define ID_X 264
#define ID_VARIOGRAM 265
#define ID_PREDICTIONS 266
#define ID_VARIANCES 267
#define ID_COVARIANCES 268
#define ID_OUTPUT 269
#define ID_MASKS 270
#define ID_EDGES 271
#define ID_SET 272
#define ID_MERGE 273
#define ID_AREA 274
#define ID_BLOCK 275
#define ID_METHOD 276
#define ID_BOUNDS 277
#define ID_MARGINALS 278




#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef union YYSTYPE
{

/* Line 214 of yacc.c  */
#line 118 "parse.y"

	int ival;
	unsigned int uval;
	double dval;
	char *sval;



/* Line 214 of yacc.c  */
#line 290 "y.tab.c"
} YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
#endif


/* Copy the second part of user declarations.  */


/* Line 264 of yacc.c  */
#line 302 "y.tab.c"

#ifdef short
# undef short
#endif

#ifdef YYTYPE_UINT8
typedef YYTYPE_UINT8 yytype_uint8;
#else
typedef unsigned char yytype_uint8;
#endif

#ifdef YYTYPE_INT8
typedef YYTYPE_INT8 yytype_int8;
#elif (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
typedef signed char yytype_int8;
#else
typedef short int yytype_int8;
#endif

#ifdef YYTYPE_UINT16
typedef YYTYPE_UINT16 yytype_uint16;
#else
typedef unsigned short int yytype_uint16;
#endif

#ifdef YYTYPE_INT16
typedef YYTYPE_INT16 yytype_int16;
#else
typedef short int yytype_int16;
#endif

#ifndef YYSIZE_T
# ifdef __SIZE_TYPE__
#  define YYSIZE_T __SIZE_TYPE__
# elif defined size_t
#  define YYSIZE_T size_t
# elif ! defined YYSIZE_T && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
#  include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  define YYSIZE_T size_t
# else
#  define YYSIZE_T unsigned int
# endif
#endif

#define YYSIZE_MAXIMUM ((YYSIZE_T) -1)

#ifndef YY_
# if YYENABLE_NLS
#  if ENABLE_NLS
#   include <libintl.h> /* INFRINGES ON USER NAME SPACE */
#   define YY_(msgid) dgettext ("bison-runtime", msgid)
#  endif
# endif
# ifndef YY_
#  define YY_(msgid) msgid
# endif
#endif

/* Suppress unused-variable warnings by "using" E.  */
#if ! defined lint || defined __GNUC__
# define YYUSE(e) ((void) (e))
#else
# define YYUSE(e) /* empty */
#endif

/* Identity function, used to suppress warnings about constant conditions.  */
#ifndef lint
# define YYID(n) (n)
#else
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static int
YYID (int yyi)
#else
static int
YYID (yyi)
    int yyi;
#endif
{
  return yyi;
}
#endif

#if ! defined yyoverflow || YYERROR_VERBOSE

/* The parser invokes alloca or malloc; define the necessary symbols.  */

# ifdef YYSTACK_USE_ALLOCA
#  if YYSTACK_USE_ALLOCA
#   ifdef __GNUC__
#    define YYSTACK_ALLOC __builtin_alloca
#   elif defined __BUILTIN_VA_ARG_INCR
#    include <alloca.h> /* INFRINGES ON USER NAME SPACE */
#   elif defined _AIX
#    define YYSTACK_ALLOC __alloca
#   elif defined _MSC_VER
#    include <malloc.h> /* INFRINGES ON USER NAME SPACE */
#    define alloca _alloca
#   else
#    define YYSTACK_ALLOC alloca
#    if ! defined _ALLOCA_H && ! defined _STDLIB_H && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
#     include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#     ifndef _STDLIB_H
#      define _STDLIB_H 1
#     endif
#    endif
#   endif
#  endif
# endif

# ifdef YYSTACK_ALLOC
   /* Pacify GCC's `empty if-body' warning.  */
#  define YYSTACK_FREE(Ptr) do { /* empty */; } while (YYID (0))
#  ifndef YYSTACK_ALLOC_MAXIMUM
    /* The OS might guarantee only one guard page at the bottom of the stack,
       and a page size can be as small as 4096 bytes.  So we cannot safely
       invoke alloca (N) if N exceeds 4096.  Use a slightly smaller number
       to allow for a few compiler-allocated temporary stack slots.  */
#   define YYSTACK_ALLOC_MAXIMUM 4032 /* reasonable circa 2006 */
#  endif
# else
#  define YYSTACK_ALLOC YYMALLOC
#  define YYSTACK_FREE YYFREE
#  ifndef YYSTACK_ALLOC_MAXIMUM
#   define YYSTACK_ALLOC_MAXIMUM YYSIZE_MAXIMUM
#  endif
#  if (defined __cplusplus && ! defined _STDLIB_H \
       && ! ((defined YYMALLOC || defined malloc) \
	     && (defined YYFREE || defined free)))
#   include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#   ifndef _STDLIB_H
#    define _STDLIB_H 1
#   endif
#  endif
#  ifndef YYMALLOC
#   define YYMALLOC malloc
#   if ! defined malloc && ! defined _STDLIB_H && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
void *malloc (YYSIZE_T); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
#  ifndef YYFREE
#   define YYFREE free
#   if ! defined free && ! defined _STDLIB_H && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
void free (void *); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
# endif
#endif /* ! defined yyoverflow || YYERROR_VERBOSE */


#if (! defined yyoverflow \
     && (! defined __cplusplus \
	 || (defined YYSTYPE_IS_TRIVIAL && YYSTYPE_IS_TRIVIAL)))

/* A type that is properly aligned for any stack member.  */
union yyalloc
{
  yytype_int16 yyss_alloc;
  YYSTYPE yyvs_alloc;
};

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAXIMUM (sizeof (union yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# define YYSTACK_BYTES(N) \
     ((N) * (sizeof (yytype_int16) + sizeof (YYSTYPE)) \
      + YYSTACK_GAP_MAXIMUM)

/* Copy COUNT objects from FROM to TO.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if defined __GNUC__ && 1 < __GNUC__
#   define YYCOPY(To, From, Count) \
      __builtin_memcpy (To, From, (Count) * sizeof (*(From)))
#  else
#   define YYCOPY(To, From, Count)		\
      do					\
	{					\
	  YYSIZE_T yyi;				\
	  for (yyi = 0; yyi < (Count); yyi++)	\
	    (To)[yyi] = (From)[yyi];		\
	}					\
      while (YYID (0))
#  endif
# endif

/* Relocate STACK from its old location to the new one.  The
   local variables YYSIZE and YYSTACKSIZE give the old and new number of
   elements in the stack, and YYPTR gives the new location of the
   stack.  Advance YYPTR to a properly aligned location for the next
   stack.  */
# define YYSTACK_RELOCATE(Stack_alloc, Stack)				\
    do									\
      {									\
	YYSIZE_T yynewbytes;						\
	YYCOPY (&yyptr->Stack_alloc, Stack, yysize);			\
	Stack = &yyptr->Stack_alloc;					\
	yynewbytes = yystacksize * sizeof (*Stack) + YYSTACK_GAP_MAXIMUM; \
	yyptr += yynewbytes / sizeof (*yyptr);				\
      }									\
    while (YYID (0))

#endif

/* YYFINAL -- State number of the termination state.  */
#define YYFINAL  7
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   298

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  36
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  45
/* YYNRULES -- Number of rules.  */
#define YYNRULES  133
/* YYNRULES -- Number of states.  */
#define YYNSTATES  253

/* YYTRANSLATE(YYLEX) -- Bison symbol number corresponding to YYLEX.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   278

#define YYTRANSLATE(YYX)						\
  ((unsigned int) (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[YYLEX] -- Bison symbol number corresponding to YYLEX.  */
static const yytype_uint8 yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,    32,     2,
      30,    31,     2,    33,    27,    34,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,    26,    24,
       2,    25,     2,     2,    35,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,    28,     2,    29,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     1,     2,     3,     4,
       5,     6,     7,     8,     9,    10,    11,    12,    13,    14,
      15,    16,    17,    18,    19,    20,    21,    22,    23
};

#if YYDEBUG
/* YYPRHS[YYN] -- Index of the first RHS symbol of rule number YYN in
   YYRHS.  */
static const yytype_uint16 yyprhs[] =
{
       0,     0,     3,     4,     7,     9,    11,    14,    17,    20,
      23,    26,    29,    32,    35,    38,    41,    44,    47,    49,
      51,    53,    55,    57,    59,    63,    65,    69,    72,    74,
      76,    78,    80,    82,    84,    86,    88,    90,    92,    94,
      96,    98,   100,   102,   104,   106,   110,   112,   117,   121,
     125,   127,   131,   132,   136,   138,   142,   146,   150,   154,
     156,   158,   162,   164,   166,   168,   170,   174,   176,   180,
     184,   186,   190,   192,   196,   202,   210,   214,   220,   224,
     229,   233,   238,   245,   247,   250,   255,   261,   268,   275,
     277,   279,   283,   289,   297,   309,   323,   325,   328,   331,
     333,   336,   340,   347,   351,   358,   362,   371,   374,   376,
     380,   384,   388,   392,   394,   398,   402,   404,   408,   412,
     414,   418,   423,   434,   445,   449,   453,   455,   458,   462,
     466,   468,   472,   474
};

/* YYRHS -- A `-1'-separated list of the rules' RHS.  */
static const yytype_int8 yyrhs[] =
{
      37,     0,    -1,    -1,    37,    38,    -1,    43,    -1,    24,
      -1,    46,    24,    -1,    59,    24,    -1,    76,    24,    -1,
      72,    24,    -1,    74,    24,    -1,    53,    24,    -1,    57,
      24,    -1,    67,    24,    -1,    71,    24,    -1,    77,    24,
      -1,    79,    24,    -1,    68,    24,    -1,     3,    -1,     5,
      -1,    25,    -1,    26,    -1,    26,    -1,    27,    -1,    28,
      43,    29,    -1,    44,    -1,    43,    27,    44,    -1,    43,
      44,    -1,    39,    -1,     7,    -1,    19,    -1,    20,    -1,
      22,    -1,    13,    -1,     8,    -1,    23,    -1,    15,    -1,
      21,    -1,    14,    -1,    11,    -1,    17,    -1,    12,    -1,
      10,    -1,     9,    -1,    16,    -1,    47,    26,    48,    -1,
      47,    -1,     8,    30,    45,    31,    -1,     8,    30,    31,
      -1,     8,    30,     1,    -1,    49,    -1,    48,    27,    49,
      -1,    -1,     9,    25,    50,    -1,     6,    -1,    52,    25,
       3,    -1,    52,    25,     5,    -1,    52,    25,     6,    -1,
      52,    25,    42,    -1,    52,    -1,    51,    -1,    50,    32,
      51,    -1,     7,    -1,     3,    -1,     7,    -1,    20,    -1,
      20,    26,    54,    -1,    55,    -1,    54,    27,    55,    -1,
      56,    25,    39,    -1,     7,    -1,    58,    26,    48,    -1,
      19,    -1,    19,    30,    31,    -1,    60,    26,     6,    41,
      61,    -1,    60,    26,     6,    41,     6,    41,    61,    -1,
      60,    26,     6,    -1,    60,    26,     6,    41,     6,    -1,
      60,    26,    61,    -1,    60,    26,     1,    24,    -1,    10,
      30,    31,    -1,    10,    30,    45,    31,    -1,    10,    30,
      45,    27,    45,    31,    -1,    62,    -1,    61,    62,    -1,
      65,    63,    30,    31,    -1,    65,    63,    30,    64,    31,
      -1,    33,    65,    63,    30,    64,    31,    -1,    34,    65,
      63,    30,    64,    31,    -1,     7,    -1,    66,    -1,    66,
      27,    39,    -1,    66,    27,    39,    27,    39,    -1,    66,
      27,    39,    27,    39,    27,    39,    -1,    66,    27,    39,
      27,    39,    27,    39,    27,    39,    27,    39,    -1,    66,
      27,    39,    27,    39,    27,    39,    27,    39,    27,    39,
      27,    39,    -1,    39,    -1,    35,    39,    -1,    39,    35,
      -1,    39,    -1,    35,    39,    -1,    14,    25,     6,    -1,
      11,    30,    45,    31,    26,     6,    -1,    11,    26,     6,
      -1,    12,    30,    45,    31,    26,     6,    -1,    12,    26,
       6,    -1,    13,    30,    45,    27,    45,    31,    26,     6,
      -1,    17,    69,    -1,    69,    -1,    70,    40,     3,    -1,
      70,    40,     4,    -1,    70,    40,     5,    -1,    70,    40,
       6,    -1,     7,    -1,    21,    26,     7,    -1,    15,    26,
      73,    -1,     6,    -1,    73,    27,     6,    -1,    16,    26,
      75,    -1,     6,    -1,    75,    27,     6,    -1,    18,    45,
       7,    45,    -1,    18,    45,    30,     3,    31,     7,    45,
      30,     3,    31,    -1,    18,     3,    30,     3,    31,     7,
       3,    30,     3,    31,    -1,    22,    26,     6,    -1,    22,
      26,    78,    -1,    39,    -1,    78,    39,    -1,    78,    27,
      39,    -1,    23,    26,    80,    -1,    39,    -1,    80,    27,
      39,    -1,     6,    -1,    80,    27,     6,    -1
};

/* YYRLINE[YYN] -- source line where rule number YYN was defined.  */
static const yytype_uint16 yyrline[] =
{
       0,   142,   142,   143,   144,   147,   148,   149,   150,   151,
     152,   153,   154,   155,   156,   157,   158,   159,   162,   163,
     165,   166,   169,   170,   173,   176,   177,   178,   181,   189,
     189,   189,   189,   189,   190,   190,   190,   190,   190,   191,
     191,   191,   191,   191,   192,   195,   196,   199,   205,   209,
     212,   213,   216,   217,   218,   219,   227,   235,   242,   253,
     261,   262,   265,   279,   282,   290,   294,   297,   298,   301,
     307,   310,   313,   317,   323,   324,   325,   326,   327,   328,
     334,   340,   346,   358,   359,   362,   366,   369,   372,   378,
     385,   390,   395,   403,   412,   421,   433,   434,   435,   438,
     439,   442,   443,   448,   456,   461,   469,   481,   482,   485,
     493,   500,   507,   516,   519,   533,   536,   537,   540,   543,
     544,   547,   566,   587,   610,   611,   614,   615,   616,   619,
     622,   623,   624,   625
};
#endif

#if YYDEBUG || YYERROR_VERBOSE || YYTOKEN_TABLE
/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "INT", "UINT", "REAL", "QSTR", "IDENT",
  "ID_DATA", "ID_X", "ID_VARIOGRAM", "ID_PREDICTIONS", "ID_VARIANCES",
  "ID_COVARIANCES", "ID_OUTPUT", "ID_MASKS", "ID_EDGES", "ID_SET",
  "ID_MERGE", "ID_AREA", "ID_BLOCK", "ID_METHOD", "ID_BOUNDS",
  "ID_MARGINALS", "';'", "'='", "':'", "','", "'['", "']'", "'('", "')'",
  "'&'", "'+'", "'-'", "'@'", "$accept", "input", "command", "val",
  "assign", "comcol", "d_vector", "d_list", "d_val", "any_id", "data_cmd",
  "data_decl", "data_cont", "data_exp", "data_X", "data_X_what",
  "data_what", "block_cmd", "block_cont", "block_exp", "block_lhs",
  "area_cmd", "area_decl", "vgm_cmd", "vgm_decl", "vgm_cont", "vgm_model",
  "vgm_model_type", "vgm_range", "sill_val", "range_val", "output_cmd",
  "set_cmd", "set_exp", "set_lhs", "method_cmd", "mask_cmd", "mask_cont",
  "edges_cmd", "edges_cont", "merge_cmd", "bounds_cmd", "bounds_exp",
  "marginals_cmd", "marginals_cont", 0
};
#endif

# ifdef YYPRINT
/* YYTOKNUM[YYLEX-NUM] -- Internal token number corresponding to
   token YYLEX-NUM.  */
static const yytype_uint16 yytoknum[] =
{
       0,   256,   257,   258,   259,   260,   261,   262,   263,   264,
     265,   266,   267,   268,   269,   270,   271,   272,   273,   274,
     275,   276,   277,   278,    59,    61,    58,    44,    91,    93,
      40,    41,    38,    43,    45,    64
};
# endif

/* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint8 yyr1[] =
{
       0,    36,    37,    37,    37,    38,    38,    38,    38,    38,
      38,    38,    38,    38,    38,    38,    38,    38,    39,    39,
      40,    40,    41,    41,    42,    43,    43,    43,    44,    45,
      45,    45,    45,    45,    45,    45,    45,    45,    45,    45,
      45,    45,    45,    45,    45,    46,    46,    47,    47,    47,
      48,    48,    49,    49,    49,    49,    49,    49,    49,    49,
      50,    50,    51,    51,    52,    53,    53,    54,    54,    55,
      56,    57,    58,    58,    59,    59,    59,    59,    59,    59,
      60,    60,    60,    61,    61,    62,    62,    62,    62,    63,
      64,    64,    64,    64,    64,    64,    65,    65,    65,    66,
      66,    67,    67,    67,    67,    67,    67,    68,    68,    69,
      69,    69,    69,    70,    71,    72,    73,    73,    74,    75,
      75,    76,    76,    76,    77,    77,    78,    78,    78,    79,
      80,    80,    80,    80
};

/* YYR2[YYN] -- Number of symbols composing right hand side of rule YYN.  */
static const yytype_uint8 yyr2[] =
{
       0,     2,     0,     2,     1,     1,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     1,     1,
       1,     1,     1,     1,     3,     1,     3,     2,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     3,     1,     4,     3,     3,
       1,     3,     0,     3,     1,     3,     3,     3,     3,     1,
       1,     3,     1,     1,     1,     1,     3,     1,     3,     3,
       1,     3,     1,     3,     5,     7,     3,     5,     3,     4,
       3,     4,     6,     1,     2,     4,     5,     6,     6,     1,
       1,     3,     5,     7,    11,    13,     1,     2,     2,     1,
       2,     3,     6,     3,     6,     3,     8,     2,     1,     3,
       3,     3,     3,     1,     3,     3,     1,     3,     3,     1,
       3,     4,    10,    10,     3,     3,     1,     2,     3,     3,
       1,     3,     1,     3
};

/* YYDEFACT[STATE-NAME] -- Default rule to reduce with in state
   STATE-NUM when YYTABLE doesn't specify something else to do.  Zero
   means the default is an error.  */
static const yytype_uint8 yydefact[] =
{
       2,    18,    19,     0,    28,     4,    25,     1,   113,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,    72,
      65,     0,     0,     0,     5,     3,     0,    46,     0,     0,
       0,     0,     0,     0,     0,   108,     0,     0,     0,     0,
       0,     0,     0,     0,    27,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   107,     0,    29,    34,    43,
      42,    39,    41,    33,    38,    36,    44,    40,    30,    31,
      37,    32,    35,     0,     0,     0,     0,     0,     0,     6,
      52,    11,    12,    52,     7,     0,    13,    17,    20,    21,
       0,    14,     9,    10,     8,    15,    16,    26,    49,    48,
       0,    80,     0,   103,     0,   105,     0,     0,   101,   116,
     115,   119,   118,     0,     0,     0,    73,    70,    66,    67,
       0,   114,   124,   126,   125,   132,   130,   129,    54,    64,
       0,    45,    50,    59,    71,     0,    76,     0,     0,     0,
      96,    78,    83,     0,   109,   110,   111,   112,    47,     0,
      81,     0,     0,     0,     0,     0,     0,   121,     0,     0,
       0,     0,   127,     0,     0,    52,     0,    79,    22,    23,
       0,     0,     0,    97,    98,    84,    89,     0,     0,     0,
       0,     0,   117,   120,     0,     0,    68,    69,   128,   133,
     131,    63,    62,    53,    60,    51,    55,    56,    57,     0,
      58,    77,    74,     0,     0,     0,    82,   102,   104,     0,
       0,     0,     0,     0,     0,     0,     0,    85,     0,    99,
       0,    90,     0,     0,     0,    61,    24,    75,     0,     0,
     100,    86,     0,   106,     0,     0,    87,    88,    91,     0,
       0,     0,   123,   122,    92,     0,    93,     0,     0,     0,
      94,     0,    95
};

/* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int16 yydefgoto[] =
{
      -1,     3,    25,   140,    90,   170,   200,     5,     6,    73,
      26,    27,   131,   132,   193,   194,   133,    28,   118,   119,
     120,    29,    30,    31,    32,   141,   142,   177,   220,   143,
     221,    33,    34,    35,    36,    37,    38,   110,    39,   112,
      40,    41,   124,    42,   127
};

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
#define YYPACT_NINF -164
static const yytype_int16 yypact[] =
{
     115,  -164,  -164,   130,  -164,    20,  -164,  -164,  -164,   -20,
     -15,    31,    38,    37,    46,    54,    59,    93,   214,    78,
      80,    89,    97,   107,  -164,  -164,    95,   108,   105,   112,
     109,   131,   132,   133,   135,  -164,    96,   138,   140,   141,
     142,   143,   144,   115,  -164,    82,   163,   150,   245,   175,
     245,   245,   181,   182,   183,  -164,   160,  -164,  -164,  -164,
    -164,  -164,  -164,  -164,  -164,  -164,  -164,  -164,  -164,  -164,
    -164,  -164,  -164,    12,   161,   184,   186,    76,    81,  -164,
     103,  -164,  -164,   103,  -164,    11,  -164,  -164,  -164,  -164,
      70,  -164,  -164,  -164,  -164,  -164,  -164,  -164,  -164,  -164,
     164,  -164,    39,  -164,   165,  -164,   166,   171,  -164,  -164,
     173,  -164,   174,   199,   245,   200,  -164,  -164,   177,  -164,
     185,  -164,  -164,  -164,    33,  -164,  -164,   179,  -164,  -164,
     187,   180,  -164,   189,   180,   195,    14,    23,    23,   115,
     176,    19,  -164,   201,  -164,  -164,  -164,  -164,  -164,   245,
    -164,   194,   212,   245,   233,   234,   211,  -164,   213,   184,
     115,   115,  -164,   111,    62,   103,    28,  -164,  -164,  -164,
      15,   201,   201,  -164,  -164,  -164,  -164,   216,   217,   237,
     244,   232,  -164,  -164,   262,   263,  -164,  -164,  -164,  -164,
    -164,  -164,  -164,   239,  -164,  -164,  -164,  -164,  -164,   115,
    -164,    14,    19,   242,   243,    24,  -164,  -164,  -164,   248,
     272,   245,    62,     8,    19,    27,    27,  -164,   115,  -164,
     246,   249,   273,   250,   251,  -164,  -164,    19,   247,   252,
    -164,  -164,   115,  -164,   279,   281,  -164,  -164,   258,   255,
     256,   115,  -164,  -164,   261,   115,   264,   115,   265,   115,
     266,   115,  -164
};

/* YYPGOTO[NTERM-NUM].  */
static const yytype_int16 yypgoto[] =
{
    -164,  -164,  -164,     0,  -164,    88,  -164,    91,    -4,   -42,
    -164,  -164,   215,   129,  -164,    83,  -164,  -164,  -164,   137,
    -164,  -164,  -164,  -164,  -164,  -163,  -139,   -46,   -88,    -6,
    -164,  -164,  -164,   280,  -164,  -164,  -164,  -164,  -164,  -164,
    -164,  -164,  -164,  -164,  -164
};

/* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule which
   number is the opposite.  If zero, do what YYDEFACT says.
   If YYTABLE_NINF, syntax error.  */
#define YYTABLE_NINF -1
static const yytype_uint8 yytable[] =
{
       4,    44,   175,   100,   102,     4,   104,   202,   106,   107,
      45,     1,   135,     2,     1,    46,     2,   136,     1,   114,
       2,   201,     1,     1,     2,     2,     1,     1,     2,     2,
       1,   196,     2,   197,   198,    43,     1,   226,     2,    97,
     168,   169,   115,     4,   137,   138,   139,    43,   137,   138,
     139,   227,   137,   138,   139,   217,   199,    47,   139,   218,
     161,    48,   218,   175,    49,   191,   149,    51,    50,   192,
     150,    52,   157,   144,   145,   146,   147,   123,   126,     1,
      53,     2,   122,    98,     1,    54,     2,   125,   175,    57,
      58,    59,    60,    61,    62,    63,    64,    65,    66,    67,
       8,    68,    69,    70,    71,    72,    75,   178,    74,   128,
     129,   181,   130,    99,     1,    76,     2,   189,     1,    79,
       2,    88,    89,    77,   162,   203,   204,   228,   229,    81,
       7,   171,   172,    78,    80,    83,    82,     8,     9,   173,
      10,    11,    12,    13,    14,    15,    16,    17,    18,    19,
      20,    21,    22,    23,    24,    84,   103,    86,    85,    87,
     187,   188,    91,   190,    92,    93,    94,    95,    96,   224,
      57,    58,    59,    60,    61,    62,    63,    64,    65,    66,
      67,   105,    68,    69,    70,    71,    72,   108,   109,   111,
     113,   117,   116,   121,   101,   148,   151,   152,   153,     4,
     154,   155,   156,   158,   159,   219,   163,   165,   176,    44,
     160,   174,   164,     4,   166,   219,   219,    56,   230,   167,
     179,    57,    58,    59,    60,    61,    62,    63,    64,    65,
      66,    67,   238,    68,    69,    70,    71,    72,   180,   182,
     183,   244,   184,   207,   185,   246,   205,   248,   206,   250,
     208,   252,    57,    58,    59,    60,    61,    62,    63,    64,
      65,    66,    67,   209,    68,    69,    70,    71,    72,   210,
     211,   212,   215,   216,   222,   223,   232,   231,   236,   233,
     234,   235,   239,   237,   240,   241,   242,   243,   245,   214,
     213,   247,   249,   251,   195,   225,   186,    55,   134
};

static const yytype_uint8 yycheck[] =
{
       0,     5,   141,    45,    46,     5,    48,   170,    50,    51,
      30,     3,     1,     5,     3,    30,     5,     6,     3,     7,
       5,     6,     3,     3,     5,     5,     3,     3,     5,     5,
       3,     3,     5,     5,     6,    27,     3,    29,     5,    43,
      26,    27,    30,    43,    33,    34,    35,    27,    33,    34,
      35,   214,    33,    34,    35,    31,    28,    26,    35,    35,
      27,    30,    35,   202,    26,     3,    27,    30,    30,     7,
      31,    25,   114,     3,     4,     5,     6,    77,    78,     3,
      26,     5,     6,     1,     3,    26,     5,     6,   227,     7,
       8,     9,    10,    11,    12,    13,    14,    15,    16,    17,
       7,    19,    20,    21,    22,    23,    26,   149,    30,     6,
       7,   153,     9,    31,     3,    26,     5,     6,     3,    24,
       5,    25,    26,    26,   124,   171,   172,   215,   216,    24,
       0,   137,   138,    26,    26,    26,    24,     7,     8,   139,
      10,    11,    12,    13,    14,    15,    16,    17,    18,    19,
      20,    21,    22,    23,    24,    24,     6,    24,    26,    24,
     160,   161,    24,   163,    24,    24,    24,    24,    24,   211,
       7,     8,     9,    10,    11,    12,    13,    14,    15,    16,
      17,     6,    19,    20,    21,    22,    23,     6,     6,     6,
      30,     7,    31,     7,    31,    31,    31,    31,    27,   199,
      27,    27,     3,     3,    27,   205,    27,    27,     7,   213,
      25,    35,    25,   213,    25,   215,   216,     3,   218,    24,
      26,     7,     8,     9,    10,    11,    12,    13,    14,    15,
      16,    17,   232,    19,    20,    21,    22,    23,    26,     6,
       6,   241,    31,     6,    31,   245,    30,   247,    31,   249,
       6,   251,     7,     8,     9,    10,    11,    12,    13,    14,
      15,    16,    17,    31,    19,    20,    21,    22,    23,     7,
       7,    32,    30,    30,    26,     3,    27,    31,    31,     6,
      30,    30,     3,    31,     3,    27,    31,    31,    27,   201,
     199,    27,    27,    27,   165,   212,   159,    17,    83
};

/* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
   symbol of state STATE-NUM.  */
static const yytype_uint8 yystos[] =
{
       0,     3,     5,    37,    39,    43,    44,     0,     7,     8,
      10,    11,    12,    13,    14,    15,    16,    17,    18,    19,
      20,    21,    22,    23,    24,    38,    46,    47,    53,    57,
      58,    59,    60,    67,    68,    69,    70,    71,    72,    74,
      76,    77,    79,    27,    44,    30,    30,    26,    30,    26,
      30,    30,    25,    26,    26,    69,     3,     7,     8,     9,
      10,    11,    12,    13,    14,    15,    16,    17,    19,    20,
      21,    22,    23,    45,    30,    26,    26,    26,    26,    24,
      26,    24,    24,    26,    24,    26,    24,    24,    25,    26,
      40,    24,    24,    24,    24,    24,    24,    44,     1,    31,
      45,    31,    45,     6,    45,     6,    45,    45,     6,     6,
      73,     6,    75,    30,     7,    30,    31,     7,    54,    55,
      56,     7,     6,    39,    78,     6,    39,    80,     6,     7,
       9,    48,    49,    52,    48,     1,     6,    33,    34,    35,
      39,    61,    62,    65,     3,     4,     5,     6,    31,    27,
      31,    31,    31,    27,    27,    27,     3,    45,     3,    27,
      25,    27,    39,    27,    25,    27,    25,    24,    26,    27,
      41,    65,    65,    39,    35,    62,     7,    63,    45,    26,
      26,    45,     6,     6,    31,    31,    55,    39,    39,     6,
      39,     3,     7,    50,    51,    49,     3,     5,     6,    28,
      42,     6,    61,    63,    63,    30,    31,     6,     6,    31,
       7,     7,    32,    43,    41,    30,    30,    31,    35,    39,
      64,    66,    26,     3,    45,    51,    29,    61,    64,    64,
      39,    31,    27,     6,    30,    30,    31,    31,    39,     3,
       3,    27,    31,    31,    39,    27,    39,    27,    39,    27,
      39,    27,    39
};

#define yyerrok		(yyerrstatus = 0)
#define yyclearin	(yychar = YYEMPTY)
#define YYEMPTY		(-2)
#define YYEOF		0

#define YYACCEPT	goto yyacceptlab
#define YYABORT		goto yyabortlab
#define YYERROR		goto yyerrorlab


/* Like YYERROR except do call yyerror.  This remains here temporarily
   to ease the transition to the new meaning of YYERROR, for GCC.
   Once GCC version 2 has supplanted version 1, this can go.  */

#define YYFAIL		goto yyerrlab

#define YYRECOVERING()  (!!yyerrstatus)

#define YYBACKUP(Token, Value)					\
do								\
  if (yychar == YYEMPTY && yylen == 1)				\
    {								\
      yychar = (Token);						\
      yylval = (Value);						\
      yytoken = YYTRANSLATE (yychar);				\
      YYPOPSTACK (1);						\
      goto yybackup;						\
    }								\
  else								\
    {								\
      yyerror (YY_("syntax error: cannot back up")); \
      YYERROR;							\
    }								\
while (YYID (0))


#define YYTERROR	1
#define YYERRCODE	256


/* YYLLOC_DEFAULT -- Set CURRENT to span from RHS[1] to RHS[N].
   If N is 0, then set CURRENT to the empty location which ends
   the previous symbol: RHS[0] (always defined).  */

#define YYRHSLOC(Rhs, K) ((Rhs)[K])
#ifndef YYLLOC_DEFAULT
# define YYLLOC_DEFAULT(Current, Rhs, N)				\
    do									\
      if (YYID (N))                                                    \
	{								\
	  (Current).first_line   = YYRHSLOC (Rhs, 1).first_line;	\
	  (Current).first_column = YYRHSLOC (Rhs, 1).first_column;	\
	  (Current).last_line    = YYRHSLOC (Rhs, N).last_line;		\
	  (Current).last_column  = YYRHSLOC (Rhs, N).last_column;	\
	}								\
      else								\
	{								\
	  (Current).first_line   = (Current).last_line   =		\
	    YYRHSLOC (Rhs, 0).last_line;				\
	  (Current).first_column = (Current).last_column =		\
	    YYRHSLOC (Rhs, 0).last_column;				\
	}								\
    while (YYID (0))
#endif


/* YY_LOCATION_PRINT -- Print the location on the stream.
   This macro was not mandated originally: define only if we know
   we won't break user code: when these are the locations we know.  */

#ifndef YY_LOCATION_PRINT
# if YYLTYPE_IS_TRIVIAL
#  define YY_LOCATION_PRINT(File, Loc)			\
     fprintf (File, "%d.%d-%d.%d",			\
	      (Loc).first_line, (Loc).first_column,	\
	      (Loc).last_line,  (Loc).last_column)
# else
#  define YY_LOCATION_PRINT(File, Loc) ((void) 0)
# endif
#endif


/* YYLEX -- calling `yylex' with the right arguments.  */

#ifdef YYLEX_PARAM
# define YYLEX yylex (YYLEX_PARAM)
#else
# define YYLEX yylex ()
#endif

/* Enable debugging if requested.  */
#if YYDEBUG

# ifndef YYFPRINTF
#  include <stdio.h> /* INFRINGES ON USER NAME SPACE */
#  define YYFPRINTF fprintf
# endif

# define YYDPRINTF(Args)			\
do {						\
  if (yydebug)					\
    YYFPRINTF Args;				\
} while (YYID (0))

# define YY_SYMBOL_PRINT(Title, Type, Value, Location)			  \
do {									  \
  if (yydebug)								  \
    {									  \
      YYFPRINTF (stderr, "%s ", Title);					  \
      yy_symbol_print (stderr,						  \
		  Type, Value); \
      YYFPRINTF (stderr, "\n");						  \
    }									  \
} while (YYID (0))


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

/*ARGSUSED*/
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_symbol_value_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep)
#else
static void
yy_symbol_value_print (yyoutput, yytype, yyvaluep)
    FILE *yyoutput;
    int yytype;
    YYSTYPE const * const yyvaluep;
#endif
{
  if (!yyvaluep)
    return;
# ifdef YYPRINT
  if (yytype < YYNTOKENS)
    YYPRINT (yyoutput, yytoknum[yytype], *yyvaluep);
# else
  YYUSE (yyoutput);
# endif
  switch (yytype)
    {
      default:
	break;
    }
}


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_symbol_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep)
#else
static void
yy_symbol_print (yyoutput, yytype, yyvaluep)
    FILE *yyoutput;
    int yytype;
    YYSTYPE const * const yyvaluep;
#endif
{
  if (yytype < YYNTOKENS)
    YYFPRINTF (yyoutput, "token %s (", yytname[yytype]);
  else
    YYFPRINTF (yyoutput, "nterm %s (", yytname[yytype]);

  yy_symbol_value_print (yyoutput, yytype, yyvaluep);
  YYFPRINTF (yyoutput, ")");
}

/*------------------------------------------------------------------.
| yy_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (included).                                                   |
`------------------------------------------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_stack_print (yytype_int16 *yybottom, yytype_int16 *yytop)
#else
static void
yy_stack_print (yybottom, yytop)
    yytype_int16 *yybottom;
    yytype_int16 *yytop;
#endif
{
  YYFPRINTF (stderr, "Stack now");
  for (; yybottom <= yytop; yybottom++)
    {
      int yybot = *yybottom;
      YYFPRINTF (stderr, " %d", yybot);
    }
  YYFPRINTF (stderr, "\n");
}

# define YY_STACK_PRINT(Bottom, Top)				\
do {								\
  if (yydebug)							\
    yy_stack_print ((Bottom), (Top));				\
} while (YYID (0))


/*------------------------------------------------.
| Report that the YYRULE is going to be reduced.  |
`------------------------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_reduce_print (YYSTYPE *yyvsp, int yyrule)
#else
static void
yy_reduce_print (yyvsp, yyrule)
    YYSTYPE *yyvsp;
    int yyrule;
#endif
{
  int yynrhs = yyr2[yyrule];
  int yyi;
  unsigned long int yylno = yyrline[yyrule];
  YYFPRINTF (stderr, "Reducing stack by rule %d (line %lu):\n",
	     yyrule - 1, yylno);
  /* The symbols being reduced.  */
  for (yyi = 0; yyi < yynrhs; yyi++)
    {
      YYFPRINTF (stderr, "   $%d = ", yyi + 1);
      yy_symbol_print (stderr, yyrhs[yyprhs[yyrule] + yyi],
		       &(yyvsp[(yyi + 1) - (yynrhs)])
		       		       );
      YYFPRINTF (stderr, "\n");
    }
}

# define YY_REDUCE_PRINT(Rule)		\
do {					\
  if (yydebug)				\
    yy_reduce_print (yyvsp, Rule); \
} while (YYID (0))

/* Nonzero means print parse trace.  It is left uninitialized so that
   multiple parsers can coexist.  */
int yydebug;
#else /* !YYDEBUG */
# define YYDPRINTF(Args)
# define YY_SYMBOL_PRINT(Title, Type, Value, Location)
# define YY_STACK_PRINT(Bottom, Top)
# define YY_REDUCE_PRINT(Rule)
#endif /* !YYDEBUG */


/* YYINITDEPTH -- initial size of the parser's stacks.  */
#ifndef	YYINITDEPTH
# define YYINITDEPTH 200
#endif

/* YYMAXDEPTH -- maximum size the stacks can grow to (effective only
   if the built-in stack extension method is used).

   Do not make this value too large; the results are undefined if
   YYSTACK_ALLOC_MAXIMUM < YYSTACK_BYTES (YYMAXDEPTH)
   evaluated with infinite-precision integer arithmetic.  */

#ifndef YYMAXDEPTH
# define YYMAXDEPTH 10000
#endif



#if YYERROR_VERBOSE

# ifndef yystrlen
#  if defined __GLIBC__ && defined _STRING_H
#   define yystrlen strlen
#  else
/* Return the length of YYSTR.  */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static YYSIZE_T
yystrlen (const char *yystr)
#else
static YYSIZE_T
yystrlen (yystr)
    const char *yystr;
#endif
{
  YYSIZE_T yylen;
  for (yylen = 0; yystr[yylen]; yylen++)
    continue;
  return yylen;
}
#  endif
# endif

# ifndef yystpcpy
#  if defined __GLIBC__ && defined _STRING_H && defined _GNU_SOURCE
#   define yystpcpy stpcpy
#  else
/* Copy YYSRC to YYDEST, returning the address of the terminating '\0' in
   YYDEST.  */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static char *
yystpcpy (char *yydest, const char *yysrc)
#else
static char *
yystpcpy (yydest, yysrc)
    char *yydest;
    const char *yysrc;
#endif
{
  char *yyd = yydest;
  const char *yys = yysrc;

  while ((*yyd++ = *yys++) != '\0')
    continue;

  return yyd - 1;
}
#  endif
# endif

# ifndef yytnamerr
/* Copy to YYRES the contents of YYSTR after stripping away unnecessary
   quotes and backslashes, so that it's suitable for yyerror.  The
   heuristic is that double-quoting is unnecessary unless the string
   contains an apostrophe, a comma, or backslash (other than
   backslash-backslash).  YYSTR is taken from yytname.  If YYRES is
   null, do not copy; instead, return the length of what the result
   would have been.  */
static YYSIZE_T
yytnamerr (char *yyres, const char *yystr)
{
  if (*yystr == '"')
    {
      YYSIZE_T yyn = 0;
      char const *yyp = yystr;

      for (;;)
	switch (*++yyp)
	  {
	  case '\'':
	  case ',':
	    goto do_not_strip_quotes;

	  case '\\':
	    if (*++yyp != '\\')
	      goto do_not_strip_quotes;
	    /* Fall through.  */
	  default:
	    if (yyres)
	      yyres[yyn] = *yyp;
	    yyn++;
	    break;

	  case '"':
	    if (yyres)
	      yyres[yyn] = '\0';
	    return yyn;
	  }
    do_not_strip_quotes: ;
    }

  if (! yyres)
    return yystrlen (yystr);

  return yystpcpy (yyres, yystr) - yyres;
}
# endif

/* Copy into YYRESULT an error message about the unexpected token
   YYCHAR while in state YYSTATE.  Return the number of bytes copied,
   including the terminating null byte.  If YYRESULT is null, do not
   copy anything; just return the number of bytes that would be
   copied.  As a special case, return 0 if an ordinary "syntax error"
   message will do.  Return YYSIZE_MAXIMUM if overflow occurs during
   size calculation.  */
static YYSIZE_T
yysyntax_error (char *yyresult, int yystate, int yychar)
{
  int yyn = yypact[yystate];

  if (! (YYPACT_NINF < yyn && yyn <= YYLAST))
    return 0;
  else
    {
      int yytype = YYTRANSLATE (yychar);
      YYSIZE_T yysize0 = yytnamerr (0, yytname[yytype]);
      YYSIZE_T yysize = yysize0;
      YYSIZE_T yysize1;
      int yysize_overflow = 0;
      enum { YYERROR_VERBOSE_ARGS_MAXIMUM = 5 };
      char const *yyarg[YYERROR_VERBOSE_ARGS_MAXIMUM];
      int yyx;

# if 0
      /* This is so xgettext sees the translatable formats that are
	 constructed on the fly.  */
      YY_("syntax error, unexpected %s");
      YY_("syntax error, unexpected %s, expecting %s");
      YY_("syntax error, unexpected %s, expecting %s or %s");
      YY_("syntax error, unexpected %s, expecting %s or %s or %s");
      YY_("syntax error, unexpected %s, expecting %s or %s or %s or %s");
# endif
      char *yyfmt;
      char const *yyf;
      static char const yyunexpected[] = "syntax error, unexpected %s";
      static char const yyexpecting[] = ", expecting %s";
      static char const yyor[] = " or %s";
      char yyformat[sizeof yyunexpected
		    + sizeof yyexpecting - 1
		    + ((YYERROR_VERBOSE_ARGS_MAXIMUM - 2)
		       * (sizeof yyor - 1))];
      char const *yyprefix = yyexpecting;

      /* Start YYX at -YYN if negative to avoid negative indexes in
	 YYCHECK.  */
      int yyxbegin = yyn < 0 ? -yyn : 0;

      /* Stay within bounds of both yycheck and yytname.  */
      int yychecklim = YYLAST - yyn + 1;
      int yyxend = yychecklim < YYNTOKENS ? yychecklim : YYNTOKENS;
      int yycount = 1;

      yyarg[0] = yytname[yytype];
      yyfmt = yystpcpy (yyformat, yyunexpected);

      for (yyx = yyxbegin; yyx < yyxend; ++yyx)
	if (yycheck[yyx + yyn] == yyx && yyx != YYTERROR)
	  {
	    if (yycount == YYERROR_VERBOSE_ARGS_MAXIMUM)
	      {
		yycount = 1;
		yysize = yysize0;
		yyformat[sizeof yyunexpected - 1] = '\0';
		break;
	      }
	    yyarg[yycount++] = yytname[yyx];
	    yysize1 = yysize + yytnamerr (0, yytname[yyx]);
	    yysize_overflow |= (yysize1 < yysize);
	    yysize = yysize1;
	    yyfmt = yystpcpy (yyfmt, yyprefix);
	    yyprefix = yyor;
	  }

      yyf = YY_(yyformat);
      yysize1 = yysize + yystrlen (yyf);
      yysize_overflow |= (yysize1 < yysize);
      yysize = yysize1;

      if (yysize_overflow)
	return YYSIZE_MAXIMUM;

      if (yyresult)
	{
	  /* Avoid sprintf, as that infringes on the user's name space.
	     Don't have undefined behavior even if the translation
	     produced a string with the wrong number of "%s"s.  */
	  char *yyp = yyresult;
	  int yyi = 0;
	  while ((*yyp = *yyf) != '\0')
	    {
	      if (*yyp == '%' && yyf[1] == 's' && yyi < yycount)
		{
		  yyp += yytnamerr (yyp, yyarg[yyi++]);
		  yyf += 2;
		}
	      else
		{
		  yyp++;
		  yyf++;
		}
	    }
	}
      return yysize;
    }
}
#endif /* YYERROR_VERBOSE */


/*-----------------------------------------------.
| Release the memory associated to this symbol.  |
`-----------------------------------------------*/

/*ARGSUSED*/
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yydestruct (const char *yymsg, int yytype, YYSTYPE *yyvaluep)
#else
static void
yydestruct (yymsg, yytype, yyvaluep)
    const char *yymsg;
    int yytype;
    YYSTYPE *yyvaluep;
#endif
{
  YYUSE (yyvaluep);

  if (!yymsg)
    yymsg = "Deleting";
  YY_SYMBOL_PRINT (yymsg, yytype, yyvaluep, yylocationp);

  switch (yytype)
    {

      default:
	break;
    }
}

/* Prevent warnings from -Wmissing-prototypes.  */
#ifdef YYPARSE_PARAM
#if defined __STDC__ || defined __cplusplus
int yyparse (void *YYPARSE_PARAM);
#else
int yyparse ();
#endif
#else /* ! YYPARSE_PARAM */
#if defined __STDC__ || defined __cplusplus
int yyparse (void);
#else
int yyparse ();
#endif
#endif /* ! YYPARSE_PARAM */


/* The lookahead symbol.  */
int yychar;

/* The semantic value of the lookahead symbol.  */
YYSTYPE yylval;

/* Number of syntax errors so far.  */
int yynerrs;



/*-------------------------.
| yyparse or yypush_parse.  |
`-------------------------*/

#ifdef YYPARSE_PARAM
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
int
yyparse (void *YYPARSE_PARAM)
#else
int
yyparse (YYPARSE_PARAM)
    void *YYPARSE_PARAM;
#endif
#else /* ! YYPARSE_PARAM */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
int
yyparse (void)
#else
int
yyparse ()

#endif
#endif
{


    int yystate;
    /* Number of tokens to shift before error messages enabled.  */
    int yyerrstatus;

    /* The stacks and their tools:
       `yyss': related to states.
       `yyvs': related to semantic values.

       Refer to the stacks thru separate pointers, to allow yyoverflow
       to reallocate them elsewhere.  */

    /* The state stack.  */
    yytype_int16 yyssa[YYINITDEPTH];
    yytype_int16 *yyss;
    yytype_int16 *yyssp;

    /* The semantic value stack.  */
    YYSTYPE yyvsa[YYINITDEPTH];
    YYSTYPE *yyvs;
    YYSTYPE *yyvsp;

    YYSIZE_T yystacksize;

  int yyn;
  int yyresult;
  /* Lookahead token as an internal (translated) token number.  */
  int yytoken;
  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE yyval;

#if YYERROR_VERBOSE
  /* Buffer for error messages, and its allocated size.  */
  char yymsgbuf[128];
  char *yymsg = yymsgbuf;
  YYSIZE_T yymsg_alloc = sizeof yymsgbuf;
#endif

#define YYPOPSTACK(N)   (yyvsp -= (N), yyssp -= (N))

  /* The number of symbols on the RHS of the reduced rule.
     Keep to zero when no symbol should be popped.  */
  int yylen = 0;

  yytoken = 0;
  yyss = yyssa;
  yyvs = yyvsa;
  yystacksize = YYINITDEPTH;

  YYDPRINTF ((stderr, "Starting parse\n"));

  yystate = 0;
  yyerrstatus = 0;
  yynerrs = 0;
  yychar = YYEMPTY; /* Cause a token to be read.  */

  /* Initialize stack pointers.
     Waste one element of value and location stack
     so that they stay on the same level as the state stack.
     The wasted elements are never initialized.  */
  yyssp = yyss;
  yyvsp = yyvs;

  goto yysetstate;

/*------------------------------------------------------------.
| yynewstate -- Push a new state, which is found in yystate.  |
`------------------------------------------------------------*/
 yynewstate:
  /* In all cases, when you get here, the value and location stacks
     have just been pushed.  So pushing a state here evens the stacks.  */
  yyssp++;

 yysetstate:
  *yyssp = yystate;

  if (yyss + yystacksize - 1 <= yyssp)
    {
      /* Get the current used size of the three stacks, in elements.  */
      YYSIZE_T yysize = yyssp - yyss + 1;

#ifdef yyoverflow
      {
	/* Give user a chance to reallocate the stack.  Use copies of
	   these so that the &'s don't force the real ones into
	   memory.  */
	YYSTYPE *yyvs1 = yyvs;
	yytype_int16 *yyss1 = yyss;

	/* Each stack pointer address is followed by the size of the
	   data in use in that stack, in bytes.  This used to be a
	   conditional around just the two extra args, but that might
	   be undefined if yyoverflow is a macro.  */
	yyoverflow (YY_("memory exhausted"),
		    &yyss1, yysize * sizeof (*yyssp),
		    &yyvs1, yysize * sizeof (*yyvsp),
		    &yystacksize);

	yyss = yyss1;
	yyvs = yyvs1;
      }
#else /* no yyoverflow */
# ifndef YYSTACK_RELOCATE
      goto yyexhaustedlab;
# else
      /* Extend the stack our own way.  */
      if (YYMAXDEPTH <= yystacksize)
	goto yyexhaustedlab;
      yystacksize *= 2;
      if (YYMAXDEPTH < yystacksize)
	yystacksize = YYMAXDEPTH;

      {
	yytype_int16 *yyss1 = yyss;
	union yyalloc *yyptr =
	  (union yyalloc *) YYSTACK_ALLOC (YYSTACK_BYTES (yystacksize));
	if (! yyptr)
	  goto yyexhaustedlab;
	YYSTACK_RELOCATE (yyss_alloc, yyss);
	YYSTACK_RELOCATE (yyvs_alloc, yyvs);
#  undef YYSTACK_RELOCATE
	if (yyss1 != yyssa)
	  YYSTACK_FREE (yyss1);
      }
# endif
#endif /* no yyoverflow */

      yyssp = yyss + yysize - 1;
      yyvsp = yyvs + yysize - 1;

      YYDPRINTF ((stderr, "Stack size increased to %lu\n",
		  (unsigned long int) yystacksize));

      if (yyss + yystacksize - 1 <= yyssp)
	YYABORT;
    }

  YYDPRINTF ((stderr, "Entering state %d\n", yystate));

  if (yystate == YYFINAL)
    YYACCEPT;

  goto yybackup;

/*-----------.
| yybackup.  |
`-----------*/
yybackup:

  /* Do appropriate processing given the current state.  Read a
     lookahead token if we need one and don't already have one.  */

  /* First try to decide what to do without reference to lookahead token.  */
  yyn = yypact[yystate];
  if (yyn == YYPACT_NINF)
    goto yydefault;

  /* Not known => get a lookahead token if don't already have one.  */

  /* YYCHAR is either YYEMPTY or YYEOF or a valid lookahead symbol.  */
  if (yychar == YYEMPTY)
    {
      YYDPRINTF ((stderr, "Reading a token: "));
      yychar = YYLEX;
    }

  if (yychar <= YYEOF)
    {
      yychar = yytoken = YYEOF;
      YYDPRINTF ((stderr, "Now at end of input.\n"));
    }
  else
    {
      yytoken = YYTRANSLATE (yychar);
      YY_SYMBOL_PRINT ("Next token is", yytoken, &yylval, &yylloc);
    }

  /* If the proper action on seeing token YYTOKEN is to reduce or to
     detect an error, take that action.  */
  yyn += yytoken;
  if (yyn < 0 || YYLAST < yyn || yycheck[yyn] != yytoken)
    goto yydefault;
  yyn = yytable[yyn];
  if (yyn <= 0)
    {
      if (yyn == 0 || yyn == YYTABLE_NINF)
	goto yyerrlab;
      yyn = -yyn;
      goto yyreduce;
    }

  /* Count tokens shifted since error; after three, turn off error
     status.  */
  if (yyerrstatus)
    yyerrstatus--;

  /* Shift the lookahead token.  */
  YY_SYMBOL_PRINT ("Shifting", yytoken, &yylval, &yylloc);

  /* Discard the shifted token.  */
  yychar = YYEMPTY;

  yystate = yyn;
  *++yyvsp = yylval;

  goto yynewstate;


/*-----------------------------------------------------------.
| yydefault -- do the default action for the current state.  |
`-----------------------------------------------------------*/
yydefault:
  yyn = yydefact[yystate];
  if (yyn == 0)
    goto yyerrlab;
  goto yyreduce;


/*-----------------------------.
| yyreduce -- Do a reduction.  |
`-----------------------------*/
yyreduce:
  /* yyn is the number of a rule to reduce with.  */
  yylen = yyr2[yyn];

  /* If YYLEN is nonzero, implement the default value of the action:
     `$$ = $1'.

     Otherwise, the following line sets YYVAL to garbage.
     This behavior is undocumented and Bison
     users should not rely upon it.  Assigning to YYVAL
     unconditionally makes the parser a bit smaller, and it avoids a
     GCC warning that YYVAL may be used uninitialized.  */
  yyval = yyvsp[1-yylen];


  YY_REDUCE_PRINT (yyn);
  switch (yyn)
    {
        case 2:

/* Line 1455 of yacc.c  */
#line 142 "parse.y"
    { ; }
    break;

  case 3:

/* Line 1455 of yacc.c  */
#line 143 "parse.y"
    { reset_parser(); }
    break;

  case 4:

/* Line 1455 of yacc.c  */
#line 144 "parse.y"
    { vector_only = 1; }
    break;

  case 5:

/* Line 1455 of yacc.c  */
#line 147 "parse.y"
    { ; }
    break;

  case 6:

/* Line 1455 of yacc.c  */
#line 148 "parse.y"
    { ; }
    break;

  case 7:

/* Line 1455 of yacc.c  */
#line 149 "parse.y"
    { update_variogram(v); }
    break;

  case 8:

/* Line 1455 of yacc.c  */
#line 150 "parse.y"
    { ; }
    break;

  case 9:

/* Line 1455 of yacc.c  */
#line 151 "parse.y"
    { ; }
    break;

  case 10:

/* Line 1455 of yacc.c  */
#line 152 "parse.y"
    { ; }
    break;

  case 11:

/* Line 1455 of yacc.c  */
#line 153 "parse.y"
    { ; }
    break;

  case 12:

/* Line 1455 of yacc.c  */
#line 154 "parse.y"
    { ; }
    break;

  case 13:

/* Line 1455 of yacc.c  */
#line 155 "parse.y"
    { ; }
    break;

  case 14:

/* Line 1455 of yacc.c  */
#line 156 "parse.y"
    { ; }
    break;

  case 15:

/* Line 1455 of yacc.c  */
#line 157 "parse.y"
    { ; }
    break;

  case 16:

/* Line 1455 of yacc.c  */
#line 158 "parse.y"
    { ; }
    break;

  case 17:

/* Line 1455 of yacc.c  */
#line 159 "parse.y"
    { ; }
    break;

  case 18:

/* Line 1455 of yacc.c  */
#line 162 "parse.y"
    { (yyval.dval) = (double) (yyvsp[(1) - (1)].ival); }
    break;

  case 20:

/* Line 1455 of yacc.c  */
#line 165 "parse.y"
    { ; }
    break;

  case 21:

/* Line 1455 of yacc.c  */
#line 166 "parse.y"
    { ; }
    break;

  case 22:

/* Line 1455 of yacc.c  */
#line 169 "parse.y"
    { ; }
    break;

  case 23:

/* Line 1455 of yacc.c  */
#line 170 "parse.y"
    { ; }
    break;

  case 24:

/* Line 1455 of yacc.c  */
#line 173 "parse.y"
    { ; }
    break;

  case 28:

/* Line 1455 of yacc.c  */
#line 181 "parse.y"
    { 
			if (d == NULL)
				sd_vector = push_d_vector((yyvsp[(1) - (1)].dval), sd_vector);
			else
				d->beta = push_d_vector((yyvsp[(1) - (1)].dval), d->beta);
		}
    break;

  case 44:

/* Line 1455 of yacc.c  */
#line 192 "parse.y"
    { ; }
    break;

  case 45:

/* Line 1455 of yacc.c  */
#line 195 "parse.y"
    { verify_data(d); }
    break;

  case 46:

/* Line 1455 of yacc.c  */
#line 196 "parse.y"
    { d->dummy = 1; }
    break;

  case 47:

/* Line 1455 of yacc.c  */
#line 199 "parse.y"
    {
			id = which_identifier((yyvsp[(3) - (4)].sval));
			dpp = get_gstat_data();
			d = dpp[id];
			d->id = id;
		}
    break;

  case 48:

/* Line 1455 of yacc.c  */
#line 205 "parse.y"
    {
			d = get_dataval();
			d->id = ID_OF_VALDATA;
		}
    break;

  case 49:

/* Line 1455 of yacc.c  */
#line 209 "parse.y"
    { ErrMsg(ER_SYNTAX, "invalid identifier"); }
    break;

  case 52:

/* Line 1455 of yacc.c  */
#line 216 "parse.y"
    { ; }
    break;

  case 54:

/* Line 1455 of yacc.c  */
#line 218 "parse.y"
    { d->fname = (yyvsp[(1) - (1)].sval); }
    break;

  case 55:

/* Line 1455 of yacc.c  */
#line 219 "parse.y"
    {
			switch (expr.what) { 
				case IS_INT: *((int *)expr.ptr) = (yyvsp[(3) - (3)].ival); break;
				case IS_REAL: *((double *)expr.ptr) = (double) (yyvsp[(3) - (3)].ival); break;
				default: lex_error(); YYERROR; break;
			}
			check_assign_expr(&expr);
		}
    break;

  case 56:

/* Line 1455 of yacc.c  */
#line 227 "parse.y"
    {
			if (expr.what != IS_REAL) {
				lex_error();
				YYERROR;
			}
			*((double *)expr.ptr) = (yyvsp[(3) - (3)].dval);
			check_assign_expr(&expr);
		}
    break;

  case 57:

/* Line 1455 of yacc.c  */
#line 235 "parse.y"
    {
			if (expr.what != IS_STRING) {
				lex_error();
				YYERROR;
			}
			*((char **)expr.ptr) = (yyvsp[(3) - (3)].sval);
		}
    break;

  case 58:

/* Line 1455 of yacc.c  */
#line 242 "parse.y"
    {
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
    break;

  case 59:

/* Line 1455 of yacc.c  */
#line 253 "parse.y"
    {
			if (expr.what != NO_ARG) {
				lex_error();
				YYERROR;
			}
		}
    break;

  case 62:

/* Line 1455 of yacc.c  */
#line 265 "parse.y"
    {
			for (id = 0; id < N_POLY; id++) {
				if (almost_equals((yyvsp[(1) - (1)].sval), polynomial[id].name)) {
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
    break;

  case 63:

/* Line 1455 of yacc.c  */
#line 279 "parse.y"
    { push_data_X(d, (yyvsp[(1) - (1)].ival)); }
    break;

  case 64:

/* Line 1455 of yacc.c  */
#line 282 "parse.y"
    {
			if (! is_data_expr(d, &expr, (yyvsp[(1) - (1)].sval))) {
				lex_error();
				YYERROR;
			}
		}
    break;

  case 65:

/* Line 1455 of yacc.c  */
#line 290 "parse.y"
    {
			bp = get_block_p();
			bp->x = -1.0; /* will be set to grid cell size in predict.c */
		}
    break;

  case 69:

/* Line 1455 of yacc.c  */
#line 301 "parse.y"
    {
			*((double *)expr.ptr) = (yyvsp[(3) - (3)].dval);
			check_assign_expr(&expr);
		}
    break;

  case 70:

/* Line 1455 of yacc.c  */
#line 307 "parse.y"
    { if (! is_block_expr(&expr, (yyvsp[(1) - (1)].sval))) { lex_error(); YYERROR; }}
    break;

  case 71:

/* Line 1455 of yacc.c  */
#line 310 "parse.y"
    { ; }
    break;

  case 72:

/* Line 1455 of yacc.c  */
#line 313 "parse.y"
    {
			d = create_data_area();
			d->id = ID_OF_AREA;
		}
    break;

  case 73:

/* Line 1455 of yacc.c  */
#line 317 "parse.y"
    {
			d = create_data_area();
			d->id = ID_OF_AREA;
		}
    break;

  case 74:

/* Line 1455 of yacc.c  */
#line 323 "parse.y"
    { v->fname = (yyvsp[(3) - (5)].sval); }
    break;

  case 75:

/* Line 1455 of yacc.c  */
#line 324 "parse.y"
    { v->fname = (yyvsp[(3) - (7)].sval); v->fname2 = (yyvsp[(5) - (7)].sval); }
    break;

  case 76:

/* Line 1455 of yacc.c  */
#line 325 "parse.y"
    {v->fname = (yyvsp[(3) - (3)].sval); }
    break;

  case 77:

/* Line 1455 of yacc.c  */
#line 326 "parse.y"
    {v->fname = (yyvsp[(3) - (5)].sval); v->fname2 = (yyvsp[(5) - (5)].sval); }
    break;

  case 79:

/* Line 1455 of yacc.c  */
#line 328 "parse.y"
    { 
			/* this will eat the ';' as well, but we're bailing out anyway: */
			YYERROR; 
		}
    break;

  case 80:

/* Line 1455 of yacc.c  */
#line 334 "parse.y"
    { 
			/* only allow this when called through read_variogram(): */
			assert(parse_variogram != NULL);
			v = parse_variogram;
			v->n_models = v->n_fit = 0;
		}
    break;

  case 81:

/* Line 1455 of yacc.c  */
#line 340 "parse.y"
    {
			id = which_identifier((yyvsp[(3) - (4)].sval));
			v = get_vgm(LTI(id,id));
			v->id = v->id1 = v->id2 = id;
			v->n_models = v->n_fit = 0;
		}
    break;

  case 82:

/* Line 1455 of yacc.c  */
#line 346 "parse.y"
    {
			id1 = which_identifier((yyvsp[(3) - (6)].sval));
			id2 = which_identifier((yyvsp[(5) - (6)].sval));
			id = LTI(id1,id2);
			v = get_vgm(id);
			v->id = id;
			v->id1 = id1;
			v->id2 = id2;
			v->n_models = v->n_fit = 0;
		}
    break;

  case 85:

/* Line 1455 of yacc.c  */
#line 362 "parse.y"
    {
			range[0] = 0.0;
			push_to_v(v, (yyvsp[(2) - (4)].sval), (yyvsp[(1) - (4)].dval), range, 1, NULL, fit_sill, fit_range);
		}
    break;

  case 86:

/* Line 1455 of yacc.c  */
#line 366 "parse.y"
    {
			push_to_v(v, (const char *) (yyvsp[(2) - (5)].sval), (yyvsp[(1) - (5)].dval), range, nrangepars, anis, fit_sill, fit_range);
		}
    break;

  case 87:

/* Line 1455 of yacc.c  */
#line 369 "parse.y"
    {
			push_to_v(v, (yyvsp[(3) - (6)].sval), (yyvsp[(2) - (6)].dval), range, nrangepars, anis, fit_sill, fit_range);
		}
    break;

  case 88:

/* Line 1455 of yacc.c  */
#line 372 "parse.y"
    {
			push_to_v(v, (yyvsp[(3) - (6)].sval), -1.0 * (yyvsp[(2) - (6)].dval), range, nrangepars, anis, 
				fit_sill, fit_range);
		}
    break;

  case 89:

/* Line 1455 of yacc.c  */
#line 378 "parse.y"
    { 
			if (which_variogram_model((yyvsp[(1) - (1)].sval)) == NOT_SP) {
				lex_error(); YYERROR;
			}
	}
    break;

  case 90:

/* Line 1455 of yacc.c  */
#line 385 "parse.y"
    { 
			range[0] = (yyvsp[(1) - (1)].dval); 
			nrangepars = 1;
			anis[0] = -9999.0; 
		}
    break;

  case 91:

/* Line 1455 of yacc.c  */
#line 390 "parse.y"
    {
			range[0] = (yyvsp[(1) - (3)].dval);
			range[1] = (yyvsp[(3) - (3)].dval);
			nrangepars = 2;
		}
    break;

  case 92:

/* Line 1455 of yacc.c  */
#line 395 "parse.y"
    {
			range[0] = (yyvsp[(1) - (5)].dval);
			nrangepars = 1;
			anis[0] = (yyvsp[(3) - (5)].dval);
			anis[3] = (yyvsp[(5) - (5)].dval);
			anis[1] = anis[2] = 0.0;
			anis[4] = 1.0;
		}
    break;

  case 93:

/* Line 1455 of yacc.c  */
#line 403 "parse.y"
    {
			range[0] = (yyvsp[(1) - (7)].dval);
			range[1] = (yyvsp[(3) - (7)].dval);
			nrangepars = 2;
			anis[0] = (yyvsp[(5) - (7)].dval);
			anis[3] = (yyvsp[(7) - (7)].dval);
			anis[1] = anis[2] = 0.0;
			anis[4] = 1.0;
		}
    break;

  case 94:

/* Line 1455 of yacc.c  */
#line 412 "parse.y"
    {
			range[0] = (yyvsp[(1) - (11)].dval);
			nrangepars = 1;
			anis[0] = (yyvsp[(3) - (11)].dval);
			anis[1] = (yyvsp[(5) - (11)].dval);
			anis[2] = (yyvsp[(7) - (11)].dval);
			anis[3] = (yyvsp[(9) - (11)].dval);
			anis[4] = (yyvsp[(11) - (11)].dval);
		}
    break;

  case 95:

/* Line 1455 of yacc.c  */
#line 421 "parse.y"
    {
			range[0] = (yyvsp[(1) - (13)].dval);
			range[1] = (yyvsp[(3) - (13)].dval);
			nrangepars = 2;
			anis[0] = (yyvsp[(5) - (13)].dval);
			anis[1] = (yyvsp[(7) - (13)].dval);
			anis[2] = (yyvsp[(9) - (13)].dval);
			anis[3] = (yyvsp[(11) - (13)].dval);
			anis[4] = (yyvsp[(13) - (13)].dval);
		}
    break;

  case 96:

/* Line 1455 of yacc.c  */
#line 433 "parse.y"
    { fit_sill = 1; }
    break;

  case 97:

/* Line 1455 of yacc.c  */
#line 434 "parse.y"
    { fit_sill = 0; (yyval.dval) = (yyvsp[(2) - (2)].dval); }
    break;

  case 98:

/* Line 1455 of yacc.c  */
#line 435 "parse.y"
    { fit_sill = 0; (yyval.dval) = (yyvsp[(1) - (2)].dval); }
    break;

  case 99:

/* Line 1455 of yacc.c  */
#line 438 "parse.y"
    { fit_range = 1; }
    break;

  case 100:

/* Line 1455 of yacc.c  */
#line 439 "parse.y"
    { fit_range = 0; (yyval.dval) = (yyvsp[(2) - (2)].dval); }
    break;

  case 101:

/* Line 1455 of yacc.c  */
#line 442 "parse.y"
    { o_filename = (yyvsp[(3) - (3)].sval); }
    break;

  case 102:

/* Line 1455 of yacc.c  */
#line 443 "parse.y"
    { 
			id = which_identifier((yyvsp[(3) - (6)].sval));
			ofn = (char **) get_outfile_name();
			ofn[2 * id] = (yyvsp[(6) - (6)].sval);
		}
    break;

  case 103:

/* Line 1455 of yacc.c  */
#line 448 "parse.y"
    { 
			if (get_n_vars() == 0) {
				lex_error();
				ErrMsg(ER_SYNTAX, "define data first");
			}
			ofn = (char **) get_outfile_name();
			ofn[0] = (yyvsp[(3) - (3)].sval);
		}
    break;

  case 104:

/* Line 1455 of yacc.c  */
#line 456 "parse.y"
    { 
			id = which_identifier((yyvsp[(3) - (6)].sval));
			ofn = (char **) get_outfile_name();
			ofn[2 * id + 1] = (yyvsp[(6) - (6)].sval);
		}
    break;

  case 105:

/* Line 1455 of yacc.c  */
#line 461 "parse.y"
    { 
			if (get_n_vars() == 0) {
				lex_error();
				ErrMsg(ER_SYNTAX, "define data first");
			}
			ofn = (char **) get_outfile_name();
			ofn[1] = (yyvsp[(3) - (3)].sval);
		}
    break;

  case 106:

/* Line 1455 of yacc.c  */
#line 469 "parse.y"
    { 
			id = get_n_vars();
			id1 = which_identifier((yyvsp[(3) - (8)].sval));
			id2 = which_identifier((yyvsp[(5) - (8)].sval));
			if (id != get_n_vars())	
				ErrMsg(ER_SYNTAX, "define all data(..) before covariances(..,..)");
			ofn = (char **) get_outfile_name();
			id = 2 * id + LTI2(id1, id2);
			ofn[id] = (yyvsp[(8) - (8)].sval);
		}
    break;

  case 109:

/* Line 1455 of yacc.c  */
#line 485 "parse.y"
    {
			switch (expr.what) { 
				case IS_INT: *((int *)expr.ptr) = (yyvsp[(3) - (3)].ival); break;
				case IS_REAL: *((double *)expr.ptr) = (double) (yyvsp[(3) - (3)].ival); break;
				default: lex_error(); YYERROR;
			}
			check_assign_expr(&expr);
		}
    break;

  case 110:

/* Line 1455 of yacc.c  */
#line 493 "parse.y"
    {
			switch (expr.what) { 
				case IS_UINT: *((unsigned int *)expr.ptr) = (yyvsp[(3) - (3)].uval); break;
				default: lex_error(); YYERROR;
			}
			check_assign_expr(&expr);
		}
    break;

  case 111:

/* Line 1455 of yacc.c  */
#line 500 "parse.y"
    {
			switch (expr.what) {
				case IS_REAL: *((double *)expr.ptr) = (yyvsp[(3) - (3)].dval); break;
				default: lex_error(); YYERROR; break;
			}
			check_assign_expr(&expr);
		}
    break;

  case 112:

/* Line 1455 of yacc.c  */
#line 507 "parse.y"
    {
			if (expr.what != IS_STRING) {
				lex_error();
				YYERROR;
			}
			*((char **) expr.ptr) = (yyvsp[(3) - (3)].sval);
		}
    break;

  case 113:

/* Line 1455 of yacc.c  */
#line 516 "parse.y"
    { if (! is_set_expr(&expr, (yyvsp[(1) - (1)].sval))) { lex_error(); YYERROR; }}
    break;

  case 114:

/* Line 1455 of yacc.c  */
#line 519 "parse.y"
    {
			for (id = 1; methods[id].name != NULL; id++) {
				if (almost_equals((yyvsp[(3) - (3)].sval), methods[id].name)) {
					set_method(methods[id].m);
					break; /* id-loop */
				}
			}
			if (methods[id].m == NSP) {
				lex_error();
				YYERROR;
			}
		}
    break;

  case 116:

/* Line 1455 of yacc.c  */
#line 536 "parse.y"
    { push_mask_name((yyvsp[(1) - (1)].sval)); }
    break;

  case 117:

/* Line 1455 of yacc.c  */
#line 537 "parse.y"
    { push_mask_name((yyvsp[(3) - (3)].sval)); }
    break;

  case 119:

/* Line 1455 of yacc.c  */
#line 543 "parse.y"
    { push_edges_name((yyvsp[(1) - (1)].sval)); }
    break;

  case 120:

/* Line 1455 of yacc.c  */
#line 544 "parse.y"
    { push_edges_name((yyvsp[(3) - (3)].sval)); }
    break;

  case 121:

/* Line 1455 of yacc.c  */
#line 547 "parse.y"
    {
			if (!almost_equals((yyvsp[(3) - (4)].sval), "w$ith"))
				lex_error();
			id1 = which_identifier((yyvsp[(2) - (4)].sval));
			id2 = which_identifier((yyvsp[(4) - (4)].sval));
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
    break;

  case 122:

/* Line 1455 of yacc.c  */
#line 566 "parse.y"
    {
			if (!almost_equals((yyvsp[(6) - (10)].sval), "w$ith"))
				lex_error();
			id1 = which_identifier((yyvsp[(2) - (10)].sval));
			id2 = which_identifier((yyvsp[(7) - (10)].sval));
			dpp = get_gstat_data();
			if (dpp[id1]->id != id1 || dpp[id2]->id != id2) {
				lex_error();
				ErrMsg(ER_IMPOSVAL, "define data before attempting to merge");
			}
			col1 = (yyvsp[(4) - (10)].ival);
			col2 = (yyvsp[(9) - (10)].ival);
			if (id1 < id2) { /* swap id and col */
				id = id1; id1 = id2; id2 = id;
				id = col1; col1 = col2; col2 = id;
			}
			if (push_to_merge_table(dpp[id1], id2, col1, col2)) {
				lex_error();
				ErrMsg(ER_IMPOSVAL, "attempt to merge failed");
			}
		}
    break;

  case 123:

/* Line 1455 of yacc.c  */
#line 587 "parse.y"
    {
			if (!almost_equals((yyvsp[(6) - (10)].sval), "w$ith"))
				lex_error();
			id1 = (yyvsp[(2) - (10)].ival);
			id2 = (yyvsp[(7) - (10)].ival);
			if (id1 >= get_n_vars() || id2 >= get_n_vars() || id1 < 0 || id2 < 0) {
				lex_error();
				ErrMsg(ER_IMPOSVAL, "id values out of range");
			}
			col1 = (yyvsp[(4) - (10)].ival);
			col2 = (yyvsp[(9) - (10)].ival);
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
    break;

  case 124:

/* Line 1455 of yacc.c  */
#line 610 "parse.y"
    { boundary_file = (yyvsp[(3) - (3)].sval); }
    break;

  case 126:

/* Line 1455 of yacc.c  */
#line 614 "parse.y"
    { push_bound((yyvsp[(1) - (1)].dval)); }
    break;

  case 127:

/* Line 1455 of yacc.c  */
#line 615 "parse.y"
    { push_bound((yyvsp[(2) - (2)].dval)); }
    break;

  case 128:

/* Line 1455 of yacc.c  */
#line 616 "parse.y"
    { push_bound((yyvsp[(3) - (3)].dval)); }
    break;

  case 130:

/* Line 1455 of yacc.c  */
#line 622 "parse.y"
    { push_marginal(NULL, (yyvsp[(1) - (1)].dval)); }
    break;

  case 131:

/* Line 1455 of yacc.c  */
#line 623 "parse.y"
    { push_marginal(NULL, (yyvsp[(3) - (3)].dval)); }
    break;

  case 132:

/* Line 1455 of yacc.c  */
#line 624 "parse.y"
    { push_marginal((yyvsp[(1) - (1)].sval), -1.0); }
    break;

  case 133:

/* Line 1455 of yacc.c  */
#line 625 "parse.y"
    { push_marginal((yyvsp[(3) - (3)].sval), -1.0); }
    break;



/* Line 1455 of yacc.c  */
#line 2687 "y.tab.c"
      default: break;
    }
  YY_SYMBOL_PRINT ("-> $$ =", yyr1[yyn], &yyval, &yyloc);

  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);

  *++yyvsp = yyval;

  /* Now `shift' the result of the reduction.  Determine what state
     that goes to, based on the state we popped back to and the rule
     number reduced by.  */

  yyn = yyr1[yyn];

  yystate = yypgoto[yyn - YYNTOKENS] + *yyssp;
  if (0 <= yystate && yystate <= YYLAST && yycheck[yystate] == *yyssp)
    yystate = yytable[yystate];
  else
    yystate = yydefgoto[yyn - YYNTOKENS];

  goto yynewstate;


/*------------------------------------.
| yyerrlab -- here on detecting error |
`------------------------------------*/
yyerrlab:
  /* If not already recovering from an error, report this error.  */
  if (!yyerrstatus)
    {
      ++yynerrs;
#if ! YYERROR_VERBOSE
      yyerror (YY_("syntax error"));
#else
      {
	YYSIZE_T yysize = yysyntax_error (0, yystate, yychar);
	if (yymsg_alloc < yysize && yymsg_alloc < YYSTACK_ALLOC_MAXIMUM)
	  {
	    YYSIZE_T yyalloc = 2 * yysize;
	    if (! (yysize <= yyalloc && yyalloc <= YYSTACK_ALLOC_MAXIMUM))
	      yyalloc = YYSTACK_ALLOC_MAXIMUM;
	    if (yymsg != yymsgbuf)
	      YYSTACK_FREE (yymsg);
	    yymsg = (char *) YYSTACK_ALLOC (yyalloc);
	    if (yymsg)
	      yymsg_alloc = yyalloc;
	    else
	      {
		yymsg = yymsgbuf;
		yymsg_alloc = sizeof yymsgbuf;
	      }
	  }

	if (0 < yysize && yysize <= yymsg_alloc)
	  {
	    (void) yysyntax_error (yymsg, yystate, yychar);
	    yyerror (yymsg);
	  }
	else
	  {
	    yyerror (YY_("syntax error"));
	    if (yysize != 0)
	      goto yyexhaustedlab;
	  }
      }
#endif
    }



  if (yyerrstatus == 3)
    {
      /* If just tried and failed to reuse lookahead token after an
	 error, discard it.  */

      if (yychar <= YYEOF)
	{
	  /* Return failure if at end of input.  */
	  if (yychar == YYEOF)
	    YYABORT;
	}
      else
	{
	  yydestruct ("Error: discarding",
		      yytoken, &yylval);
	  yychar = YYEMPTY;
	}
    }

  /* Else will try to reuse lookahead token after shifting the error
     token.  */
  goto yyerrlab1;


/*---------------------------------------------------.
| yyerrorlab -- error raised explicitly by YYERROR.  |
`---------------------------------------------------*/
yyerrorlab:

  /* Pacify compilers like GCC when the user code never invokes
     YYERROR and the label yyerrorlab therefore never appears in user
     code.  */
  if (/*CONSTCOND*/ 0)
     goto yyerrorlab;

  /* Do not reclaim the symbols of the rule which action triggered
     this YYERROR.  */
  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);
  yystate = *yyssp;
  goto yyerrlab1;


/*-------------------------------------------------------------.
| yyerrlab1 -- common code for both syntax error and YYERROR.  |
`-------------------------------------------------------------*/
yyerrlab1:
  yyerrstatus = 3;	/* Each real token shifted decrements this.  */

  for (;;)
    {
      yyn = yypact[yystate];
      if (yyn != YYPACT_NINF)
	{
	  yyn += YYTERROR;
	  if (0 <= yyn && yyn <= YYLAST && yycheck[yyn] == YYTERROR)
	    {
	      yyn = yytable[yyn];
	      if (0 < yyn)
		break;
	    }
	}

      /* Pop the current state because it cannot handle the error token.  */
      if (yyssp == yyss)
	YYABORT;


      yydestruct ("Error: popping",
		  yystos[yystate], yyvsp);
      YYPOPSTACK (1);
      yystate = *yyssp;
      YY_STACK_PRINT (yyss, yyssp);
    }

  *++yyvsp = yylval;


  /* Shift the error token.  */
  YY_SYMBOL_PRINT ("Shifting", yystos[yyn], yyvsp, yylsp);

  yystate = yyn;
  goto yynewstate;


/*-------------------------------------.
| yyacceptlab -- YYACCEPT comes here.  |
`-------------------------------------*/
yyacceptlab:
  yyresult = 0;
  goto yyreturn;

/*-----------------------------------.
| yyabortlab -- YYABORT comes here.  |
`-----------------------------------*/
yyabortlab:
  yyresult = 1;
  goto yyreturn;

#if !defined(yyoverflow) || YYERROR_VERBOSE
/*-------------------------------------------------.
| yyexhaustedlab -- memory exhaustion comes here.  |
`-------------------------------------------------*/
yyexhaustedlab:
  yyerror (YY_("memory exhausted"));
  yyresult = 2;
  /* Fall through.  */
#endif

yyreturn:
  if (yychar != YYEMPTY)
     yydestruct ("Cleanup: discarding lookahead",
		 yytoken, &yylval);
  /* Do not reclaim the symbols of the rule which action triggered
     this YYABORT or YYACCEPT.  */
  YYPOPSTACK (yylen);
  YY_STACK_PRINT (yyss, yyssp);
  while (yyssp != yyss)
    {
      yydestruct ("Cleanup: popping",
		  yystos[*yyssp], yyvsp);
      YYPOPSTACK (1);
    }
#ifndef yyoverflow
  if (yyss != yyssa)
    YYSTACK_FREE (yyss);
#endif
#if YYERROR_VERBOSE
  if (yymsg != yymsgbuf)
    YYSTACK_FREE (yymsg);
#endif
  /* Make sure YYID is used.  */
  return YYID (yyresult);
}



/* Line 1675 of yacc.c  */
#line 628 "parse.y"


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

