
/* A Bison parser, made by GNU Bison 2.4.1.  */

/* Skeleton interface for Bison's Yacc-like parsers in C
   
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

/* Line 1676 of yacc.c  */
#line 118 "parse.y"

	int ival;
	unsigned int uval;
	double dval;
	char *sval;



/* Line 1676 of yacc.c  */
#line 107 "y.tab.h"
} YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
#endif

extern YYSTYPE gstat_yylval;


