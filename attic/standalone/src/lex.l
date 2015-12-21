%{
  /* NOTE CEES REMOVED section from distro */
/*
    Gstat, a program for geostatistical modelling, prediction and simulation
    Copyright 1992, 2003 (C) Edzer J. Pebesma

    Edzer J. Pebesma, e.pebesma@geog.uu.nl
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
 * lex.l: lexical analyzer, to create lex.c using flex or lex.
 * provides yylex(), a function that splits the input into lexical tokens.
 * note that white space is analyzed, but skipped by the yylex()
 * function by recursing into yylex(): yylex() will never return 
 * a token denoting white space.
 */

#include <stdio.h>
#include <string.h>
#include "defs.h"
#include "userio.h"
#include "utils.h"
#include "read.h"
#include "parse.h" /* YYSTYPE and token values */
#include "lex.h"

#ifdef HAVE_UNISTD_H
# include <unistd.h> /* isatty() */
#endif

#define yylval gstat_yylval

/* 
 * local variables: 
 */
static int yy_Lexpos = 0, yy_Posafternl = 1, My_yy_lineno = 1,
	new_source = 0;
static const char *yy_Source = NULL, *lastnl = NULL;
static const char *yy_File_name = NULL;
static char *cp = NULL;

#ifdef yywrap
#undef yywrap
#endif

int MY_input(void);
int MY_output(int c);
int MY_unput(int c);
int yywrap(void);

static char *unquote(const char *txt);
static char *bquote(char *command);
static void reset_lex(void);

#ifdef FLEX_SCANNER /* copied this from the flex man pages: */
#	undef YY_INPUT
#	define YY_INPUT(buf,result,max_size) { \
		int c = MY_input(); \
		result = (c == '\0') ? YY_NULL : (buf[0] = c, 1); \
	}
#else /* other, non-flex lex: undef i/o macro's */
# ifdef input
#  undef input
# endif
# ifdef unput
#  undef unput 
# endif
# ifdef output
#  undef output
# endif /* and re-define them: */
# define input()    MY_input()
# define unput(c)   MY_unput(c)
# define output(c)  MY_output(c)
#endif /* else FLEX_SCANNER */


#ifdef LDEBUG
#	define SCANNED(a) { if (*a) printlog("scanned %s [%s]\n", a, (char *) yytext); \
	else printlog("scanned %s\n", (char *) yytext); }
#else
#	define SCANNED(a)
#endif

%}


D				[0-9]
INT				[+-]?{D}+
UINT			{D}+[U]
EXP				[eE][-+]?{D}+
REAL1			[-+]?{D}+"."{D}*({EXP})?
REAL2			[-+]?{D}*"."{D}+({EXP})?
REAL3			[-+]?{D}+{EXP}
CHAR			[A-Za-z]
IDENT			{CHAR}({D}|{CHAR}|[][~!@$%^_.\/-])*

DQSTRING		[\"][^\"\n]+[\"]
SQSTRING		['][^'\n]+[']
BQSTRING		[`][^`\n]+[`]

COMMENT			"#"[^\n]*"\n" 
WS				[ \t\032\015]
NL				"\n"

	/*
	 * following is the RULES section:
	 */
%%
{INT}       { 	SCANNED("INT");
				if (read_int((char *) yytext, &(yylval.ival))) {
					lex_error();
					ErrMsg(ER_RDINT, (const char *) yytext); 
				}
				return(INT);
			}
{UINT}      { 	SCANNED("UINT");
				cp = strchr((char *) yytext, 'U');
				*cp = '\0';
				if (read_uint((char *) yytext, &(yylval.uval))) {
					lex_error();
					ErrMsg(ER_RDINT, (const char *) yytext); 
				}
				return(UINT);
			}
{REAL1}     |
{REAL2}     |
{REAL3}     {	SCANNED("REAL");
				if (read_double((char *) yytext, &(yylval.dval))) {
					lex_error();
					ErrMsg(ER_RDFLT, (const char *) yytext); 
				}
				return(REAL);
			}
{DQSTRING}  |
{SQSTRING}  {	SCANNED("[S|D]QSTRING");
				yylval.sval = unquote((char *) yytext);
				return(QSTR);
			}
{BQSTRING}  {	SCANNED("BQSTRING");
				yylval.sval = bquote(unquote((char *) yytext));
				return(QSTR);
			}
{IDENT}     {	/* SCANNED("IDENT"); */
#ifdef USING_R /* don't leak memory */
				yylval.sval = (char *) yytext;
#else /* avoid trouble if we don't duplicate the string: */
				yylval.sval = string_dup((char *) yytext);
#endif
				if (almost_equals((char *) yytext, "data") ||
						almost_equals((char *) yytext, "points")) {
					SCANNED("data"); return(ID_DATA);
				} else if (almost_equals((char *) yytext, "vario$gram"))
					return(ID_VARIOGRAM);
				else if (almost_equals((char *) yytext, "pred$ictions"))
					return(ID_PREDICTIONS);
				else if (almost_equals((char *) yytext, "varia$nces"))
					return(ID_VARIANCES);
				else if (almost_equals((char *) yytext, "cov$ariances"))
					return(ID_COVARIANCES);
				else if (almost_equals((char *) yytext, "mas$ks"))
					return(ID_MASKS);
				else if (almost_equals((char *) yytext, "edges"))
					return(ID_EDGES);
				else if (almost_equals((char *) yytext, "set"))
					return(ID_SET);
				else if (almost_equals((char *) yytext, "mer$ge"))
					return(ID_MERGE);
				else if (almost_equals((char *) yytext, "met$hod"))
					return(ID_METHOD);
				else if (almost_equals((char *) yytext, "area"))
					return(ID_AREA);
				else if (almost_equals((char *) yytext, "bl$ocksize"))
					return(ID_BLOCK);
				else if (almost_equals((char *) yytext, "bo$unds"))
					return(ID_BOUNDS);
				else if (almost_equals((char *) yytext, "mar$ginals"))
					return(ID_MARGINALS);
				else if (almost_equals((char *) yytext, "X"))
					return(ID_X);
				else /* no match */
					return(IDENT);
			}
{WS}        {	SCANNED("WS"); return(yylex()); }
{NL}        {	SCANNED("NL"); My_yy_lineno++; return(yylex()); }
{COMMENT}   {	SCANNED("COMMENT (#)"); My_yy_lineno++; return(yylex()); }
.           {	SCANNED(""); return(*yytext); }

%%
/* User supplied subroutines: */
int MY_input(void) {
	static int add_at_next_call = 0;
	int c;
	
	if (yy_Source == NULL)
		ErrMsg(ER_NULL, "MY_input()");

	if (new_source) {
		yy_Lexpos = 0;
		lastnl = yy_Source;
		new_source = 0;
	}
	if (add_at_next_call == 1) {
		lastnl = yy_Source;
		yy_Posafternl = 1;
		add_at_next_call = 0;
	}
	c = *yy_Source;
	if (c == '\0')
		return 0;
	if (c == '\n') 
		add_at_next_call = 1;
	else
		yy_Posafternl++;
	yy_Source++;
	yy_Lexpos++;
	return c;
}

int MY_output(int c) {
	printlog("%c", c);
	return 0;
}

int MY_unput(int c) {
	if (yy_Source == NULL)
		ErrMsg(ER_NULL, "MY_unput()");
	if (yy_Lexpos == 0)
		ErrMsg(ER_IMPOSVAL, "unput(): probably empty source file");
	if (c != 0) {
		yy_Source--; 
		yy_Lexpos--;
		yy_Posafternl--;
	}
	*((char *) yy_Source) = c;
	return 0;
}

int yywrap(void) {
/* clean up: */
	reset_lex();
	yy_Source = NULL;
	yy_File_name = NULL;
/* and tell yylex() there's no input left: */
	return 1;
}

#define ADD_STR(fmt, value) {sprintf(tmp, fmt, value);strcat(lexerr, tmp);}
#define MAXCHAR 72

void lex_error(void) {
/* 
 * make error position indication message like:
 * points(zn): xx=1, y=2, v=3;
 *             ^^
 */
	int i = 0, j;
	char lexerr[ERROR_BUFFER_SIZE] = "", tmp[ERROR_BUFFER_SIZE] = "";

	/* print current line, \t -> ' ' */
	if (lastnl != NULL && strlen(lastnl) > 0) {
		while (lastnl[i] != '\n' && lastnl[i] != '\0' && i < MAXCHAR) {
			ADD_STR("%c", lastnl[i] == '\t' ? ' ' : lastnl[i]);
			i++;
		}
		ADD_STR("%s", "\n");
		for (j = 0; j < yy_Posafternl - yyleng - 2; j++)
			ADD_STR("%s", " ");
#ifndef FLEX_SCANNER
		if (j > 0)
			ADD_STR("%s", " ");
#endif
		for (j = 0; j < MIN(yyleng, MAXCHAR); j++)
			ADD_STR("%s", "^");
		ADD_STR("%s", "\n");
	}
	ADD_STR("gstat:%s:", yy_File_name ? yy_File_name : "?");
	ADD_STR("%d: ", My_yy_lineno); 
	switch (*yytext) {
		case '\0':
			ADD_STR("%s", "unexcpected end of file\n");
			break;
		case '\n':
			ADD_STR("%s", "unexpected end of line\n");
			break;
		default:
			ADD_STR("`%s' unexpected\n", (char *) yytext);
	}
	message("%s", lexerr);
	reset_lex();
	return;
}

void set_lex_source(const char *source, const char *fname) {
/*
 * initialize sources:
 */
	reset_lex();
	yy_Source = source;
	yy_File_name = fname;
}

static void reset_lex(void) {
	yy_Lexpos = 0;
	new_source = 1;
	yy_Posafternl = 1;
	My_yy_lineno = 1;
}

static char *unquote(const char *txt) {
/* return duplicate of txt, without the quotes (first & last char) */
	int len;
	char *cp;

	len = strlen(txt);
	cp = (char *) emalloc ((len - 1) * sizeof (char));
	strncpy(cp, (const char *) &(txt[1]), len - 2);
	cp[len-2] = '\0';
	return cp;
}

static char *bquote(char *command) {
/*
 * substitute back-quoted string with result from shell command
 * uses pipe or temporary file 
 */
	int size = 0; 
	char *cp = NULL;
	FILE *f;

#ifndef HAVE_POPEN
	char *fn;

	fn = tmpnam(NULL);
	cp = (char *) emalloc((strlen(command) + strlen(fn) + 5) * sizeof(char));
	sprintf(cp, "%s > %s", command, fn);
	esystem(cp);
	size = file_size(fn);
	cp = (char *) erealloc(cp, (size + 1) * sizeof(char));
	f = efopen(fn, "r");
	if (fread(cp, 1, size, f) != size)
		ErrMsg(ER_READ, fn);
	efclose(f);
	eremove(fn);
#else
#define STEP 1000
	int where = 0;

	f = epopen(command, "r");
	do {
		cp = (char *) erealloc(cp, where + STEP);
		where += fread(&(cp[where]), 1, STEP-1, f);
	} while (! feof(f));
	epclose(f);
	cp[size = where] = '\0';
#endif
	/* clean up: */
	if (cp[size-1] == '\n')
		cp[size-1] = '\0';
	else
		cp[size] = '\0';
	efree(command);
	return cp;
}
