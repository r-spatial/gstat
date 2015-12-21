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
 * sem_main.c: handle "gstat -e [semivariogram|covariogram] ..."
 * command line interface
 */
#include <stdio.h>
#include <stdlib.h> 
#include <string.h>

#include "defs.h"
#ifdef HAVE_GETOPT_H
# include <getopt.h>
#endif
#ifdef HAVE_UNISTD_H
# include <unistd.h>
#endif
#ifndef HAVE_GETOPT
# include "getopt.h"
#endif

#include "userio.h"
#include "debug.h"
#include "read.h"
#include "lex.h" /* -B option: uses parser */
#include "data.h"
#include "utils.h"
#include "vario.h"
#include "glvars.h"
#include "defaults.h"
#include "direct.h"
#include "sem.h"
#include "sem_main.h"


#define ARGFLT ErrMsg(ER_RDFLT, optarg)
#define ARGINT ErrMsg(ER_RDINT, optarg)

static char *progname = NULL;
static void do_help(void);

#ifndef MAIN
int main_sem(int argc, char *argv[]) 
#else
int main(int argc, char *argv[]) 
#endif
{
	char *ofname = "-", *cp = NULL; 
	int c, var2 = 0, i, j;
	double Jcutoff, iwidth, cutoff;
	DATA **d = NULL, *data1 = NULL, *data2 = NULL;
	VARIOGRAM *vp = NULL;

	set_gstat_log_file(stderr);
	progname = string_dup(argv[0]);
	if (*progname == '-')
		progname++;
	set_mv_double(&Jcutoff);

	if (which_identifier("var1") != 0)
		ErrMsg(ER_IMPOSVAL, "sem_main: Whichsample_var() != 0");
	if (which_identifier("var2") != 1)
		ErrMsg(ER_IMPOSVAL, "sem_main: Whichsample_var() != 1");
	d = get_gstat_data();
	data1 = d[0];
	set_mv_double(&iwidth);
	set_mv_double(&cutoff);
	data1->id = 0;
	data1->colnx = 1;
	data1->colny = 2;
	data1->colnvalue = 3;
	data1->average = 0;
	data1->fname = "-";
	set_mv_double(&(data1->mv));
	if (argc == 1) {
		printlog("usage: gstat -e %s [options] file\n", progname);
		printlog("gstat -e %s -h for help on options\n", progname);
	}
	/* parse options: */ 
	opterr = 0;
	while ((c = getopt(argc, argv, 
			"Aa:B:b:Cc:hI:i:J:lm:o:st:u:v:V:w:X:x:y:z:Z:")) != EOF) {
		switch (c) {
			case 'A': data1->average = 1; break;
			case 'a': if (read_double(optarg, &gl_alpha)) ARGFLT; break;
			case 'b': if (read_double(optarg, &gl_beta)) ARGFLT; break;
			case 'B': parse_file(optarg); break;
			case 'C': gl_cressie = 1; break;
			case 'c': if (read_double(optarg, &cutoff)) ARGFLT; 
				if (cutoff <= 0.0)
					ErrMsg(ER_RANGE, "positive cutoff required");
				break;
			case 'd': if (read_double(optarg, &(data1->dX))) ARGFLT; break;
			case 'h': do_help(); exit(0); 
			case 'I': if (read_double(optarg, &(data1->Icutoff))) ARGFLT; break;
			case 'J': if (read_double(optarg, &Jcutoff)) ARGFLT; break;
			case 'i': if (read_double(optarg, &iwidth)) ARGFLT;
				if (iwidth < 0.0)
					ErrMsg(ER_RANGE, "nonnegative interval width required");
				break;
			case 'l': data1->log = 1; break;
			case 'm': if (read_double(optarg, &(data1->mv))) ARGFLT; break;
			case 'o': ofname = optarg; break;
			case 's': debug_level = 0; break;
			case 't': if (read_double(optarg, &gl_tol_hor)) ARGFLT; break;
			case 'u': if (read_double(optarg, &gl_tol_ver)) ARGFLT; break;
			case 'v': if (read_int(optarg, &(data1->colnvalue))) ARGINT; break;
			case 'V': if (read_int(optarg, &(data1->colnvariance))) ARGINT;
				break;
			case 'w': if (read_int(optarg, &(var2))) ARGINT; break;
			case 'x': if (read_int(optarg, &(data1->colnx))) ARGINT; break;
			case 'X': i = 0; 
				while ((cp = strtok((i == 0) ? optarg : NULL, ",")) != NULL) {
					if (read_int(cp, &j)) 
						ARGINT; 
					if (j == -1) {
						if (data1->n_X != 1)
							ARGINT;
						data1->n_X = 0;
					} else {
						data1->n_X++;
						data1->colX = (int *) 
							erealloc(data1->colX, (data1->n_X) * sizeof(int));
						data1->colX[data1->n_X - 1] = j;
					}
					i++;
				} 
				break;
			case 'y': if (read_int(optarg, &(data1->colny))) ARGINT; break;
			case 'z': if (read_int(optarg, &(data1->colnz))) ARGINT; break; 
			case 'Z': if (read_int(optarg, &gl_zero_est)) ARGINT; break; 
			case '?': do_help(); ErrClo(optopt); break;
		} /* switch */
	} /* while getopt */
/* 
 * if there is one unparsed argument AND it's not "-",
 * it must be the data file name:
 */
	if (optind < argc && ! almost_equals(argv[optind], "-"))
		data1->fname = string_dup(argv[optind]);
	read_gstat_data(data1);
	if (var2 > 0) {
		data2 = d[1];
		data2->id = 1;
		data2->colnx = data1->colnx;
		data2->colny = data1->colny;
		data2->colnz = data1->colnz;
		data2->colnvalue = var2;
		data2->average = 0;
		data2->fname = string_dup(data1->fname);
		/* copy X structure: */
		data2->n_X = data1->n_X;
		data2->colX = (int *) 
			erealloc(data2->colX, data1->n_X * sizeof(int));
		for (i = 0; i < data2->n_X; i++)
			data2->colX[i] = data1->colX[i];
		if (!is_mv_double(&Jcutoff))
			data2->Icutoff = Jcutoff;
		read_gstat_data(data2);
		vp = get_vgm(LTI(0,1));
		vp->id1 = 0;
		vp->id2 = 1;
	} else {
		vp = get_vgm(LTI(0,0));
		vp->id1 = 0;
		vp->id2 = 0;
	}
	report_data(data1);
	if (var2 > 0)
		report_data(data2);
	if (!is_mv_double(&cutoff))
		vp->ev->cutoff = cutoff;
	if (!is_mv_double(&iwidth))
		vp->ev->iwidth = iwidth;
	vp->ev->zero = zero_int2enum(gl_zero_est);
	if (almost_equals(argv[0], "semivariogram")) {
		vp->ev->evt = var2 ? CROSSVARIOGRAM : SEMIVARIOGRAM;
		if (vp->ev->zero == ZERO_DEFAULT)
			vp->ev->zero = ZERO_INCLUDE;
	} else { 
		if (almost_equals(argv[0], "covariogram")) {
			vp->ev->evt = var2 ? CROSSCOVARIOGRAM : COVARIOGRAM;
			if (vp->ev->zero == ZERO_DEFAULT)
				vp->ev->zero = ZERO_SPECIAL;
	 	} else
			ErrMsg(ER_IMPOSVAL, "sem(): unknown mode");
	}
	calc_variogram(vp, ofname);
	free_data(data1);
	if (data2)
		free_data(data2);
	return 0;
}

static void do_help(void) {
printlog(
"usage: gstat -e %s [options] file:\ncalculate sample %s from file\n%s",
progname, *progname == 'c' ? "covariogram" : "semivariogram",
"options: [default value]\n\
-c $ cutoff value: max. dist to calculate semivariances [max_dim/3]\n\
-i $ width of interval for semivariance calculation [cutoff/15]\n\
-o f output to file f [stdout]\n\
-x # -y # -z # x-, y- and z-coordinate column number [1, 2, 0]\n\
-v # attribute column number [3];\n\
-w # second attribute column number for cross (co)variogram [0]\n\
-m $ missing value flag [NA]\n\
-l   take log-values of attributes\n\
-I $ indicator cutoff value (if v <= value then 1 else 0) [none]\n\
-J $ indicator cutoff value second attribute [none]\n\
-a $ direction in <x,y> plane, in positive degrees clockwise from y [0.0]\n\
-b $ direction in z, in positive degrees up from <x,y> plane [0.0]\n\
-t $ horizontal tolerance angle in degrees [90.0]\n\
-u $ vertical tolerance angle in degrees [90.0]\n\
-Z # pairs at lag zero: 1 in first interval, 2 omit, 3 special [0: choose]\n\
-A   average all records with equal coordinates\n\
-B f parse file f (with variogram interval boundaries using bounds command)\n\
-C   use Cressie's robust variogram estimator (2.4.12)\n\
-X $,$,$,... enter X-columns;  -d $ max. dX value for including pairs\n\
     option arguments: $ = value  # = number  f = file \n");
}
