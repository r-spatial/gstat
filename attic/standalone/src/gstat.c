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
 * gstat.c: program flow, init/exit procedures
 */
#include <stdio.h>
#include <stdlib.h> 
#include <string.h> /* memcpy */
#include <math.h>

#include "defs.h"

#ifdef HAVE_GETOPT_H
# include <getopt.h>
#endif
#ifdef HAVE_UNISTD_H
# include <unistd.h> /* isatty() */
#endif
#ifndef HAVE_GETOPT
# include "getopt.h"
#endif

#ifdef HAVE_LIBCSF
# include "csf.h"
#endif

#include "userio.h"
#include "debug.h"
#include "data.h"
#include "utils.h"
#include "vario.h"
#include "glvars.h"
#include "xvalid.h"
#include "read.h"
#include "lex.h"
#include "sem_main.h"
#include "sem.h"
#include "reml.h"
#include "fit.h"
#include "plot.h"
#include "stat.h"
#include "random.h"
#include "ui.h"
#include "predict.h"
#include "mapio.h"
#include "maputils.h"
#include "polygon.h"
#include "map2fig.h"
#include "map2gd.h"
#include "palet.h"
#include "sample.h"
#include "ossfim.h"
#include "version.h"

/*
 * gstat, a program for modelling, prediction and simulation of spatial data
 * the program is command file driven, it makes use of two PD libraries:
 * the meschach matrix library (available from netlib sites), and the csf
 * binary raster map library (available from http://www.geog.uu.nl/pcraster/).
 * input: ascii (table, GeoEAS) files or raster maps (csf, idrisi, arcinfo);
 * output: ascii (GeoEAS) files or maps (input map format);
 *
 * for prediction, the main modules are in 
 *   predict.c (loops over prediction locations)
 *   select.c (neighbourhood selection) and 
 *   getest.c (calls the estimation functions)
 *   gls.c (kriging prediction, estimation)
 *   sim.c (additional simulation routines)
 *   mapio.c (map i/o)
 * for sample (co)variogram calculation: sem(_main).c; 
 * for variogram models: vario.h (all structures)
 *   vario.c (general routines)
 *   vario_io.c (point-to block, block-to-block)
 *   vario_rd.c (parse a model from a string)
 *   vario_fn.c (functions for basic models)
 * the rest is organization.
 */

static void parse_options(int argc, char *argv[]);
static void read_all_data(DATA **data, DATA *valdata, int n_vars);
static void do_variogram(int nvars, METHOD m);
static int exec_action(int argc, char **argv);
static int calc_stats(int argc, char *argv[]);

void close_gstat_log_file(void);

#ifdef CYGWIN_ERRNOBUG /* seems to be a bug with -mno-cygwin; cygwin b20 */
int __errno;
#endif

static int check_only = 0;

#define HELP "options:\n\
-C | --copyright    print copyright notice\n\
-W | --warranty     print no warranty claim\n\
-i | --interactive  start interactive mode\n\
-v | --version      print version information\n\
-c | --check        check command file syntax and exit\n\
-e a | --execute a  execute action a\n\n\
the following options have an equivalent command -> syntax:\n\
-d n | --debug=n    debug (sum) level: -1 list options -> set debug = n;\n\
-o f | --output=f   write output results to file f     -> set output = 'f';\n\
-l f | --logfile=f  write debug log to file f          -> set logfile = 'f';\n\
-p f | --plotfile=f write plot commands to file f      -> set plotfile = 'f';\n\
-S   | --secure     secure mode                        -> set secure = 1;\n\
                    (this implies no system(), popen() or remove() calls)\n\
-x   | --xvalid     cross validation                   -> set xvalid = 1;\n\
option arguments: f file name, n number, a action\n"

#ifndef LIBGSTAT
int main(int argc, char *argv[]) {
#else
int gstat_main(int argc, char *argv[]) {
#endif

	DATA      **data = NULL, *valdata = NULL;

/*
 * initialise some global variables:
 */
	atexit(close_gstat_log_file);
	init_userio(1);
	init_global_variables();
	argv0 = argv[0];
#ifdef HAVE_LIBGDAL
	GDALAllRegister();
#endif
/*
 * register command line arguments on command_line:
 */
	command_line = store_argv(argc, argv);
	parse_gstatrc();
/*
 * INPUT: command line options;
 */
	parse_options(argc, argv); /* exits on -e options */

/*
 * start with program heading:
 */
	printlog("%s: %s %s version %s\n", GSTAT_NAME, GSTAT_OS, TARGET, VERSION);
	printlog("%s\n", GSTAT_CR);
	gstat_start();

/*
 * INPUT: Parse command files: 
 */
	if (optind == argc) { /* there's no command file name left */
		if (get_method() != UIF) { /* -i or -m command line option */
			/* no arguments */
			printlog("Updates, manuals, GPL sources: %s\n", 
				GSTAT_HOME);
			printlog("%s\n", USAGE);
			ErrMsg(ER_NOCMD, "");
		} else {
			start_ui();
			exit(0);
		}
	} else { /* we have one or more command files to be processed */
		for ( ; optind < argc; optind++) {
			command_file_name = argv[optind];
			parse_file(command_file_name);
			if (logfile_name != NULL)
				set_gstat_log_file(efopen(logfile_name, "w"));
/* 
 * get global variables locally: 
 */
			data = 			get_gstat_data();
			valdata = 		get_dataval();
			set_seed(gl_seed);

/* 
 * check variable settings and next
 * INPUT: read data values from file: 
 */
			read_all_data(data, valdata, get_n_vars());
			if (get_method() == NSP) /* Still no answer to this: choose default */
				set_method(get_default_method());
			set_mode();
			check_global_variables();
			setup_meschach_error_handler(0);
			if (DEBUG_DUMP)
				dump_all();
			if (get_method() != NSP)
				printlog("[%s]\n", method_string(get_method()));
			if (check_only)
				set_method(NSP);

/*
 * start calculations && OUTPUT routines:
 */
			switch (get_method()) {
				case UIF:
					start_ui();
					break;
				case NSP:
					break;
				case COV: 
				case SEM:
					do_variogram(get_n_vars(), get_method());
					break;
        		case POLY:
            		setup_poly_method();
            		/*FALLTHROUGH*/
				default: 
					if (gl_xvalid) /* validation/cross validation */
						cross_valid(data);
					else
						predict_all(data); /* or prediction: */
					break;
			} /* switch get_method() */
			remove_all(); /* free all data etc. */
			init_global_variables(); /* re-init for next round */
		}
	}

	if (DEBUG_DUMP)
		atexit(print_file_record);

	if (get_method() != UIF)
		elapsed();
/* 
 * file closing & data freeing time:
 */
	if (plotfile != NULL)
		efclose(plotfile);

	exit(0);
} /* end of main() */

static void parse_options(int argc, char *argv[]) {
	int c;
#ifdef HAVE_GETOPT_LONG
	static struct option long_options[] = {
		{ "copyright",  0, 0, 'C' },
		{ "check",      0, 0, 'c' },
		{ "debug",      1, 0, 'd' },
		{ "help",       0, 0, 'h' },
		{ "interactive",0, 0, 'i' },
		{ "logfile",    1, 0, 'l' },
		{ "plotfile",   1, 0, 'p' },
		{ "output",     1, 0, 'o' },
		{ "silent",     0, 0, 's' },
		{ "secure",     0, 0, 'S' },
		{ "version",    0, 0, 'v' },
		{ "warranty",   0, 0, 'W' },
		{ "xvalid",     0, 0, 'x' },
		{ NULL,         0, 0, 0 }
	};
	int option_index = 0;
#endif

	/* who am i, some -e option? quick to getopt() of the -e option */
	if (argc > 1 && (almost_equals(argv[1], "-e$xecute") ||
			almost_equals(argv[1], "--execute")))
		exit(exec_action(argc - 2, argv + 2));

	/* no, we're plain gstat. Parse the command line: */
	opterr = 0;
#ifdef HAVE_GETOPT_LONG
	while ((c = getopt_long(argc, argv, "Ccd:ehil:mo:p:sSvWxV",
					long_options, &option_index)) != EOF)
#else
	while ((c = getopt(argc, argv, "Ccd:ehil:mo:p:sSvWxV")) != EOF)
#endif
	{
		switch (c) {
			case 'd':
				if (read_int(optarg, &debug_level) || debug_level < 0) {
					debug_level = 1;
					message(DEBUG_OPTIONS);
					ErrMsg(ER_ARGOPT, "d");
				}
				break;
			case 'e':
				ErrMsg(ER_ARGOPT, "option -e only allowed as first option");
			case 'i': 
			case 'm': 
				set_method(UIF); 
				break;
			case 'l':
				set_gstat_log_file(efopen(logfile_name = optarg, "w"));
				break;
			case 'o': o_filename = optarg; break;
			case 'p': gl_plotfile = optarg; break;
			case 'S': gl_secure = 1; break;
			case 's': debug_level = 0; break;
			case 'c': check_only = 1; break;
			case 'x': gl_xvalid = 1; break;
			case 'C': printlog("%s\n", COPYRIGHT); break;
			case 'W': printlog("\n%s\n\n", NOWARRANTY); break;
			case 'V': case 'v':
				printlog("compiled on:          %s\n", __DATE__);
				printlog("with libraries:       ");
#ifdef HAVE_LIBCSF
				printlog("csf ");
#endif
#ifdef HAVE_LIBCURSES
				printlog("curses ");
#endif
#ifdef HAVE_LIBGD
				printlog("gd ");
# ifdef HAVE_GDIMAGEGIF
				printlog("(gif) ");
# else
				printlog("(png) ");
# endif
#endif
#ifdef HAVE_LIBGIS
				printlog("grass ");
#endif
#ifdef HAVE_LIBGDAL
				printlog("gdal ");
#endif
#ifdef HAVE_LIBGSL
				printlog("gsl ");
#endif
#ifdef HAVE_LIBNCURSES
				printlog("ncurses ");
#endif
#ifdef HAVE_LIBNETCDF
				printlog("netcdf ");
#endif
#ifdef LIBGSTAT
				printlog("qt ");
#endif
#ifdef HAVE_LIBREADLINE
				printlog("readline ");
#endif
				printlog("\n");
				printlog("last modified on:     %s\n", LASTMOD);
				printlog("gstat home page:      %s\n", GSTAT_HOME);
				printlog("questions, bugs etc.  mailto:%s\n\n", GSTAT_INFO);
				break;
			case 'h':
				printlog("%s\n\n%s\n", USAGE, HELP);
				break;
			case '?':
			default:
				ErrClo(optopt);
		}
	}
	return;
}

static void read_all_data(DATA **data, DATA *valdata, int n_vars) {
	int i;
	DATA *area;

	init_data_minmax();
	area = get_data_area();
	for (i = 0; i < n_vars; i++)  {
		if (get_mode() == STRATIFY)
			printlog("stratum # %d:\n", i + strata_min);
		printlog("data(%s): ", name_identifier(i));
		if (data[i]->id < 0) {
			message("data(%s) was not specified\n", name_identifier(i));
			ErrMsg(ER_SYNTAX, "data specification error");
		}
		read_gstat_data(data[i]);
		report_data(data[i]);
	} /* for i */

/*
 * what to do when area is specified, but no masks or data()?
 * default prediction to `area'. Create a valdata with one point at
 * centre of area (for select()); and centre area at (0,0,0)
 */
	if (area && get_n_masks() <= 0 && valdata->id == -1) {
		valdata->id = ID_OF_VALDATA;
		valdata->centre = area->centre = 1;
	}
/* 
 * read data() data:
 */
	if (valdata->id > -1) {
		setup_valdata_X(valdata);
		if (! valdata->centre)
			valdata = read_gstat_data(valdata);
	}
/*
 * read area, if existed
 */
	if (area != NULL && get_method() != POLY) {
		read_gstat_data(area);
		/* now, before centring area: */
		if (valdata->centre)
			valdata = get_area_centre(area, valdata);
		if (area->centre)
			centre_area(area);
		printlog("area:%s\n", area->centre ? " (centred around 0)" : "");
		report_data(area);
		if (DEBUG_DATA) 
			print_data_list(area);
	}
/*
 * read edges, if existed
 */
    if (get_n_edges() > 0) {
        read_edges();
        report_edges();
/*         setup_visibility_graph(); */
        
        /*setup_planar_subdivisions();*/
    }
/*
 * setup and report data
 */

	if (valdata->id > -1) {
		printlog("data():%s ",
			valdata->centre ? " [at area centre]" : "");
		report_data(valdata);
	}
	for (i = 0; i < n_vars; i++) 
		setup_data_minmax(data[i]);
	if (valdata->id > -1)
		setup_data_minmax(valdata);
	for (i = 0; i < n_vars; i++) 
		calc_polynomials(data[i]);
	if (valdata->id > -1)
		calc_polynomials(valdata);
	if (DEBUG_DATA) {
		for (i = 0; i < n_vars; i++) 
			print_data_list(data[i]);
		if (valdata->id > -1)
			print_data_list(valdata);
	}
}

void close_gstat_log_file(void) {
	if (logfile_name != NULL)
		set_gstat_log_file(NULL); /* closes file */
}

static void do_variogram(int nvars, METHOD m) {
	int i, j;
	VARIOGRAM *vp = NULL;

	if (nvars == 0)
		return;

	for (i = 0; i < nvars; i++) {
		for (j = i; j >= 0; j--) {
			vp = get_vgm(LTI(i,j)); /* */
			vp->id1 = j;
			vp->id2 = i;
			if (m == COV)
				vp->ev->evt = (i != j) ? CROSSCOVARIOGRAM : COVARIOGRAM;
			else
				vp->ev->evt = (i != j) ? CROSSVARIOGRAM : SEMIVARIOGRAM;
			if (vp->fname != NULL || o_filename != NULL) {
				calc_variogram(vp, vp->fname ? vp->fname : o_filename);
				if (vp->n_models > 0 && gl_fit) {
					vp->ev->fit = fit_int2enum(gl_fit);
					if (fit_variogram(vp))
						pr_warning("error during variogram fit");
					else
						logprint_variogram(vp, 1);
				}
			}
		}
	}

	if (plotfile) {
		if (nvars > 1)
			ErrMsg(ER_IMPOSVAL, "plot file only works for single variable");
		if (vp->ev->map)
			ErrMsg(ER_IMPOSVAL, "cannot make plot file for variogram map");
		if (gl_jgraph)
			fprint_jgraph_variogram(plotfile, vp);
		else 
			fprint_gnuplot_variogram(plotfile, vp, "gnuplot.out", GNUPLOT, 0);
	}
}

static int exec_action(int argc, char *argv[]) {
	int i;
	const struct {
		const char *name, *ab_name;
		int (*prog_fn)(int c, char *v[]);
		const char *desc;
	} exec_actions[] = {
		{ "convert",       "con$vert",      map_convert, "convert map format" },
		{ "cover",         "cov$er",        map_cover,   "combine non-missing valued areas from maps" },
		{ "cut",           "cu$t",          map_cut,     "cut square area from map" },
		{ "lnh",           "lnh",           map_lnh,     "" },
		{ "mapdiff",       "mapdiff",       map_diff,    "check difference between two maps" },
		{ "map2png",       "map2png",       map2gd,      "create png image from map" },
		{ "map2gif",       "map2gif",       map2gd,      "create gif image from map" },
		{ "map2fig",       "map2fig",       map2fig,     "create fig (xfig/fig2dev) from map or data" },
		{ "nominal",       "no$minal",      map_nominal, "convert n 0-1 maps into one 0...n-1 map" },
		{ "ossfim",        "o$ssfim",       ossfim,       "kriging errors as function of sample spacing and block size" },
		{ "palet",         "pa$let",        palet,       "print out a colour palet" },
		{ "semivariance",  "semivariance",  vario,       "semivariance for variogram model" },
		{ "covariance",    "covariance",    vario,       "covariance for variogram model" },
		{ "semivariogram", "semivariogram", main_sem,    "sample variogram from data" },
		{ "covariogram",   "covariogram",   main_sem,    "sample covariogram from data" },
		{ "statistics",    "st$atistics",   calc_stats,  "sample summary statistics" },
		{ "sample",        "sa$mple",       sample_main, "map sampling strategy realisations" },
		{ "random",        "ra$ndom",       e_random,    "print random numbers" },
		{ "q",             "q",             map_q,       "" },
		{ NULL, NULL, NULL, NULL }
	};

	for (i = 0; exec_actions[i].name; i++)
		if (almost_equals(argv[0], exec_actions[i].ab_name))
			return exec_actions[i].prog_fn(argc, argv);

	/* not returned -- on error: */
	if (argc)
		printlog("error: unknown action: %s\n", argv[0]);
	printlog("gstat -execute [action] [arguments]\n");
	printlog("valid actions are:\n");
	for (i = 0; exec_actions[i].name; i++)
		if (exec_actions[i].desc[0] != '\0')
			printlog("\t%-16s[%s]\n", exec_actions[i].name, exec_actions[i].desc);
	return 1;
}

static int calc_stats(int argc, char *argv[]) {
/*
 * read double from stdin until EOF, print summary statistics
 * print error on non-white space; gl_mv_string is the missing value
 */
#define CSUSAGE \
"[options] [file [file ...]]\n\
options:\n\
-s don't print header line\n\
-q quantile (0.25)\n"

	int i, c, silent = 0;
	double q = 0.25;

	opterr = 0;
	while ((c = getopt(argc, argv, "f:hq:s?")) != EOF) {
		switch (c) {
			case '?': 
			case 'h': 
				printf("%s %s", argv[0], CSUSAGE); 
				ErrClo(optopt); 
				break;
			case 's': silent = 1; break;
			case 'q':
				if (read_double(optarg, &q))
					ErrMsg(ER_ARGOPT, "reading float");
				if (q < 0.0 || q > 1.0)
					ErrMsg(ER_ARGOPT, "q value out of range [0,1]");
				if (q > 0.5)
					q = 1.0 - q;
				break;
			default:
				ErrMsg(ER_ARGOPT, "unknown option");
				break;
		}
	}

	if (optind == argc) { /* no file arguments left: */
		if (isatty(fileno(stdin))) {
			printf("%s %s", argv[0], CSUSAGE);
			exit(0);
		} else
			stats(NULL, silent, q);
	} else for (i = optind; i < argc; i++)
		stats(argv[i], silent || (i != optind), q);
	return 0;
}
