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
 * writecmd.c: write command file from current settings 
 */
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <limits.h>

#include "defs.h"
#include "debug.h"
#include "version.h"
#include "data.h"
#include "utils.h"
#include "vario.h"
#include "defaults.h"
#include "glvars.h"
#include "random.h"
#include "userio.h"
#include "writecmd.h"

void fprint_cmd(FILE *f) { 
	fprintf(f, "%s", sprint_cmd());
}

void logprint_cmd() { 
	printlog("%s", sprint_cmd());
}

const char *sprint_cmd(void) {
	int i, j;
	DATA **d = NULL, *dv;
	char *cp = NULL;
	const char *ccp, *name;
	char tmp[ERROR_BUFFER_SIZE];
	static char *s = NULL;
	time_t tp;

	s = (char *) erealloc(s, (5 + get_n_vars()) * ERROR_BUFFER_SIZE 
			* sizeof(char));
	s[0] = '\0';

	d = get_gstat_data();
	tp = time(&tp);

	strcat(s, "#\n");
	sprintf(tmp, "# gstat command file, %s version %s\n", GSTAT_OS, VERSION);
	strcat(s, tmp);

	sprintf(tmp, "# %s", asctime(localtime(&tp))); /* includes the \n */
	strcat(s, tmp);
	strcat(s, "#\n");

	/* data */
	for (i = 0; i < get_n_vars(); i++) {
		sprintf(tmp, "data(%s): %s\n", name_identifier(i),
			cp = print_data_line(d[i], &cp));
		strcat(s, tmp);
	}
	strcat(s, "\n");

	/* variograms */
	for (i = 0; i < get_n_vars(); i++)
		for (j = i; j >= 0; j--)
			strcat(s, sprint_variogram(get_vgm(LTI(i,j)), 0));
	strcat(s, "\n");

	/* the method: */
	if ((i = get_method()) != get_default_method()) {
		strcat(s, "method: ");
		ccp = methods[i].name;
		while (*ccp) {
			if (*ccp != '$') {
				sprintf(tmp, "%c", *ccp);
				strcat(s, tmp);
			}
			ccp++;
		}
		strcat(s, ";\n");
	}

	/* write bounds: */
	if (gl_bounds != DEF_bounds) {
		strcat(s, "bounds: ");
		for (i = 0; gl_bounds[i] >= 0.0; i++) {
			sprintf(tmp, "%g", gl_bounds[i]);
			strcat(s, tmp);
			if (gl_bounds[i+1] <= 0.0) /* last one */
				strcat(s, ";\n");
			else if (i > 0 && i % 10 == 0)
				strcat(s, ",\n\t");
			else
				strcat(s, ", ");
		}
	}

	if (get_n_masks()) {
		strcat(s, "mask: ");
		for (i = 0; i < get_n_masks(); i++) {
			sprintf(tmp, "'%s'%s", get_mask_name(i),
				i < get_n_masks()-1 ? ", " : ";\n");
			strcat(s, tmp);
		}
	}

	dv = get_dataval();
	if (dv->id > -1) {
		sprintf(tmp, "data(): %s\n", cp = print_data_line(dv, &cp));
		strcat(s, tmp);
	}

	for (i = 0; i < get_n_vars(); i++) {
		if ((name = get_outfile_namei(2 * i))) {
			sprintf(tmp, "predictions(%s): '%s';\n", name_identifier(i), name);
			strcat(s, tmp);
		}
		if ((name = get_outfile_namei(2 * i + 1))) {
			sprintf(tmp, "variances(%s): '%s';\n", name_identifier(i), name);
			strcat(s, tmp);
		}
		for (j = 0; j < i; j++)
			if ((name = get_outfile_namei(2 * get_n_vars() + LTI2(i,j)))) {
				sprintf(tmp, "covariances(%s, %s): '%s';\n",
					name_identifier(i), name_identifier(j), name);
				strcat(s, tmp);
			}
	}

	d = get_gstat_data();
	for (i = 0; i < get_n_vars(); i++)
		for (j = 0; j < d[i]->n_merge; j++) {
			sprintf(tmp, "merge %s(%d) with %s(%d);\n",
				name_identifier(i), d[i]->mtbl[j].col_this_X,
				name_identifier(d[i]->mtbl[j].to_var), d[i]->mtbl[j].col_other_X);
			strcat(s, tmp);
		}

	strcat(s, sprint_glvars(0));
	if (cp != NULL)
		efree(cp);
	
	return s;
}

const char *sprint_glvars(int anyway) {
	char tmp[ERROR_BUFFER_SIZE] = "";
	static char *s = NULL;

	s = (char *) erealloc(s, 5 * ERROR_BUFFER_SIZE * sizeof(char));
	s[0] = '\0';

#define FWRITE(var, default, fmt) \
	if (var != default) { \
		sprintf(tmp, fmt, var); \
		strcat(s, tmp); \
	} else if (anyway) { \
		strcat(s, "# "); \
		sprintf(tmp, fmt, var); \
		strcat(s, tmp); \
	}
#define FSWRITE(var, default, fmt) \
	if ((var != NULL && default != NULL && var != default && strcmp(var, default) != 0)) { \
		sprintf(tmp, fmt, var == NULL ? "[not set]" : var); \
		strcat(s, tmp); \
	} else if (anyway) { \
		strcat(s, "# "); \
		sprintf(tmp, fmt, var == NULL ? "[not set]" : var); \
		strcat(s, tmp); \
	}

	FWRITE(gl_alpha, DEF_alpha, "set alpha = %g;\n");
	FWRITE(gl_beta, DEF_beta, "set beta = %g;\n");
	FWRITE(gl_cressie, DEF_cressie, "set Cressie = %d;\n");
	FWRITE(gl_cutoff, DEF_cutoff, "set cutoff = %g;\n");
	FWRITE(debug_level, DB_NORMAL, "set debug = %d;\n");
	FWRITE(gl_dots, DEF_dots, "set dots = %d;\n");
	FWRITE(gl_fit, DEF_fit, "set fit = %d;\n");
	FWRITE(gl_fit_limit, DEF_fit_limit, "set fit_limit = %g;\n");
	FSWRITE(gl_format, DEF_format, "set format = '%s';\n");
	FWRITE(gl_fraction, DEF_fraction, "set fraction = %g;\n");
	FWRITE(gl_gcv, DEF_gcv, "set gcv = %g;\n");
	FWRITE(gl_gls_residuals, DEF_gls_residuals, "set gls_residuals = %d;\n");
	FSWRITE(gl_gnuplot, DEF_gnuplot, "set gnuplot = '%s';\n");
	FSWRITE(gl_gnuplot35, NULL, "set gnuplot35 = '%s';\n");
	FSWRITE(gl_gpterm, DEF_gpterm /* NULL */, "set gpterm = '%s';\n");
	FWRITE(gl_cn_max, DEF_cn_max, "set cn_max = %g;\n");
	FWRITE(gl_idp, DEF_idp, "set idp = %g;\n");
	FWRITE(gl_n_intervals, DEF_intervals, "set intervals = %d;\n");
	FWRITE(gl_iter, DEF_iter, "set iter = %d;\n");
	FWRITE(gl_jgraph, DEF_jgraph, "set jgraph = %d;\n");
	FSWRITE(gl_mv_string, DEF_mv_string, "set mv = '%s';\n");
	FWRITE(gl_n_uk, DEF_n_uk, "set n_uk = %d;\n");
	FWRITE(gl_nblockdiscr, DEF_nblockdiscr, "set nblockdiscr = %d;\n");
	FWRITE(gl_nsim, DEF_nsim, "set nsim = %d;\n");
	FWRITE(gl_numbers, DEF_numbers, "set numbers = %d;\n");
	FWRITE(gl_nocheck, DEF_nocheck, "set nocheck = %d;\n");
	FSWRITE(o_filename, DEF_ofilename, "set output = '%s';\n");
	FWRITE(gl_order, DEF_order, "set order = %d;\n");
	FSWRITE(gl_pager, DEF_pager, "set pager = '%s';\n");
	FWRITE(gl_rp, DEF_rp, "set rp = %d;\n");
	FSWRITE(gl_plotfile, DEF_plotfile, "set plotfile = '%s';\n");
	FSWRITE(logfile_name, NULL, "set logfile = '%s';\n");
	FWRITE(gl_quantile, DEF_quantile, "set quantile = %g;\n");
	if (gl_seed > INT_MAX) {
		FWRITE(gl_seed, DEF_seed, "set useed = %uU;\n");
	} else {
		FWRITE(gl_seed, DEF_seed, "set seed = %u;\n");
	}
	FWRITE(gl_sym_ev, DEF_sym_ev, "set sym = %d;\n");
	FWRITE(gl_tol_hor, DEF_tol_hor, "set tol_hor = %g;\n");
	FWRITE(gl_tol_ver, DEF_tol_ver, "set tol_ver = %g;\n");
	FWRITE(gl_iwidth, DEF_iwidth, "set width = %g;\n");
	FWRITE(gl_xvalid, DEF_xvalid, "set xvalid = %d;\n");
	FWRITE(gl_zero_est, DEF_zero_est, "set zero_dist = %d;\n");
	FWRITE(gl_zero, DEF_zero, "set zero = %g;\n");
	if (!is_mv_double(&gl_zmap)) {
		sprintf(tmp, "set zmap = %g;\n", gl_zmap);
		strcat(s, tmp);
	}
	if (gl_seed == DEF_seed && anyway) {
		sprintf(tmp, "# seed used: %lu\n", get_seed());
		strcat(s, tmp);
	}
	return s; 
}
