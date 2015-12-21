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
 * ui.c: variogram modelling user interface (uses curses)
 */
#include <stdio.h> /* printf(), fileno() */
#include <stdlib.h> /* getenv() etc. */
#include <math.h> /* floor(), sin(), cos() */
#include <string.h> /* strchr() */
#include "defs.h" /* may define HAVE_UNISTD_H */
#ifdef HAVE_UNISTD_H
# include <unistd.h> /* isatty() */
#endif

#include "userio.h"

#if !defined(HAVE_LIBCURSES) && !defined(HAVE_LIBNCURSES)
int start_ui(void) {
	ErrMsg(ER_NOCURSES, "curses");
	return 1;
}
#else

#include "debug.h"
#include "data.h"
#include "vario.h"
#include "glvars.h"
#include "read.h"
#include "sem.h"
#include "fit.h"
#include "plot.h"
#include "defaults.h"
#include "direct.h"
#include "version.h"
#include "gls.h"
#include "writecmd.h"
#include "xvalid.h"
#include "predict.h"
#include "lex.h"
#include "ui.h"

#include "curses.h"

#include "utils.h" /* MIN()/MAX() may be #defined through curses.h */

#define LENGTHOF(a) (sizeof(a)/sizeof(char *))
#define SHORTHELP0 \
"Use arrow keys to move, <return> to choose, press `Q' to quit"
#define SHORTHELP1 \
"Press `?' for help on current field, `H' for other commands"

int data_not_valid(DATA *d);
void fix_toggle(VARIOGRAM *v, int fix);
void fit_method(int key);

char *fit_type_str[] = { 
	"no fit", 
	"WLS, weights n(h)", 
	"WLS, weights n(h)/mod(h)^2", 
	"WLS [gnuplot], weights n(h)",
	"WLS [gnuplot], weights n(h)/mod(h)^2", 
	"REML sills",
	"OLS (unwweighted)",
	"WLS, weights n(h)/(h * h)"
};

#define N_WIND 8
static const char *winds[N_WIND] = { "N","NE","E","SE","S","SW","W","NW" },
	*warning_msg = NULL;
static char *edf = "Enter data first";
char *header = NULL, *msg = NULL, *gnuplot_name = NULL;
static int var_ids[2] = { -1, -1 }, redraw_scr = 0, window_nr = -1,
	warning_set = 0, curses_open = 0;
static VARIOGRAM *v = NULL, *v_tmp = NULL;
SAMPLE_VGM_TYPE vgm_type = NOTSPECIFIED;

#if HAVE_LIBNCURSES || __PDCURSES__
# define GET_STRING(a,b) wgetnstr(stdscr, a, b)
#else
# define GET_STRING(a,b) wgetstr(stdscr, a)
#endif

#ifdef __hpux /* something's wrong on hp 10.20 */
# undef attrset
# define attrset(x)
# undef attron
# define attron(x)
# undef attroff
# define attroff(x)
#endif

#ifndef KEY_MOUSE
# define KEY_MOUSE -9999
#endif

#include "ui_help.h"

void enter_data(int key);
void choose_variable(int key);
void choose_what(int key);
void cutoff_width(int key);
void direction(int key);
void vgm_model(int key);
void show_plot(int key);
int save_as(int key);
void help_menu(int option);
void display_menu(int option);
void display_msg(const char *msg);
int calc_vgm(VARIOGRAM *v, char *name);
int fit_vgm(VARIOGRAM *v);
void calc_progress(unsigned int step, unsigned int max_step);
void fit_progress(unsigned int step, unsigned int max_step);
void more_file(const char *name);
void close_curses(void);
void open_curses(void);
void get_neighbourhood(int key);
void get_locations(int key);
void set_output(int key);
void display_map(int key);
PLOT_TYPE key2plot(int key);
void curses_printlog(const char *s);
void set_mouse_on(void);
int get_mouse_y(void);
int get_mouse_button(void);
void curses_warning(const char *msg);
void curses_error(const char *msg, int errno);

VGM_MODEL_TYPE toggle_model(int i, VGM_MODEL_TYPE m);

typedef enum PR_TYPE { PR_INT, PR_DOUBLE, PR_STRING, PR_YN } PR_TYPE ;
int prompt_for(const char *pr, PR_TYPE what, void *def_ret);

static enum { VGM = 0, PRED, XVALID, CONDSIM } what = VGM;

const char *pr_type_str[] = {
	"",
	"prediction",
	"cross validation",
	"conditional simulation" };

#define is_left_key(key) (key == KEY_LEFT || key == 'h')
#define ct_mvaddstr(y, str) mvaddstr(y,(COLS-strlen(str))/2,str)
#define STRLEN 256

#ifdef DJGPP /* save & restore start-up screen (from: gnuplot/term/djsvga.trm) */
#include <pc.h>
static char *dj_textsave=NULL;
static int dj_cursorx, dj_cursory;
#endif

#define MAX_OPTION 8
static char *text_a[MAX_OPTION] = {
	"enter/modify data", 
	"choose variable  ", 
	"calculate what   ", 
	"cutoff, width    ", 
	"direction        ", 
	"variogram model  ", 
	"fit method       ", 
	"show plot <Tab>  " };

static char *text_b[MAX_OPTION] = {
	"enter/modify data", 
	"choose variable  ", 
	"calculate what   ", 
	" neighbourhood   ", 
	" where           ", 
	"variogram model  ", 
	" output maps     ", 
	" show map <Tab>  " };

/*
 * Command structure:
 */
typedef struct {
	char *text, **help, *entry;
	void (*function)(int);
} Command;

Command command[MAX_OPTION] = {
 { NULL, ui0_help0, NULL, enter_data },
 { NULL, ui0_help1, NULL, choose_variable },
 { NULL, ui0_help2, NULL, choose_what },
 { NULL, ui0_help3, NULL, cutoff_width },
 { NULL, ui0_help4, NULL, direction },
 { NULL, ui0_help5, NULL, vgm_model },
 { NULL, ui0_help6, NULL, fit_method },
 { NULL, ui0_help7, NULL, show_plot }
};

int wind = 0, cmd_written = 0, use_pipe = 0;
static FILE *gnu_stream = NULL;

int start_ui(void) {
	int key = 0, id, option = 0, quit = 0, len = 0, i, j, my_fit, fix = 0, y;
	DO_AT_ZERO zero = ZERO_DEFAULT;
	char *s = NULL, shell_cmd[STRLEN];
	void fill_entries(void), setup_command(void);

	if (!(isatty(fileno(stdin)) && isatty(fileno(stdout))))
		ErrMsg(ER_IO, "user interface can only be run interactively");
/*
 * initialise variables 
 */
	if (gl_gnuplot35)
		gnuplot_name = gl_gnuplot35;
	else
		gnuplot_name = gl_gnuplot;

	if (gl_zero_est >= 0 && gl_zero_est <= ZERO_SPECIAL)
		zero = zero_int2enum(gl_zero_est);
	if (get_n_vars() > 0) 
		option = MAX_OPTION - 1; /* -> Show plot */
	command[1].entry = (char *) emalloc(STRLEN * sizeof(char));
	command[3].entry = (char *) emalloc(STRLEN * sizeof(char));
	command[4].entry = (char *) emalloc(STRLEN * sizeof(char));
/* use pipe? */
#ifdef HAVE_POPEN
	use_pipe = (getenv("DISPLAY") != NULL); /* so far: X11 only */
#endif
	v_tmp = init_variogram(v_tmp);
	if (gl_fit < NO_FIT || gl_fit >= LENGTHOF(fit_type_str)) /* out of range */
		gl_fit = NO_FIT;
	my_fit = gl_fit;
	len = 10 + strlen((char *) VERSION) +
		(command_file_name ? strlen(command_file_name) : 0);
	header = (char *) emalloc(len);
	if (command_file_name)
		sprintf(header, "gstat %s, %s", VERSION, command_file_name);
	else
		sprintf(header, "gstat %s", VERSION);
	if (get_n_masks() > 0) 
		ErrMsg(ER_IMPOSVAL, "delete mask definition");
	if (get_n_vars() > 0) {
		printf("press return to continue...");
		get_line(&s, &len, stdin);
		efree(s);
	}
	initscr();
	open_curses();
	atexit(close_curses);
	while (! quit) { /* enter event loop */
		if (get_n_vars() && var_ids[0] < 0) {
			var_ids[0] = var_ids[1] = 0;
			if (get_method() == COV)
				vgm_type = COVARIOGRAM;
			else
				vgm_type = SEMIVARIOGRAM;
		}
		if (var_ids[0] > -1) {
			if (my_fit != NO_FIT) { /* First time only  */
				for (i = 0; i < get_n_vars(); i++) {
					for (j = 0; j <= i; j++) {
						v = get_vgm(LTI(i,j));
						if (v->n_models > 0)
							v->ev->fit = fit_int2enum(my_fit);
						v->ev->refit = 1;
					}
				}
				my_fit = NO_FIT;
			}
			id = LTI(var_ids[0],var_ids[1]);
			v = get_vgm(id);
			v->ev->plot_numbers = gl_numbers;
			v->id = id;
			v->id1 = var_ids[0];
			v->id2 = var_ids[1];
			if (v->ev->evt != vgm_type) {
				v->ev->evt = vgm_type;
				v->ev->recalc = 1;
			}
			if (gl_gls_residuals != 0)
				v->ev->recalc = 1;
			if (!is_direct(v) && v->ev->fit == MIVQUE_FIT)
				v->ev->fit = NO_FIT;
		} 
		setup_command();
		fill_entries();
		display_menu(option);
		switch (key = getch()) { /* wait for a key press... */
			case KEY_MOUSE:
				/*
				sprintf(shell_cmd, "line %d", y);
				msg = shell_cmd;
				*/
				y = get_mouse_y();
				if (y >= 4 && y < (4 + MAX_OPTION)) {
					option = y - 4;
					display_menu(option);
					switch (get_mouse_button()) {
						case 2:
							help_menu(option);
							break;
						case 3:
							key = KEY_LEFT; /* it is right, but I want left */
							(command[option].function)(key);
							break;
						default:
							key = KEY_RIGHT;
							(command[option].function)(key);
							break;
					}
				}
				break;
			case 10:                                 /* do it */
			case 13:
			case KEY_ENTER:
			case KEY_LEFT:
			case KEY_RIGHT:
			case 'l':
			case 'h':
				(command[option].function)(key);
				break;
			case 'k':                                /* go down */
			case KEY_UP:
				if (option > 0)
					option--;
				break;
			case 'j':                                /* go up */
			case KEY_DOWN:
				if (option < MAX_OPTION - 1)
					option++;
				break;
			case KEY_HOME:                           /* go top */
			case KEY_PPAGE:
				option = 0;
				break;
			case KEY_NPAGE:
				option = MAX_OPTION - 1;
				break;
			case '\t':                               /* show plot */
				show_plot(0);
				break;
			case 'c': 
			case 'e': /* Save to file... */
			case 'E':
			case 'G': 
			case 'g': 
			case 'L': 
			case 'p': 
			case 'P': 
			case 'J': 
			case '&':
				save_as(key);
				break;
			case 'f':                                /* toggle range fixing */
				fix = !fix;
				fix_toggle(v, fix);
				break;
			case 'n':                                /* toggle number displ */
			case 'N':
				if (v == NULL) {
					msg = edf;
					break;
				}
				if (! v->ev->plot_numbers) {
					msg = "number display on";
					gl_numbers = 1;
				} else {
					msg = "number display off";
					gl_numbers = 0;
				}
				break;
			case 'S': case 's':                   /* Show estimates */
				shell_cmd[0] = '\0';
				if (v->fname || o_filename)
					strcpy(shell_cmd, v->fname ? v->fname : o_filename);
				prompt_for("\rFile to show", PR_STRING, shell_cmd);
				clear();
				more_file(shell_cmd);
				redraw_scr = 1;
				break;
			case 'w': case 'W':                   /* toggle window */
				if (window_nr <= 0)
					window_nr = 1;
				else
					window_nr++;
				prompt_for("Window number", PR_INT, &window_nr);
				break;
			case ' ':
			case '+':                             /* Increment variable */
				choose_variable('l');
				break;
			case '-':                             /* Decrement variable */
				choose_variable('h');
				break;
			case 'z':                             /* toggle zero est */
			case 'Z':
				if (v == NULL || get_n_vars() <= 0) {
					msg = edf;
					break;
				} else
					v->ev->refit = v->ev->recalc = 1;
				switch (zero) {
					case ZERO_DEFAULT: /* BREAKTHROUGH: */
					case ZERO_SPECIAL:
						zero = ZERO_INCLUDE;
						msg = "zero distances included in first interval";
						break;
					case ZERO_INCLUDE:
						zero = ZERO_AVOID;
						msg = "zero distances ignored";
						break;
					case ZERO_AVOID:
						zero = ZERO_SPECIAL;
						msg = "separate estimates at zero distance";
						break;
				}
				v->ev->zero = zero_int2enum(gl_zero_est = zero);
				break;
			case 12:
			case 'R':
			case 'r':                                  /* redraw screen */
				redraw_scr = 1;
				break;
			case 'H':                                  /* Help screen */
				help_menu(-1);
				break;
			case '?':
				help_menu(option);                      /* Help option */
				break;
			case '!':                                  /* Shell command */
				addstr("!");
				refresh();
				echo();
				keypad(stdscr, FALSE);
				GET_STRING(shell_cmd, STRLEN);
				keypad(stdscr, TRUE);
#ifdef __PDCURSES__ 
				close_curses();
#endif
				esystem(shell_cmd);
				redraw_scr = 1;
				noecho();
				break;
			case 'q':
			case 'Q':
				addstr("Quit");
				quit = 1;
				if (! cmd_written) { 
					addstr("\n");
					if (prompt_for("Current settings will be lost. Save to file?",
							PR_YN, &cmd_written) && cmd_written) { /* a Yes */
						if (save_as('c'))
							cmd_written = quit = 0;
					} else
						prompt_for("Really quit?", PR_YN, &quit);  /* a Yes */
				}
				redraw_scr = 1;
				break;
			case 'x':                                  /* Quit menu */
			case 'X':
				addstr("Exit");
				quit = 1;
				break;
			case ';':                                  /* switch mode */
				/*
				choose_what(key);
				redraw_scr = 1;
				break;
				*/
			default:
				msg = "Unknown command (use H for help)";
				break;
		} /* switch(key) */
	} /* exit event loop while (! quit) */
#ifdef HAVE_POPEN
	if (gnu_stream != NULL)
		epclose(gnu_stream);
#endif
	move(LINES - 1, 0);
	clrtoeol();
	refresh();
	close_curses();
#ifdef DJGPP /* on normal program termination, rebuild original screen */
	ScreenUpdate(dj_textsave);
	ScreenSetCursor(dj_cursory, dj_cursorx);
#endif
	return (key == 'x' || key == 'X'); /* we quit by pressing exit */
}

void enter_data(int key) {
	int yes_no, i, id = -1, n_at_start;
	char s[STRLEN];
	DATA **d = NULL, *data = NULL;
	VARIOGRAM *vp = NULL;

	s[0] = '\0';
	n_at_start = get_n_vars();
	redraw_scr = 1;
	do {
		clear();
		echo();
		if (id == -1) { /* at start */
			if (key != -1) {
				if (n_at_start) {
					sprintf(s, "%s", name_identifier(0));
					printw("Identifiers so far:\n");
					for (i = 0; i < n_at_start; i++)
						printw("\t%s%s", name_identifier(i), 
							i > 0 && i % 5 == 0 ? "\n" : "");
				}
				if (!prompt_for("\nChoose identifier or enter new one", PR_STRING, s)
						&& n_at_start == 0)
					return;
				id = which_identifier(s);
				d = get_gstat_data();
				data = d[id];
			} else {
				data = get_dataval();
				id = ID_OF_VALDATA;
			}
			data->id = id;
			data->dummy = 0;
		} else {
			yes_no = 0;
			if (warning_set)
				warning_set = 0;
			if (msg)
				addstr(msg);
			msg = NULL;
			prompt_for("\nQuit program now?", PR_YN, &yes_no);
			if (yes_no != 0)
				exit(0);
		}
		printw("definition of `%s' (variable %d):\n", name_identifier(id), id+1);
		if (id >= n_at_start)
			data->average = 0;
		sprintf(s, "%s", data->fname ? data->fname : "");
		if (prompt_for("Enter data file name", PR_STRING, s))
			data->fname = string_dup(s);
		if (!grass() && !file_exists(data->fname))
			msg = "File not found!";
		else {
			if (! grass()) {
				prompt_for("Enter x coordinate column number", PR_INT, &(data->colnx));
				prompt_for("Enter y coordinate column number", PR_INT, &(data->colny));
				if (data->colnz)
					prompt_for("Enter z coordinate column number",
							PR_INT, &(data->colnz));
			}
			prompt_for("Enter variable column number", PR_INT, &(data->colnvalue));
			prompt_for("Log transform data", PR_YN, &(data->log));
			prompt_for("Average duplicate observations", PR_YN, &(data->average));
			data->n_averaged = 0;
		}
	} while (msg || read_gstat_data(data) == NULL || data_not_valid(data));
	data->is_residual = 0;
	msg = "Data read succesfully";

	set_gstat_log_handler(curses_printlog);
	report_data(data);
	set_gstat_log_handler(default_printlog);
	addstr("Press any key to continue...");
	refresh();
	getch();

	for (i = 0; key != -1 && i < get_n_vars(); i++) {
		vp = get_vgm(LTI(i, id));
		vp->ev->recalc = vp->ev->refit = 1;
	}
	return;
}

void cutoff_width(int key) {
	double cutoff, iwidth;

	if (v == NULL) {
		msg = edf;
		return;
	}
	if (gl_bounds != NULL) {
		msg = "fixed boundaries: modify command (or boundary) file";
		return;
	}
	cutoff = v->ev->cutoff;
	iwidth = v->ev->iwidth;
	redraw_scr = 1;
	prompt_for("\rNew cutoff: ", PR_DOUBLE, &cutoff);
	prompt_for("New width:  ", PR_DOUBLE, &iwidth);
	if (cutoff > 0.0 && iwidth >= 0.0 && cutoff >= iwidth) {
		gl_cutoff = v->ev->cutoff = cutoff;
		gl_iwidth = v->ev->iwidth = iwidth;
		v->ev->recalc = 1;
		v->ev->refit = 1;
	} else
		msg = "Wrong values!";
}

void vgm_model(int key) {
	char vgm_s[STRLEN] = "";
	VGM_MODEL_TYPE m;
	VGM_MODEL *vp;
	int i;

	if (v == NULL || v->id < 0) {
		msg = edf;
		return;
	}
	redraw_scr = 1;
	prompt_for("\rNew model", PR_STRING, vgm_s);
	if (vgm_s[0] != '\0') {
		if (read_variogram(v_tmp, vgm_s)) {
			msg = "Error on reading variogram model";
			vgm_s[0] = '\0';
		} else {
			read_variogram(v, vgm_s);
			update_variogram(v);
			v->ev->refit = 1;
		}
	} else if (v->n_models > 0) {
		for (i = 0; i < v->n_models; i++) {
			vp = &(v->part[i]);
			move(LINES - 5, 0); /* place cursor */
			sprintf(vgm_s, "\rModel %d, sill: ", i+1);
			prompt_for(vgm_s, PR_DOUBLE, &(vp->sill));
			m = toggle_model(i, vp->model);
			sprintf(vgm_s, "\rModel %d, range:", i+1);
			prompt_for(vgm_s, PR_DOUBLE, &(vp->range[0]));
			if (vp->tm_range != NULL) {
				move(LINES - 5, 0); /* place cursor */
				sprintf(vgm_s, "\rModel %d, angle 1:", i+1);
				prompt_for(vgm_s, PR_DOUBLE, &(vp->tm_range->angle[0]));
				move(LINES - 5, 0); /* place cursor */
				sprintf(vgm_s, "\rModel %d, angle 2:", i+1);
				prompt_for(vgm_s, PR_DOUBLE, &(vp->tm_range->angle[1]));
				move(LINES - 5, 0); /* place cursor */
				sprintf(vgm_s, "\rModel %d, angle 3:", i+1);
				prompt_for(vgm_s, PR_DOUBLE, &(vp->tm_range->angle[2]));
				move(LINES - 5, 0); /* place cursor */
				sprintf(vgm_s, "\rModel %d, ratio 1:", i+1);
				prompt_for(vgm_s, PR_DOUBLE, &(vp->tm_range->ratio[0]));
				move(LINES - 5, 0); /* place cursor */
				sprintf(vgm_s, "\rModel %d, ratio 2:", i+1);
				prompt_for(vgm_s, PR_DOUBLE, &(vp->tm_range->ratio[1]));
			}
			if (m != NUGGET && m != LINEAR && m != INTERCEPT && m != MERROR 
					&& vp->range[0] == 0.0)
				msg = "range should be positive";
			else if ((m == NUGGET || m == INTERCEPT) && vp->range[0] != 0.0)
				msg = "range should be zero";
			else
				vp->model = m;
		}
		update_variogram(v);
		read_variogram(v, v->descr); /* sets function pointers */
		v->ev->refit = 1;
	} 
}

VGM_MODEL_TYPE toggle_model(int i, VGM_MODEL_TYPE m) {
	int key = 0;

	/* NOT_SP ... INTERCEPT */
	do {
		switch (key) {
			case 0: break; /* don't change */
			case 'n': case 'N': m = NUGGET; break;
			case 'l': case 'L': m = LINEAR; break;
			case 'c': case 'C': m = CIRCULAR; break;
			case 's': case 'S': m = SPHERICAL; break;
			case 'p': case 'P': m = PENTASPHERICAL; break;
			case 'e': case 'E': m = EXPONENTIAL; break;
			case 'g': case 'G': m = GAUSSIAN; break;
			case 'b': case 'B': m = BESSEL; break;
			case KEY_LEFT: case KEY_UP: case 'k': 
				/* m -= 1; */
				m = model_shift(m, 0);
				break;
			default: 
				/* m = (m == INTERCEPT ? NUGGET : m + 1); */
				m = model_shift(m, 1);
				break;
		}
		/* 
		if (m == NOT_SP)
			m = INTERCEPT;
		*/
		move(LINES - 5, 0); /* place cursor */
		clrtoeol(); /* and clear line */
		noecho();
		printw("Model %d, type: %s", i+1, v_models[m].name_long);
		refresh();
	} while ((key = getch()) != KEY_ENTER && key != 10 && key != 13);
	echo();
	return m;
}


void show_plot(int key) {

	char *cp = NULL;

	if (v == NULL) {
		msg = edf;
		return;
	}
	printw("Show");
	refresh();
	calc_vgm(v, v->fname ? v->fname : o_filename);
	if (fit_vgm(v)) {
		msg = "Error on fit";
		return;
	}
#ifdef HAVE_POPEN
	if (use_pipe) {
		if (gnu_stream == NULL)
			gnu_stream = epopen(gnuplot_name, "w");
		fprint_gnuplot_variogram(gnu_stream, v, 0, key2plot(key), window_nr);
	} else
#endif
	{
		gnu_stream = efopen(gl_plotfile, "w");
		fprint_gnuplot_variogram(gnu_stream, v, 0, key2plot(key), window_nr);
		efclose(gnu_stream);
		gnu_stream = NULL; /* don't epclose() it later on... */
		cp = (char *) emalloc(2 + strlen(gnuplot_name) + strlen(gl_plotfile));
		sprintf(cp, "%s %s", gnuplot_name, gl_plotfile);
		close_curses(); /* give gnuplot complete terminal control */
		if (esystem(cp)) /* do it */
			msg = "Error (non-zero exit) from gnuplot";
		efree(cp);
		if ((cp = getenv("GNUTERM")) && strstr(cp, "dumb"))
			getch();
	}
	if (v->ev->fit == WLS_GNUFIT || v->ev->fit == WLS_GNUFIT_MOD)
		redraw_scr = 1; /* messes up screen! */
	return;
}

void display_menu(int option) {
	int i;
	static int last = -1;

	if (warning_set) {
		msg = (char *) warning_msg;
		warning_set = 0;
	} else
		warning_msg = NULL;
	if (last == -1 || curses_open == 0 || redraw_scr) { /* full menu */
		open_curses();
		clear();
		attrset(A_NORMAL);
		ct_mvaddstr(1, header);
		if (LINES > 20) {
			ct_mvaddstr(LINES - 8, SHORTHELP0);
			ct_mvaddstr(LINES - 7, SHORTHELP1);
		}
		mvaddstr(LINES - 5, 0, "Command: ");
		for (i = 0; i < MAX_OPTION; i++) {
			mvaddstr(4+i, 1, command[i].text);
			if (command[i].entry != NULL && *(command[i].entry) != '\0') {
				mvaddstr(4+i, 19, ":");
				mvaddstr(4+i, 21, command[i].entry);
			}
		}
		mvaddstr(4+option, 0, ">");
		attrset(A_REVERSE);
		mvaddstr(4+option, 1, command[option].text);
		attrset(A_NORMAL);
		redraw_scr = 0;
	}
	if (last != option) { /* only a cursor move */
		if (last > -1) {
			mvaddstr(4+last, 0, " ");
			mvaddstr(4+last, 1, command[last].text);
		}
		attrset(A_REVERSE);
		mvaddstr(4+option, 0, ">");
		mvaddstr(4+option, 1, command[option].text);
		attrset(A_NORMAL);
		refresh();
		last = option;
	} else { /* no cursor move, but an action on option/last */
		for (i = 0; i < MAX_OPTION; i++) {
			move(4+i, 19);
			clrtoeol();
			if (command[i].entry != NULL && *(command[i].entry) != '\0') {
				mvaddstr(4+i, 19, ":");
				mvaddstr(4+i, 21, command[i].entry);
			}
		}
	}
	move(LINES - 1, 0);
	clrtoeol();
	display_msg(msg);
	msg = NULL;
	last = option;
	move(LINES - 5, 9); /* place cursor */
	clrtoeol(); /* and clear line */
	noecho();
	/* nl(); ?? */
	refresh();
	return;
}

void display_msg(const char *msg) {

	if (msg != NULL) {
		move(LINES - 1, 0); /* place cursor */
		clrtoeol();
		ct_mvaddstr(LINES - 1, (char *) msg);
		refresh();
	} 
}

int calc_vgm(VARIOGRAM *v, char *name) {
	set_gstat_progress_handler(calc_progress);
	return calc_variogram(v, name);
}

int fit_vgm(VARIOGRAM *v) {
	set_gstat_progress_handler(fit_progress);
	return fit_variogram(v);
}

int prompt_for(const char *pr, PR_TYPE what, void *def_ret) {

	char s[STRLEN];
	int c = 0, i, y, x;
	char *cp = NULL;
	double d;

	if (pr[0] == '\r') {
		getyx(stdscr, y, x);
		move(y, 0);
		clrtoeol();
		pr++;
	}
	printw("%s", pr);
	switch(what) {
		case PR_INT:
			printw(" [%d]", *(int *)def_ret);
			break;
		case PR_DOUBLE: 
			if (is_mv_double(def_ret)) 
				printw(" [NA]"); 
			else
				printw(" [%g]", *(double *)def_ret); 
			break;
		case PR_STRING:
			printw(" [%s]", (char *)def_ret);
			break;
		case PR_YN:
			printw(" [%c]", (*(int *)def_ret) ? 'y' : 'n');
			break;
	}
	printw(" : ");
	refresh();
	if (what == PR_YN) {
		noecho();
		c = getch();
		if (c != 10 && c != 13 && c != KEY_ENTER && c != KEY_MOUSE) 
				/* meaning some useful letter was typed: */
			*(int *)def_ret = (c == 'y' || c == 'Y');
		if (*(int *)def_ret)
			addstr("Yes\n");
		else
			addstr("No\n");
		refresh();
		echo();
		return 1;
	} else {
		echo();
		keypad(stdscr, FALSE);
		GET_STRING(s, STRLEN);
		keypad(stdscr, TRUE);
	}
	if (s[0] == '\n' || s[0] == '\0')
		return 0;
	if ((cp = strchr(s, '\n')) != NULL)
		*cp = '\0';
	switch (what) {
		case PR_INT:
			if (read_int(s, &i))
				msg = "error reading integer";
			else
				(* (int *) def_ret) = i;
			break; 
		case PR_DOUBLE:
			if (read_double(s, &d))
				msg = "error reading real value";
			else
				(* (double *) def_ret) = d;
			break;
		case PR_STRING: strcpy((char *) def_ret, s); break;
		case PR_YN: break;
	}
	return 1;
}

void choose_variable(int key) {
	int n;
	n = get_n_vars();
	if (n == 0) {
		msg = edf;
		return;
	}
	if (n == 1) {
		var_ids[0] = var_ids[1] = 0;
		msg = "Only one variable available";
		return; /* nothing to shift */
	}

	if (! is_left_key(key)) { /* not going back */
		if (var_ids[0] == var_ids[1]) {
			if (var_ids[0] < n - 1)
				var_ids[0] = var_ids[1] = var_ids[0] + 1;
			else
				var_ids[0] = var_ids[1] = 0; 
		} else {
			if (var_ids[1] == n - 1) {
				if (var_ids[0] == n - 2) {
					var_ids[0] = 0;
					var_ids[1] = 1;
				} else { 
					var_ids[0]++;
					var_ids[1] = var_ids[0] + 1;
				}
			} else {
				var_ids[1]++;
			}
		}
	} else {
		if (var_ids[0] == var_ids[1]) {
			if (var_ids[0] > 0)
				var_ids[0] = var_ids[1] = var_ids[0] - 1;
			else
				var_ids[0] = var_ids[1] = n - 1; /* last */
		} else {
			if (var_ids[1] == 1 && var_ids[0] == 0) { /* first cross */
				var_ids[0] = n - 2;
				var_ids[1] = n - 1; /* jump to last variogram */
			} else {
				if (var_ids[1] - var_ids[0] > 1)
					var_ids[1]--;
				else {
					var_ids[0]--;
					var_ids[1] = n - 1;
				}
			}
		}
	}
}

void choose_what(int key) {

	if (v == NULL) {
		vgm_type = NOTSPECIFIED;
		msg = edf;
		return;
	}
	if (key == ';') {
		what = (what == VGM) ? PRED : (get_n_masks() == 0 ? VGM : what);
		return;
	}
	switch (what) {
		case PRED: what = XVALID; gl_xvalid = 1; return;
		case XVALID: what = CONDSIM; gl_xvalid = 0; return;
		case CONDSIM: what = PRED; gl_xvalid = 0; return;
		default: break;
	}
	if (get_n_vars() == 1) {
		vgm_type = (vgm_type == SEMIVARIOGRAM) ? COVARIOGRAM : SEMIVARIOGRAM;
		v->ev->recalc = v->ev->refit = 1;
		return;
	}
	if (is_left_key(key)) { /* decrease */
		if (vgm_type == SEMIVARIOGRAM)
			vgm_type = CROSSCOVARIOGRAM;
		else
			vgm_type--;
	} else { /* increase */
		if (vgm_type == CROSSCOVARIOGRAM)
			vgm_type = SEMIVARIOGRAM;
		else
			vgm_type++;
	}
	if (var_ids[0] == var_ids[1]) {
		if (!(vgm_type == SEMIVARIOGRAM || vgm_type == COVARIOGRAM)) {
			if (var_ids[1] < get_n_vars() - 1)
				var_ids[1]++;
			else
				var_ids[0]--;
		}
	} else {
		if (vgm_type == SEMIVARIOGRAM || vgm_type == COVARIOGRAM)
			var_ids[1] = var_ids[0];
	}
	v->ev->recalc = v->ev->refit = 1;
	return;
}

void fit_method(int key) {
#define MAX_FIT (LENGTHOF(fit_type_str) - 1)
	if (v == NULL) {
		msg = edf;
		return;
	}
	if (v->n_models <= 0) {
		msg = "Enter variogram model first";
		return;
	}
	v->ev->fit = fit_shift(v->ev->fit, !is_left_key(key));
	if (! is_variogram(v)) {
		if (v->ev->fit == WLS_FIT_MOD || v->ev->fit == WLS_GNUFIT_MOD ||
				v->ev->fit == WLS_NHH)
			msg = "this is NOT a recommended fitting method for covariograms";
	}
	gl_fit = v->ev->fit;
	v->ev->refit = 1;
}

void setup_command(void) {
	int i;
	char **cpp;

	cpp = (what == VGM) ? text_a : text_b;
	for (i = 0; i < MAX_OPTION; i++)
		command[i].text = cpp[i];
 	command[0].function = enter_data;
 	command[1].function = choose_variable;
 	command[2].function = choose_what;
 	command[5].function = vgm_model;
 	if (what == VGM) {
 		command[3].function = cutoff_width;
 		command[4].function = direction;
		command[6].function = fit_method;
		command[7].function = show_plot;
 	} else {
 		command[3].function = get_neighbourhood;
 		command[4].function = get_locations;
		command[6].function = set_output;
		command[7].function = display_map;
	}
}

void fill_entries(void) {
	DATA **dpp = NULL, *d = NULL;
	int i, j;

	/*
	 * entry 1: variables 
	 */
	if (var_ids[0] >= 0) {
		dpp = get_gstat_data();
		d = dpp[var_ids[0]];
		if (var_ids[0] == var_ids[1])
			sprintf(command[1].entry, "%s", name_identifier(var_ids[0]));
		else
			sprintf(command[1].entry, "%s and %s", 
				name_identifier(var_ids[0]), 
				name_identifier(var_ids[1]));
	} else
		command[1].entry[0] = '\0';
	/*
 	 * entry 2: calculate what?
 	 */
	if (what == VGM)
		command[2].entry = (char *) vgm_type_str[vgm_type];
	else
		command[2].entry = (char *) pr_type_str[what];
	/*
 	 * entry 3: cutoff, width <--> neighbourhood
 	 */
 	if (what == VGM) {
		if (v != NULL) {
			fill_cutoff_width(d, v);
			if (gl_bounds != NULL)
				sprintf(command[3].entry, "[%g fixed boundaries]", 
					v->ev->cutoff);
			else
				sprintf(command[3].entry, "%g, %g", v->ev->cutoff, v->ev->iwidth);
		} else
			command[3].entry[0] = '\0';
	} else {
		if (d && !IS_GLOBAL(d)) {
			sprintf(command[3].entry, "min %d, max %d, radius %g",
				d->sel_min, d->sel_max, d->sel_rad);
		} else
			sprintf(command[3].entry, "(global)");
	}
	/*
 	 * entry 4: variogram model -- where
 	 */
	if (what == VGM) {
		if (v != NULL) {
			if (! is_directional(v))
				sprintf(command[4].entry, "total");
			else { 
			if (d != NULL && (d->mode & Z_BIT_SET))
				sprintf(command[4].entry,
						"horiz. %gd +/- %g (%s), vert. %gd +/- %g",
						gl_alpha, gl_tol_hor, winds[wind], gl_beta, gl_tol_ver);
			else
				sprintf(command[4].entry, "%gd +/- %g (%s)",
					gl_alpha, gl_tol_hor, winds[wind]);
			}
		} else
			command[4].entry[0] = '\0';
	} else {
		d = get_dataval();
		if (get_n_masks() > 0) {
			sprintf(command[4].entry, "%s (mask map)", get_mask_name(0));
		} else if (d->id > -1) {
			sprintf(command[4].entry, "%s (point data)", name_identifier(d->id));
		} else
			sprintf(command[4].entry, "<undefined>");
	}
	/*
 	 * entry 5: variogram model
 	 */
	if (v != NULL)
		command[5].entry = v->descr;
	else
		command[5].entry = NULL;
	/*
 	 * entry 6: fit -- output maps
 	 */
 	if (what == VGM) {
		if (v && v->n_models > 0)
			command[6].entry = fit_type_str[v->ev->fit];
		else
			command[6].entry = NULL;
 	} else {
 		for (i = j = 0; i < get_n_outfile(); i++)
 			if (get_outfile_namei(i))
 				j++;
 		if (j == get_n_outfile())
 			command[6].entry = "all set";
 		else if (j > 0)
 			command[6].entry = "some set";
 		else
 			command[6].entry = "<not set>";
 	}
	return;
}

void direction(int key) {
	DATA **data;
	int mode;

	if (v == NULL) {
		msg = edf;
		return;
	}
	data = get_gstat_data();
	mode = data[v->id1]->mode | data[v->id2]->mode;
	if (! mode & Y_BIT_SET) {
		msg = "no directions in 1-d";
		return;
	}
	redraw_scr = 1;
	prompt_for("\rHorizontal direction (deg): ", PR_DOUBLE, &gl_alpha);
	prompt_for("Horizontal tolerance (deg): ", PR_DOUBLE, &gl_tol_hor);
	if (mode & Z_BIT_SET) { /* 3D-data */
	prompt_for("Vertical direction (deg):   ", PR_DOUBLE, &gl_beta);
	prompt_for("Vertical tolerance (deg):   ", PR_DOUBLE, &gl_tol_ver);
	}
	if (gl_alpha < 0.0 || gl_alpha > 360.0 ||
			gl_beta < 0.0 || gl_beta > 360.0 ||
			gl_tol_hor < 0.0 || gl_tol_hor > 180.0 ||
			gl_tol_ver < 0.0 || gl_tol_ver > 180.0) {
		msg = "Wrong values!";
		gl_alpha = DEF_alpha;
		gl_beta = DEF_beta;
		gl_tol_hor = DEF_tol_hor;
		gl_tol_hor = DEF_tol_hor;
	}
	v->ev->recalc = 1;
	v->ev->refit = 1;
	/* sem.c does set_direction_values(); */
	v->ev->is_directional = -1; /* not yet determined */
	wind =  (int) ceil((gl_alpha/360.0) * 8 + 0.5) - 1;
	if (wind < 0 || wind > 7) {
		wind = 0;
		msg = "Wrong wind!";
	}
	return;
}

int save_as(int key) {
/*
 * return 0 on succesful save, 1 on unsuccesful
 */
	char fname[STRLEN], prompt[STRLEN], cmd[STRLEN];
	FILE *tmp = NULL;
	int yes_no = 1;

	if (v == NULL) {
		msg = edf;
		return 1;
	}
	*fname = '\0';
	if (var_ids[0] == var_ids[1] || key == 'c')
		sprintf(fname, "%s",
			name_identifier(var_ids[0]));
	else
		sprintf(fname, "%s%s", 
			name_identifier(var_ids[0]),
			name_identifier(var_ids[1]));
	switch(key) {
		case 'c': 
			strcat(fname, ".cmd");
			addstr("Save configuration as gstat command file");
			break;
		case 'e': 
			strcat(fname, ".est");
			if (v->fname)
				sprintf(fname, "%s", v->fname);
			addstr("Save from now on estimates to file");
			break;
		case 'g': 
			strcat(fname, ".gnu");
			addstr("Save from now on gnuplot commands to file");
			break;
		case 'G': 
			strcat(fname, ".gif");
			addstr("Save plot as gif file");
			break;
		case 'L':
			strcat(fname, ".tex");
			addstr("Save plot as LaTeX/postscript file");
			break;
		case 'E':
			strcat(fname, ".tex");
			addstr("Save plot as LaTeX/eepic file");
			break;
		case 'P':
			strcat(fname, ".eps");
			addstr("Save plot as encapsulated postscript file");
			break;
		case 'p':
			strcat(fname, ".png");
			addstr("Save plot as PNG file");
			break;
		case '&':
			strcat(fname, ".eps");
			addstr("Save plot as encapsulated postscript file (through jgraph)");
			break;
		case 'J':
			strcat(fname, ".jgr");
			addstr("Save as jgraph file");
			break;
	}
	redraw_scr = 1;
	addstr("\n");
	if (file_exists(fname))
		sprintf(prompt, "Overwrite file `%s' ?", fname);
	else
		sprintf(prompt, "Write to file `%s' ?", fname);
	prompt_for(prompt, PR_YN, &yes_no);
	if (yes_no == 0) {
		*fname = '\0';
		sprintf(prompt, "Enter file name");
		if (prompt_for(prompt, PR_STRING, fname) == 0 || *fname == '\0')
			return 1;
	}
	switch(key) {
		case 'c':
			fprint_cmd(tmp = efopen(fname, "w"));
			if (fclose(tmp) == 0) {
				cmd_written = 1;
				msg = "Commands written to file";
			} else
				msg = "Error on closing file";
			break;
		case 'e':
			o_filename = string_dup(fname);
			if (calc_vgm(v, fname))
				msg = "Error on closing file";
			else
				msg = "Estimates written to file";
			break;
		case 'g':
			gl_plotfile = string_dup(fname);
			tmp = efopen(fname, "w");
			if (calc_vgm(v, o_filename)) {
				msg = "Error on closing estimates file";
				return 1;
			}
			if (fit_vgm(v)) {
				msg = "Error on fit";
				return 1;
			}
			fprint_gnuplot_variogram(tmp, v, "gnuplot.out", 
							key2plot(key), window_nr);
			if (efclose(tmp) == 0)
				msg = "Gnuplot commands written to file";
			break;
		case 'p': 
		case 'P': 
		case 'G': 
		case 'L':
		case 'E':
			if (calc_vgm(v, o_filename)) {
				msg = "Error on closing file";
				return 1;
			}
			if (fit_vgm(v)) {
				msg = "Error on fit";
				return 1;
			}
			tmp = efopen(gl_plotfile, "w");
			fprint_gnuplot_variogram(tmp, v, fname, key2plot(key), window_nr);
			if (efclose(tmp)) {
				msg = "Error on closing file";
				return 1;
			}
			sprintf(cmd, "%s %s", gnuplot_name, gl_plotfile);
#ifdef __PDCURSES__ 
			close_curses();
#endif
			if (esystem(cmd))
				msg = "Error (non-zero exit) from gnuplot";
			else
				msg = "Graph written to file (through gnuplot)";
			break;
		case '&':
			calc_vgm(v, NULL);
			if (fit_vgm(v)) {
				msg = "Error on fit";
				return 1;
			}
			tmp = efopen("gstat.jgr", "w");
			fprint_jgraph_variogram(tmp, v);
			if (efclose(tmp)) {
				msg = "Error on closing file";
				return 1;
			}
			sprintf(cmd, "jgraph gstat.jgr > %s", fname);
#ifdef __PDCURSES__ 
			close_curses();
#endif
			if (esystem(cmd))
				msg = "Error (non-zero exit) from jgraph";
			else
				msg = "Graph written to eps file (via jgraph)";
			break;
		case 'J':
			if (calc_vgm(v, o_filename)) {
				msg = "Error on closing file";
				return 1;
			}
			if (fit_vgm(v)) {
				msg = "Error on fit";
				return 1;
			}
			tmp = efopen(fname, "w");
			fprint_jgraph_variogram(tmp, v);
			if (efclose(tmp))
				msg = "Error on closing file";
			else
				msg = "Graph written to jgraph file";
			break;
		default:
			ErrMsg(ER_IMPOSVAL, "in switch, save_as()");
	}
	return 0;
}

PLOT_TYPE key2plot(int key) {
	switch (key) {
		case 'L': return PSLATEX;
		case 'C': return CGM;
		case 'G': return GIF;
		case 'P': return EPS;
		case 'p': return PNG;
		case 'E': return EEPIC;
		case 'g': 
		default: return GNUPLOT;
	}
}

void help_menu(int option) {
	int i;
	char **hlp = NULL;

	redraw_scr = 1;
	clear();
	if (option > -1)
		printw("help on field: %s\n\n", command[option].text);
	switch (option) {
		case -1: hlp = helpmenu; break;
		case 0: case 1: case 2: case 3: case 4: case 5: case 6:
		case 7: hlp = command[option].help; break;
		default: ErrMsg(ER_IMPOSVAL, "help_menu()"); break;
	}
	for (i = 0; hlp[i]; i++)
		printw("%s\n", hlp[i]);
	printw("\n\npress a key to continue...");
	refresh();
	getch();
}

void more_file(const char *name) {
	char *cmd = NULL;
	static char *pager = NULL;

	if (name == NULL) {
		msg = "NULL file name";
		return;
	}
	if (! file_exists(name)) {
		msg = "file not found";
		return;
	}
	if (file_size(name) == 0) {
		msg = "empty file";
		return;
	}
	if ((pager = getenv("PAGER")) == NULL)
		pager = gl_pager;
	cmd = emalloc(strlen(pager) + strlen(name) + 2);
	sprintf(cmd, "%s %s", pager, name);
#ifdef __PDCURSES__ 
	close_curses();
#endif
	esystem(cmd);
	efree(cmd);
	return;
}

void curses_error(const char *msg, int errno) {

	close_curses();
	
	if (warning_set)
		fprintf(stderr, "%s\n", warning_msg);
	fprintf(stderr, "%s\n", msg);
	exit(errno);
}

void curses_warning(const char *msg) {
	warning_msg = msg;
	warning_set = 1;
}

void calc_progress(unsigned int step, unsigned int max_step) {
	char s[100];

	if (step < max_step)
		sprintf(s, "sample variogram: %d%% done", 
			(int) floor((100.0 * step)/ max_step));
	else
		sprintf(s, "Ready");
	display_msg(s);
}

void fit_progress(unsigned int step, unsigned int max_step) {
	char s[100];

	if (step < max_step)
		sprintf(s, "fit iteration: %d", step + 1);
	else
		sprintf(s, "Ready");
	display_msg(s);
}

void close_curses(void) {
	if (curses_open)
		endwin();
	curses_open = 0;
	return;
}

void open_curses(void) {
#ifdef DJGPP
	if (dj_textsave == NULL) {
		dj_textsave=(char *)emalloc(ScreenRows()*ScreenCols()*2*sizeof(char));
		ScreenRetrieve(dj_textsave);
		ScreenGetCursor(&dj_cursory,&dj_cursorx);
	}
#endif
	if (curses_open == 0) {
#ifndef __PDCURSES__ /* don't restart, only refresh: */
		refresh();
#else /* but here, PDCURSES needs an */
		initscr();
#endif
		cbreak(); /* enter c-break mode, don't wait for returns */
		keypad(stdscr, TRUE); /* enable use of cursor and function keys */
		set_mouse_on();
	}
	set_gstat_error_handler(curses_error);
	set_gstat_warning_handler(curses_warning);
	curses_open = 1;
	return;
}

void get_neighbourhood(int key) {
	DATA **dpp, *d;
	double a;

	dpp = get_gstat_data();
	d = dpp[var_ids[0]];
	prompt_for("\rMinimum number in neighbourhood", PR_INT, &(d->sel_min));
	prompt_for("Maximum number in neighbourhood", PR_INT, &(d->sel_max));
	a = d->sel_rad;
	prompt_for("Maximum distance tolerance     ", PR_DOUBLE, &(a));
	if (a != d->sel_rad)
		d->sel_rad = a;
	redraw_scr = 1;
	return;
}

void get_locations(int key) {
	DATA *d;
	int i = 0;
	char s[STRLEN];

	d = get_dataval();
	if (d->id > -1 || get_n_masks() > 0)
		return; /* can't change things like this */
	redraw_scr = 1;
	s[0] = '\0';
	prompt_for("\rEnter mask map name: ", PR_STRING, s);
	if (file_exists(s)) {
		push_mask_name(s);
		return;
	}
	prompt_for("\rSpecify locations as point data? ", PR_YN, &i);
	if (i)
		enter_data(-1);
	return;
}

void set_output(int key) {
	DATA *d;
	const char *cp1, *cp2;
	char **cpp, pr[STRLEN], name[STRLEN];
	int i, j, id;

	clear();
	redraw_scr = 1;
	if (get_n_masks()) {
		cpp = (char **) get_outfile_name();
		for (i = 0; i < get_n_vars(); i++) {
			cp1 = name_identifier(i);
			id = 2 * i;
			sprintf(pr, "predictions(%s)", cp1);
			sprintf(name, "%s", cpp[id] ? cpp[id] : "");
			prompt_for(pr, PR_STRING, name);
			if (name[0] != '\0')
				cpp[id] = string_dup(name);
			id++;
			sprintf(pr, "variances(%s)", cp1);
			sprintf(name, "%s", cpp[id] ? cpp[id] : "");
			prompt_for(pr, PR_STRING, name);
			if (name[0] != '\0')
				cpp[id] = string_dup(name);
			for (j = 0; j < i; j++) {
				cp2 = name_identifier(j);
				id = 2 * get_n_vars() + LTI2(i,j);
				sprintf(pr, "covariances(%s,%s)", cp1, cp2);
				sprintf(name, "%s", cpp[id] ? cpp[id] : "");
				prompt_for(pr, PR_STRING, name);
				if (name[0] != '\0')
					cpp[id] = string_dup(name);
			}
		}
	} else {
		d = get_dataval();
		if (d->id <= -1) {
			msg = "specify prediction locations first";
			return;
		}
		if (o_filename == NULL)
			name[0] = '\0';
		else
			sprintf(name, "%s", o_filename);
		prompt_for("output file name", PR_STRING, name);
		if (name[0] != '\0')
			o_filename = string_dup(name);
	}
	return;
}

void display_map(int key) {
	char cmd[STRLEN];
	const char *name;
	int i, n;

	cmd[0] = '\0';
	clear();
	refresh();
	close_curses();
	switch(what) {
		case XVALID: 
			cross_valid(get_gstat_data());
			break;
		case CONDSIM:
			set_method(GSI);
		case PRED:
			if (get_method() == NSP || get_method() == UIF)
				set_method(get_default_method());
			set_mode();
			check_global_variables();
			printf("[%s]\n", method_string(get_method()));
			predict_all(get_gstat_data());
			strcat(cmd, gl_display);
			for (i = n = 0; i < get_n_outfile(); i++) {
				if ((name = get_outfile_namei(i))) {
					strcat(cmd, " ");
					strcat(cmd, name);
					n++;
				}
			}
			if (n)
				esystem(cmd);
			break;
		case VGM:
			break;
	}
	redraw_scr = 1;
	return;
}

int data_not_valid(DATA *data) {
	if (msg != NULL)
		return 1;
	if (!(data->mode & X_BIT_SET)) {
		msg = "No valid coordinates!";
		return 1;
	}
	if (!(data->mode & V_BIT_SET) || data->n_list <= 1) {
		msg = "No valid data read!";
		return 1;
	}
	return 0;
}

void fix_toggle(VARIOGRAM *v, int fix) {
	int i;
	for (i = 0; i < v->n_models; i++)
		v->part[i].fit_range = fix ? 0 : 1;
	update_variogram(v);
	return;
}

void curses_printlog(const char *s) {
	addstr((char *) s);
}

void set_mouse_on(void) {
#if (HAVE_LIBNCURSES) && (NCURSES_MOUSE_VERSION >= 1)
	mmask_t oldmask;
	mousemask(BUTTON1_PRESSED | BUTTON2_PRESSED | BUTTON3_PRESSED, &oldmask);
#endif
#ifdef __PDCURSES__ 
	mouse_set(BUTTON1_PRESSED | BUTTON2_PRESSED | BUTTON3_PRESSED);
#endif
}

int get_mouse_y(void) {
#if (HAVE_LIBNCURSES) && (NCURSES_MOUSE_VERSION >= 1)
	MEVENT m;
	if (getmouse(&m) == OK)
		return m.y;
#endif
#ifdef __PDCURSES__ 
	request_mouse_pos();
	return MOUSE_Y_POS;
#endif
	return 0;
}

int get_mouse_button(void) {
#if (HAVE_LIBNCURSES) && (NCURSES_MOUSE_VERSION >= 1)
	MEVENT m;
	if (getmouse(&m) == OK) {
		if (m.bstate & BUTTON1_PRESSED)
			return 1;
		if (m.bstate & BUTTON2_PRESSED)
			return 2;
		if (m.bstate & BUTTON3_PRESSED)
			return 3;
	}
#endif
#ifdef __PDCURSES__ 
	request_mouse_pos();
	if (BUTTON_CHANGED(1))
		return 1;
	if (BUTTON_CHANGED(2))
		return 2;
	if (BUTTON_CHANGED(3))
		return 3;
#endif
	/* don't know: */
	return 1;
}

#endif /* else part !define(CURSES) && !defined(NCURSES) */
