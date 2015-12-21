#ifndef USERIO_H
#define USERIO_H

#if defined(__cplusplus)
extern "C" {
#endif

enum Gstat_errno {
	ER_NOERROR     =  0 /* no error */,
	ER_NULL        =  1 /* internal error: should not occur */,
	ER_VARNOTSET   =  2 /* a required variable was not set by the user */,
	ER_RANGE       =  3 /* range error (outside permitted values) */,
	ER_IMPOSVAL    =  4 /* a variable was set to an illegal value */,
	ER_NOFILE      =  5 /* no input file specified */,
	ER_WRITE       =  6 /* write error on file */,
	ER_READ        =  7 /* read error on file */,
	ER_RDFLT       =  8 /* error while converting a string to a float */,
	ER_RDINT       =  9 /* error while converting a string to an int */,
	ER_SYNTAX      = 10 /* syntax error */,
	ER_ARGOPT      = 11 /* error in command line option arguments */,
	ER_DOMAIN      = 12 /* math error */,
	ER_MEMORY      = 13 /* memory exhausted */,
	ER_IO          = 14 /* i/o conflict (e.g. redirection not permitted) */,
	ER_NOCMD       = 15 /* no command file specified */,
	ER_NOCURSES    = 16 /* no curses user interface compiled in */,
	ER_PWRITE      = 17 /* error while writing to a pipe */,
	ER_PREAD       = 18 /* error while reading from a pipe */,
	ER_SECURE      = 19 /* secure mode: operation not allowed */,
	ER_MESCHACH    = 20 /* error happened somewhere in meschach matrix lib */,
	ER_EXT_DBASE   = 21 /* error happened somewhere in extdbase.c on hooks */
};

#define MAX_ERRNO 21
extern const char *error_messages[MAX_ERRNO+1];
void message(char *fmt, ...);  /* message() calls always preceed ErrMsg() */
#define ErrMsg(a,b) gstat_error(__FILE__,__LINE__,a,b)
void gstat_error(char *fname, int line, 
	enum Gstat_errno err_nr, const char *msg);
/* command line option error: */
#define ErrClo(a) gstat_clo_error(__FILE__,__LINE__,ER_ARGOPT,a) 
void gstat_clo_error(char *f, int l, enum Gstat_errno err, int a);

void init_userio(int use_stdio);
void pr_warning(char *fmt, ...);
void printlog(const char *fmt, ...);
void print_progress(unsigned int current, unsigned int total);
void print_to_logfile_if_open(const char *mess);

enum Gstat_errno get_gstat_errno(void);
void reset_gstat_errno(void);
void setup_meschach_error_handler(int using_R);

void set_gstat_warning_handler(void (*warning_fn)(const char *message));
void set_gstat_error_handler(void (*error_fn)(const char *message, int level));
void set_gstat_log_handler(void (*logprint)(const char *str));
void set_gstat_progress_handler(
	void (*progress)(unsigned int step, unsigned int total));
void push_gstat_progress_handler(
	void (*progress)(unsigned int step, unsigned int total));
void pop_gstat_progress_handler(void);

void default_warning(const char *mess);
void default_error(const char *mess, int level);
void default_printlog(const char *mess);
void default_progress(unsigned int step, unsigned int total);
void no_progress(unsigned int current, unsigned int total);

int set_gstat_log_file(FILE *f);

const char *get_gstat_error_message(void);

#if defined(__cplusplus)
}
#endif

#endif /* USERIO_H */
