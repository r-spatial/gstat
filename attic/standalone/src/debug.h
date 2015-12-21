#ifndef DEBUG_H
#	define DEBUG_H /* avoid multiple inclusion */

extern int debug_level;

/*
 * DEBUG macro's:
 */
#define DB_HELP  (-1) /* print debug help */
#define DB_SILENT (0)
#define DB_NORMAL (1 << 0)
#define DB_DUMP   (1 << 1) /* dump global variabels */
#define DB_FIT    (1 << 2) /* fit diagnostics */
#define DB_DATA   (1 << 3) /* drop data */
#define DB_SEL    (1 << 4) /* drop selection */
#define DB_COV    (1 << 5) /* drop covariances */
#define DB_ORDER  (1 << 6) /* order relation violation */
#define DB_FORCE  (1 << 7) /* print warning if neighbourhood selection */
#define DB_TRACE  (1 << 8) /* print numbers */
#define DB_BLOCK  (1 << 9) /* block discretization diagnostics (data) */

extern void printlog(const char *fmt, ...);
#define DUMP(a); {if(debug_level & DB_DUMP) { printlog("%s", a); }}
#define DEBUG_HELP    (debug_level & DB_HELP)
#define DEBUG_SILENT  (debug_level == DB_SILENT)
#define DEBUG_NORMAL  (debug_level & DB_NORMAL)
#define DEBUG_DUMP    (debug_level & DB_DUMP)
#define DEBUG_FIT     (debug_level & DB_FIT)
#define DEBUG_DATA    (debug_level & DB_DATA)
#define DEBUG_SEL     (debug_level & DB_SEL)
#define DEBUG_COV     (debug_level & DB_COV)
#define DEBUG_ORDER   (debug_level & DB_ORDER)
#define DEBUG_VGMFIT  (debug_level & DB_ORDER)
#define DEBUG_FORCE   (debug_level & DB_FORCE)
#define DEBUG_TRACE   (debug_level & DB_TRACE)
#define DEBUG_BLOCK   (debug_level & DB_BLOCK)

#define DEBUG_OPTIONS "\
  #  gstat debug option values:\n\
  0: no output, be silent (same as -s)\n\
  1: normal output (default value)\n\
  2: print all global variables and extended error messages\n\
  4: print OLS and WLS fit diagnostics\n\
  8: print all data\n\
 16: print every neighbourhood selection\n\
 32: print all covariance matrices, solutions, design matrices etc.\n\
 64: print variogram fit diagnostics and order relation violations\n\
128: print warning on forced neighbourhoods\n\
256: print current row,column or record number\n\
512: print block discretization points (data)\n\
to combine options, sum their values -- 1023 invokes them all\n"

#endif /* DEBUG_H */
