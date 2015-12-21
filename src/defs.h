/*
 * provides some definitions
 */

#ifndef DEFS_H
#define DEFS_H /* avoid multiple inclusion */

#include <assert.h> /* but assertions are off, by default */

#define CDECL /* empty */

/*
 * several buffer sizes 
 */
#define MAX_DATA 1250 /* not a maximum, but an increment step size */
#define INIT_N_VGMM 2
/* 
 * (for glvars.c:) something, not bigger than 127 
 * because of user interface (crazy though)
 */
#define ERROR_BUFFER_SIZE 1280

#endif /* DEFS_H */
