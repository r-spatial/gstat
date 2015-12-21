/* defs.h.  Generated automatically by configure.  */
/*
 * defs.h, (c) Edzer J. Pebesma.
 * provides all platform/namespace specific definitions.
 */

#ifndef DEFS_H
#define DEFS_H /* avoid multiple inclusion */

#ifndef NDEBUG
# define NDEBUG /* turns off assert()ions */
#endif 

#include "config.h"

#ifdef USING_R
/* # include <R.h> */
# define exit(n) Rf_error("exiting with code %d", n)
# define printf Rprintf
#else
# define Rprintf printf
#endif

#ifdef SPLUS6WIN32
# define CDECL __cdecl
#else
# define CDECL /* empty */
#endif

#ifdef DMALLOC
#define efree free
#define emalloc malloc
#define ecalloc calloc
#define erealloc realloc
#endif

#ifndef __sgi /* SGI cc works, even without STDC_HEADERS detected! */
# ifndef STDC_HEADERS
# error configure could not find an ansi-c compiler (required)
# endif
#endif

#define GSTAT_NAME      "gstat"
#define GSTAT_CR        "Copyright (C) 1992, 2010 Edzer Pebesma and others"
#define GSTAT_EMAIL     "edzer.pebesma@uni-muenster.de"
#define GSTAT_INFO      "geostatitics@52north.org (subscription required)"
#define GSTAT_ANNOUNCE  "gstat-announce@geo.uu.nl"
#define GSTAT_HOME      "http://www.gstat.org/"
#define USAGE           "usage: gstat [options] [file [file ...]]"

#ifndef GSTAT_OS /* try autodetection of the platform: */
# if defined(__DJGPP) && (__DJGPP==2)
#  define GSTAT_OS "dos/dpmi"
# elif defined(__GO32) /* DJGPP 1.x */
#  define GSTAT_OS "dos/go32"
# elif defined(__CYGWIN__) || defined(CYGWIN)
#  ifndef WIN32
#   define WIN32
#  endif
#  define GSTAT_OS "Win32/Cygwin"
# elif defined (_MSC_VER) /* some Microsoft C version, might work? */
#  define GSTAT_OS "Win32/msc"
#  define SEGMENTED
# elif defined (__MINGW32__)
#  define GSTAT_OS "Win32/MinGW"
# elif defined (BORLANDC)
#  define GSTAT_OS "Win32/bcc"
# elif defined (WIN32) /* NT/98/9x/? */
#  define GSTAT_OS "Win32/unknown"
# elif defined (_AIX)
#  define GSTAT_OS "AIX"
# elif defined (__hpux) 
#  define GSTAT_OS "HP-UX"
# elif defined (__sgi) 
#  define GSTAT_OS "SGI"
# elif defined (__linux)
#  define GSTAT_OS "Linux"
# elif defined (sparc)
#  define GSTAT_OS "Sparc"
# elif defined (__alpha)
#  define GSTAT_OS "DEC/Alpha"
# elif defined (GSTAT_UNIX)
#  define GSTAT_OS "unix"
# else
#  define GSTAT_OS "unknown" /* not important */
# endif
#endif /* ifndef GSTAT_OS */

#if (__DJGPP==2) && !defined(HAVE_UNISTD_H)
# define HAVE_UNISTD_H
#endif

#define GSTATRC    "GSTATRC" /* env. var. holding the initialisation file */
#define HOMERCFILE ".gstatrc" /* rc file in home dir */

/*
 * several buffer sizes 
 */
#ifdef SEGMENTED /* segmented memory: use small buffers */
# define MAX_DATA 64 /* not a maximum, but an increment step size */
#else
# define MAX_DATA 1250 /* not a maximum, but an increment step size */
#endif
#define MAX_ID_LENGTH 40
#define INIT_N_VGMM 2
/* 
 * (for glvars.c:) something, not bigger than 127 
 * because of user interface (crazy though)
 */
#define ERROR_BUFFER_SIZE 1280

#define NOWARRANTY \
"  Because the program is licensed free of charge, there is no warranty\n\
for the program, to the extent permitted by applicable law.  Except when\n\
otherwise stated in writing the copyright holders and/or other parties\n\
provide the program \"as is\" without warranty of any kind, either expressed\n\
or implied, including, but not limited to, the implied warranties of\n\
merchantability and fitness for a particular purpose.  The entire risk as\n\
to the quality and performance of the program is with you.  Should the\n\
program prove defective, you assume the cost of all necessary servicing,\n\
repair or correction.\n\
\n\
  In no event unless required by applicable law or agreed to in writing\n\
will any copyright holder, or any other party who may modify and/or\n\
redistribute the program as permitted above, be liable to you for damages,\n\
including any general, special, incidental or consequential damages arising\n\
out of the use or inability to use the program (including but not limited\n\
to loss of data or data being rendered inaccurate or losses sustained by\n\
you or third parties or a failure of the program to operate with any other\n\
programs), even if such holder or other party has been advised of the\n\
possibility of such damages.\n"

#define COPYRIGHT "Copyright 1992-2006 (C) Edzer J. Pebesma\n\
\n\
This program is free software; you can redistribute it and/or modify\n\
it under the terms of the GNU General Public License as published by\n\
the Free Software Foundation; either version 2 of the License, or\n\
(at your option) any later version.\n\
\n\
This program is distributed in the hope that it will be useful,\n\
but WITHOUT ANY WARRANTY; without even the implied warranty of\n\
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n\
GNU General Public License for more details.\n\
\n\
You should have received a copy of the GNU General Public License\n\
along with this program; if not, write to the Free Software\n\
Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.\n"

#endif /* DEFS_H */
