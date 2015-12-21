#ifndef EXT_DBASE_INCLUDED
#define EXT_DBASE_INCLUDED

#define EXT_DBASE_FNAME_SIG "extdbase://"

extern void read_ext_dbase(DATA *d);
extern int select_ext_dbase(DATA *d, const DPOINT *whereAbs);
extern void unlink_ext_dbase(DATA *d);

#endif
