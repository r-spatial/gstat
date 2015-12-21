/* (c) cees wesseling; this is an interface to
external code, not provided with gstat open source.
*/
#ifndef INCLUDED_DEFS
#include "defs.h"
#define INCLUDED_DEFS
#endif

#ifdef HAVE_EXT_DBASE
 /* up to end else compile to 0 object */


#ifndef INCLUDED_UTILS
#include "utils.h"
#define INCLUDED_UTILS
#endif
#ifndef INCLUDED_USERIO
#include "userio.h"
#define INCLUDED_USERIO
#endif
#ifndef INCLUDED_DATA
#include "data.h"
#define INCLUDED_DATA
#endif

#ifndef INCLUDED_EXT_DBASE
#include "ext_dbase.h"
#define INCLUDED_EXT_DBASE
#endif

#ifndef INCLUDED_DEBUG
#include "debug.h"
#define INCLUDED_DEBUG
#endif

#ifndef INCLUDED_PQDB_CAPI
#include "pqdb_capi.h"
#define INCLUDED_PQDB_CAPI
#endif

typedef struct {
  float x,y,v;
} XYV;

typedef struct {
  PQDbase      *db;
  double        offset[3];
} EXTDBASE_LINK;

static double X1=1.0;

static PQDbase* pqdbase(DATA *d) {
    return ((EXTDBASE_LINK *)d->ext_dbase)->db;
}
static double* offset(DATA *d) {
    return ((EXTDBASE_LINK *)d->ext_dbase)->offset;
}

/* clean up */
void unlink_ext_dbase(DATA *d) {
  if (d->ext_dbase) {
    pcr_PQDbaseClose(pqdbase(d));
    efree(d->ext_dbase);
  }
}

static void check_status(DATA *d) {
    if (pcr_PQDbaseError(pqdbase(d)))
      ErrMsg(ER_EXT_DBASE,pcr_PQDbaseErrorMessage(pqdbase(d)));
}


/*!
 * initialise with conditions to be met by points returned
 * from select_at(); dist2 must be set
 * p is ptr to buffer that is owned by ext_dbase module,
 * other gstat code shall not call free on p!
 * not yet used in test setup
 * \todo
 *  see TODO in code
 *
 */
static DPOINT* init_dpoint(
    DPOINT *p,              /* init this one */
    const XYV           *r, /* with this relative value */
    const DPOINT        *where, /* where relative value */
    const double        *offset) /* offset of relative value */
{
  p->x=r->x+offset[0];
  p->y=r->y+offset[1];
  p->z=0;
  p->attr=r->v;
  p->u.dist2=SQR(r->x-where->x)+SQR(r->y-where->y); /* dist2 in relative */

  /* do p->X[0]=1.0, using static buffer
   *  instead of attaching to (double *)emalloc(sizeof(double))
   */
  p->X = &X1;
  /* I am a point not a block */
  SET_POINT(p);
  return p;
}

static void init_data(DATA *d) {
  const char *urlArg = d->fname+strlen(EXT_DBASE_FNAME_SIG);
  int colNrs[3];
  DUMP(d->fname);
  DUMP(": trying extdbase ... ");

  assert(d);

  colNrs[0]=d->colnx-1;
  colNrs[1]=d->colny-1;
  colNrs[2]=d->colnvalue-1;

  d->ext_dbase = (EXTDBASE_LINK *)emalloc(sizeof(EXTDBASE_LINK));
  ((EXTDBASE_LINK *)d->ext_dbase)->db = pcr_PQDbaseOpen(urlArg,colNrs);
  if (!pqdbase(d))
    ErrMsg(ER_EXT_DBASE, "Out of Memory");

  // TODO make this a test, like file not exists
  check_status(d);

  pcr_PQDbaseOffset(pqdbase(d),offset(d));

  d->type.type = DATA_EXT_DBASE;
  d->type.name = "extdbase";

  d->mode = X_BIT_SET | Y_BIT_SET | V_BIT_SET;
  d->n_X = 1;
}


/*
 * \todo see TODO in code
 */
int select_ext_dbase(
    DATA *d,
    const DPOINT *whereAbs) /* where in absolute coords */
{
  DPOINT where=*whereAbs; /* becomes relative to offset */
  int i;
  PQDbaseProximitySearch ps;

  assert(d->type.type==DATA_EXT_DBASE);

  ps.coords[0]= (where.x-=offset(d)[0]);
  ps.coords[1]= (where.y-=offset(d)[1]);
  ps.coords[2]= (where.z-=offset(d)[2]);

  ps.radius= d->sel_rad;
  ps.square= d->square;
  ps.minNr = d->sel_min;
  ps.maxNr = d->sel_max;

  pcr_PQDbaseSearch(pqdbase(d),&ps);
  check_status(d);

  d->n_sel    =ps.nrResult;
  d->n_sel_max=ps.nrResult;

  if(!ps.nrResult)
    return 0;


  d->P_base = (DPOINT *)erealloc(d->P_base,sizeof(DPOINT)*d->n_sel);
  d->sel    = (DPOINT **)erealloc(d->sel, d->n_sel * sizeof(DPOINT *));

  for(i=0; i < ps.nrResult; i++) {
    d->sel[i]= init_dpoint(d->P_base+i,ps.result[i],
                           &where,offset(d));
    assert(ps.recLen==12);
    // compute index by subtracting ptr from base
    SET_INDEX(d->sel[i],((int)((char *)ps.result[i]-(char *)ps.baseOfData))/ps.recLen);
  }
  return d->n_sel;
  /* TODO -1 MORE WORK TO DO ? */
}


/* initing DATA *d
 * zorg dat d->list=NULL
 * d->n_list # waarneming zetten, of i.i.g = 0, niet -1
 *  zodat checks weten dat er data inzit
 * report_data aanpassen, doorkijken wat zinnig is
 * NB. het werkt wel al
 *  i.i.g niet stdDev en mean enzo berekenen.
 *  simpelste is dat pqdb die info kan terug geven
 * in pqdb de stats berekenen opslaan:
 *  min/max X,Y min/max/mean/stddev V
 */
void read_ext_dbase(DATA *d)
{

  int      nrPoints;
  double   minX,  maxX, minY, maxY, meanV, sdV;

  init_data(d);
  d->id=0;

  pcr_PQDbaseInfoXYV(pqdbase(d),&nrPoints,
                     &minX,&maxX,
                     &minY,&maxY,
                     &meanV, &sdV);
  check_status(d);

  d->n_list =nrPoints;
  d->maxX = maxX;
  d->minX = minX;
  d->maxY = maxY;
  d->minY = minY;
  d->mean = meanV;
  d->std  = sdV;

}

/* HAVE_EXT_DBASE */
#endif 
