#ifndef MAPIO_H /* avoid multiple inclusion */
#define MAPIO_H

#ifdef HAVE_LIBGDAL
# include "gdal.h"
#endif

/*! \file mapio.h 
	\brief functions for reading and writing grid map in several formats
*/

/*! type of grid map format */
typedef enum {
	MT_UNKNOWN = 0,
	MT_CSF,      /*!< <A HREF="http://www.geog.uu.nl/pcraster.html">PCRaster</A>
				format (binary), API in csf directory */
	MT_GDAL,     /*!< <A HREF="http://remotesensing.org/gdal">GDAL</A> format */
	MT_ARCGRID,  /*!< <A HREF="http://www.esri.com/">ArcInfo</A> gridfile, 
				gridascii (ascii) or gridfloat (binary) */	
	MT_IDRISI,   /*!< <A HREF="http://www.clarklabs.org>Idrisi</A> .img file,
				ascii or binary */	                        
	MT_IDRISI32, /*!< <A HREF="http://www.clarklabs.org>Idrisi 32</A>image file,
				ascii or binary */	                        
	MT_GNUPLOT,  /*!< binary gnuplot 2d matrix format (cannot handle MV's) */
	MT_T2,       /*!< MIKE SHE ``T2'' 2d grid map format (www.dhi.dk) */
	MT_ERMAPPER, /*!< <A HREF="http://www.ermapper.com">ER-Mapper</A> 
				V5.0+ dataset with a single image channel */	    
	MT_GRASS,    /*!< <A HREF="http://www.geog.uni-hannover.de/grass/">
				Grass</A> raster format, uses grass' gis library */	
	MT_GMT,      /*!< <A HREF="http://www.soest.hawaii.edu/soest/gmt3.0.html">
				GMT</A> Grid format, uses NetCDF library */	
	MT_SURFER,   /*!< <A HREF="http://www.golden.com">Surfer</A> DSAA 
				(ascii grid) format */
	MT_GSLIB     /*!< <A HREF="http://www.gslib.com/">GSLIB</A> grid format 
				(only a 2D subset) */
} MAPTYPE;

/*! ER-Mapper celltypes used for read/write_binary */
typedef enum {	 
	CT_NONE			= 0,
	CT_UNKNOWN		= 1,
	CT_UINT8		= 2,
	CT_UINT16		= 3,
	CT_UINT32		= 4,
	CT_INT8			= 5,
	CT_INT16		= 6,
	CT_INT32		= 7,
	CT_IEEE4		= 8,
	CT_IEEE8		= 9
} CellType;

#define UINT8 unsigned char
#define UINT16 unsigned short int
#define UINT32 unsigned long int
#define INT8  char
#define INT16 short int
#define INT32 long int
#define IEEE4 float
#define IEEE8 double

typedef enum {
		READ_ONLY = 0,
		WRITE_ONLY = 1
} MAP_READ_STATUS;

/*! structure to hold grid map information */
typedef struct gridmap {
	MAPTYPE type;         /*!< type of grid map */
	MAP_READ_STATUS status; /*! READ_ONLY or WRITE_ONLY */
	const char *filename; /*!< name (or base name) of grid map */
	char *history,        /*!< only used for CSF maps */
		*description;     /*!< only used for some maps */
	unsigned int rows,    /*!< number of rows in map */
		cols,             /*!< number of colums in map */
		base_size,        /*!< size of malloced area (cells), in case of blocked allocation */
		current_row,      /*!< current row, if reading/writing is row-wise */
		*first_time_row;  /*!< have we been at this row before? */
	CellType celltype;    /*!< cell type */
	int is_binary,        /*!< is this a binary map format */
		swap_buf;         /*!< swap contents of base buffer? */
	double x_ul,          /*!< x-coordinate upper left corner of map area */
		y_ul,             /*!< y-coordinate of upper left corner of map area */
		cellsizex,        /*!< size of grid cells in x-direction */
		cellsizey;        /*!< size of grid cells in y-direction */
	float cellmin,        /*!< minimum value of grid map */
		cellmax,          /*!< maximum value of grid map */
		misval;           /*!< missing value flag (if present) */
	float **grid,         /*!< 2d matrix holding the values (pointer array) */
		*base;            /*!< base pointer to malloc'ed or mmap'ed area */
	void *CSF_MAP;        /*!< cast to MAP * */
#ifdef HAVE_LIBGDAL
	GDALDatasetH  hDataset;
	GDALDriverH   hDriver;
	double *GeoTransform;
#endif
	struct gridmap * (*write)(struct gridmap *m); /*!< write & close a map */
	void (*read_row)(struct gridmap *m, float *buf, unsigned int row);
	void (*write_row)(struct gridmap *m, float *buf, unsigned int row);
	/* see mapio.c remark around #define CSF_ROWWISE for the last two */
} GRIDMAP;

#define SQUARECELLSIZE(map) ((map->cellsizex != map->cellsizey ? \
	ErrMsg(ER_IMPOSVAL, "cannot deal with non-square cells"), 0.0 : \
	map->cellsizex))

GRIDMAP *map_read(GRIDMAP *m);
int map_cell_is_mv(GRIDMAP *m, unsigned int row, unsigned int col);
float map_get_cell(GRIDMAP *m, unsigned int row, unsigned int col);
int map_put_cell(GRIDMAP *m, unsigned int row, unsigned int col, float value);
int map_xy2rowcol(GRIDMAP *m, double x, double y, 
				unsigned int *row, unsigned int *col);
int map_rowcol2xy(GRIDMAP *m, 
				unsigned int row, unsigned int col, double *x, double *y);
GRIDMAP *map_dup(const char *fname, GRIDMAP *m);
GRIDMAP *new_map(MAP_READ_STATUS status);
void map_free(GRIDMAP *m);

GRIDMAP *map_switch_type(GRIDMAP *in, MAPTYPE type);

void map_name_nr(GRIDMAP *mask, const char *base, char *name, int nr, int max);
int map_equal(GRIDMAP *a, GRIDMAP *b);

#endif
