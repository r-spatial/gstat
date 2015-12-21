typedef enum {
	NOTSET = 0,
	GYR = 1, /* green-yellow-red */
	RYG, /* red-yellow-green */
	UNIRAS, /* uniras, rainbow */
	BPY, /* blue-pink-yellow */
	GYR0, /* `old' GYR scale */
	BW, 
	WB,
	USER_DEFINED
} PALETTE;

int palet(int argc, char *argv[]);
float *GetRGB(PALETTE palet, float pos, float maxcol);
PALETTE int2pal(int ipal);
