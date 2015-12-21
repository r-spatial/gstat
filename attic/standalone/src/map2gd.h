#ifndef MAP2GD_H
# define MAP2GD_H

typedef struct {
	int every, n;
	char **entries;
} TICKS;

int map2gd(int argc, char *argv[]);
int one_map2gd(GRIDMAP *map_pointer, const char *f_name,
	TICKS *xticks, TICKS *ytics, int legend);

#endif
