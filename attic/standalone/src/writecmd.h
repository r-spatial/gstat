#ifndef WRITECMD_H
#define WRITECMD_H

#if defined(__cplusplus)
extern "C" {
#endif

void logprint_cmd(void);
void fprint_cmd(FILE *f);

const char *sprint_cmd(void);
const char *sprint_glvars(int anyway);

#if defined(__cplusplus)
}
#endif

#endif
