/* parse.y & lex.l */

#if defined(__cplusplus)
extern "C" {
#endif

/* parse.y: */
int parse_file(const char *fileName);
int parse_cmd(const char *cmd, const char *fname);
void parse_gstatrc(void);
#ifdef VARIO_H
int read_variogram(VARIOGRAM *v, const char *source);
#endif

#ifdef DATA_H
int read_vector(D_VECTOR *d, char *fname);
#endif

/* lex.l: */
int yylex(void);
void lex_error(void);
void set_lex_source(const char *source, const char *fname);

#if defined(__cplusplus)
}
#endif
