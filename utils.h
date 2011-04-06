#ifndef _UTILS_
#define _UTILS_ 1


#define MAX(a, b) ((a > b) ? a : b)
#define MIN(a, b) ((a < b) ? a : b)
#define MOD(a, b) (a % b)
/* Numerical Functions */

// max|min value of an array with length given
double maxval(const double *head, const int len);
double minval(const double *head, const int len);

/* String Functions */

// string equals
bool streqls(const char *dst, const char *src);

// part string equals
bool streqls(const char *dst, const char *src, int dststart, int dstend);

// string trim
char* trim(char *str);

// string trim len
int len_trim(char *str);

// fetch date string
void sysdate(char *datestr);

// fetch time string
void systime(char *timestr);

#endif
