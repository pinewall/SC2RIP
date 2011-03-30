#ifndef _UTILS_
#define _UTILS_ 1


#define MAX(a, b) ((a > b) ? a : b)
#define MIN(a, b) ((a < b) ? a : b)

// max|min value of an array with length given
double maxval(const double *head, const int len);
double minval(const double *head, const int len);

// string equals
bool streqls(const char *dst, const char *src);

// fetch date string
void sysdate(char *datestr);

// fetch time string
void systime(char *timestr);

#endif
