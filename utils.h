#ifndef _UTILS_
#define _UTILS_ 1

#include <string.h>
#include "constants.h"

#define MAX(a, b) ((a > b) ? a : b)
#define MIN(a, b) ((a < b) ? a : b)

double maxval(const double *head, const int len);
double minval(const double *head, const int len);

bool streqls(const char *dst, const char *src);

#endif
