#include "utils.h"

double maxval(const double *head, const int len)
{
    double max = - PI2;
    for (int i = 0; i < len; i++)
        if (head[i] > max)
            max = head[i];
    return max;
}

double minval(const double *head, const int len)
{
    double min = PI2;
    for (int i = 0; i < len; i++)
        if (head[i] <  min)
            min = head[i];
    return min;
}

bool streqls(const char *dst, const char *src)
{
    return (strcmp(dst, src) == 0);
}
