#include "utils.h"
#include "constants.h"

#include <string.h>
#include <time.h>

// max value of an array with length given
double maxval(const double *head, const int len)
{
    double max = - PI2;
    for (int i = 0; i < len; i++)
        if (head[i] > max)
            max = head[i];
    return max;
}

// min value of an array with length given
double minval(const double *head, const int len)
{
    double min = PI2;
    for (int i = 0; i < len; i++)
        if (head[i] <  min)
            min = head[i];
    return min;
}

// string equals
bool streqls(const char *dst, const char *src)
{
    return (strcmp(dst, src) == 0);
}

// part string equals
bool streqls(const char *dst, const char *src, int dststart, int dstend)
{
    int dstlen = dstend - dststart + 1;
    char * dst_part = new char[dstlen + 1];
    for (int i = 0; i < dstlen; i++)
    {
        dst_part[i] = dst[dststart + i];
    }
    dst_part[dstlen] = '\0';
    return (strcmp(dst_part, src) == 0);
}

// string trim
char* trim(char *str)
{
    int len = strlen(str);  // string length
    int count = 0;
    for (int i = len - 1; i > -1; i--)
    {
        if (str[i] == ' ')
            count ++;
        else
            break;
    }
    if (count != 0)
        str[len - count] = '\0';
    return str;
}

// string trim len
int len_trim(char *str)
{
    return strlen(trim(str));
}
// fetch date string; datestr length must exceed 10
void sysdate(char *datestr)
{
    // directly use libary function is harmful, because _strdate is for MS compiler only
    //return _strdate(datestr);
    time_t t;
    struct tm * a;
    time(&t);
    a = localtime(&t);
    char day[2], mon[2], year[4];
    day[0] = a->tm_mday / 10;
    day[1] = a->tm_mday % 10;
    mon[0] = a->tm_mon / 10;
    mon[1] = a->tm_mon % 10;
    year[3] = (a->tm_year + 1900) % 10;
    year[2] = ((a->tm_year + 1900) / 10 ) % 10;
    year[1] = ((a->tm_year + 1900) / 100) % 10;
    year[0] = a->tm_year / 1000;
    strcpy(datestr, mon);
    strcat(datestr, "-");
    strcat(datestr, day);
    strcat(datestr, "-");
    strcat(datestr, year);
}

// fetch time string; timestr length must exceed 8
void systime(char *timestr)
{
    // directly use libary function is harmful, because _strtime is for MS compiler only
    //return _strtime(timestr);
    time_t t;
    struct tm *a;
    time(&t);
    a = localtime(&t);
    char hour[2], min[2], sec[2];
    hour[0] = a->tm_hour / 10;
    hour[1] = a->tm_hour % 10;
    min[0] = a->tm_hour / 10;
    min[1] = a->tm_hour % 10;
    sec[0] = a->tm_sec / 10;
    sec[1] = a->tm_sec % 10;
    strcpy(timestr, hour);
    strcat(timestr, ":");
    strcat(timestr, min);
    strcat(timestr, ":");
    strcat(timestr, sec);
}
