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
    bool ret = strcmp(dst_part, src) == 0;
    delete [] dst_part;
    return ret;
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
    day[0] = a->tm_mday / 10 + '0';             // day from 1 to 31
    day[1] = a->tm_mday % 10 + '0';
    mon[0] = (a->tm_mon + 1) / 10 + '0';        // month from 0 to 11
    mon[1] = (a->tm_mon + 1)% 10 + '0';
    year[3] = (a->tm_year + 1900) % 10 + '0';
    year[2] = ((a->tm_year + 1900) / 10 ) % 10 + '0';   // year since 1900
    year[1] = ((a->tm_year + 1900) / 100) % 10 + '0';
    year[0] = (a->tm_year + 1900) / 1000 + '0';
    int index = 0;
    datestr[index++] = mon[0];
    datestr[index++] = mon[1];
    datestr[index++] = '-';
    datestr[index++] = day[0];
    datestr[index++] = day[1];
    datestr[index++] = '-';
    datestr[index++] = year[0];
    datestr[index++] = year[1];
    datestr[index++] = year[2];
    datestr[index++] = year[3];
    datestr[index++] = '\0';
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
