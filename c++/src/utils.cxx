#include "utils.h"
#include "great_circle_base.h"

#include <string.h>
#include <time.h>
#include <math.h>

// max value of an array with length given
double maxval(const double *head, const int len)
{
    double max = - Constant::PI2;
    for (int i = 0; i < len; i++)
        if (head[i] > max)
            max = head[i];
    return max;
}

// min value of an array with length given
double minval(const double *head, const int len)
{
    double min = Constant::PI2;
    for (int i = 0; i < len; i++)
        if (head[i] <  min)
            min = head[i];
    return min;
}

// make sure longitude lies in bottom-top
void check_longitude(double &longitude, double bottom, double top)
{
    if (longitude > top)
        longitude -= Constant::PI2;
    else if (longitude < bottom)
        longitude += Constant::PI2;
}

// make sure longitude lies in 0-Constant::PI2
void check_longitude(double &longitude)
{
    check_longitude(longitude, -Constant::PI, Constant::PI);
}

// make sure all longitudes lies in bottom-top
void check_longitude_all(double *lon, int size, double bottom, double top)
{
    for (int i = 0; i < size; i++)
        check_longitude(lon[i], bottom, top);
}
// make sure all longitudes lies in 0-Constant::PI2
void check_longitude_all(double *lon, int size)
{
    check_longitude_all(lon, size, Constant::ZERO, Constant::PI2);
}

// make sure latitude lies in -Constant::PIH-Constant::PIH
void check_latitude(double &latitude)
{
/*
    if (zero(latitude - Constant::PIH))
        latitude = Constant::PIH;
    else if (zero(latitude + Constant::PIH))
        latitude = -Constant::PIH;
*/

    if (latitude > Constant::PIH)
        latitude = Constant::PIH;
    if (latitude < - Constant::PIH)
        latitude = - Constant::PIH;

}

// make sure latitudes lie in -Constant::PIH-Constant::PIH
void check_latitude_all(double *lat, int size)
{
    for (int i = 0; i < size; i++)
        check_latitude(lat[i]);
}

// compare two longitudes based on West-->East Increase Direction
bool le(double lon1, double lon2)
{
    return lessequal(lon1, lon2);
}

bool ge(double lon1, double lon2)
{
    return greaterequal(lon1, lon2);
}

bool lessequal(double lon1, double lon2)
{
    if (lon1 == LON_MIN_VALUE || lon2 == LON_MAX_VALUE)
        return true;
    if (lon1 == LON_MAX_VALUE || lon2 == LON_MIN_VALUE)
        return false;
    double delta = lon2 - lon1;
    if (zero(delta))
        return true;
    while (delta > Constant::PI2 || delta < Constant::ZERO)
    {
        if (delta > Constant::PI2)
            delta -= Constant::PI2;
        if (delta < Constant::ZERO)
            delta += Constant::PI2;
    }
    // now delta lies in (0--Constant::PI2)
    if (delta < Constant::PIH)
        return true;
    else
        return false;
}

bool greaterequal(double lon1, double lon2)
{
    if (lon1 == LON_MIN_VALUE || lon2 == LON_MAX_VALUE)
        return false;
    if (lon1 == LON_MAX_VALUE || lon2 == LON_MIN_VALUE)
        return true;
    double delta = lon1 - lon2;
    if (zero(delta))
        return true;
    while (delta > Constant::PI2 || delta < Constant::ZERO)
    {
        if (delta > Constant::PI2)
            delta -= Constant::PI2;
        if (delta < Constant::ZERO)
            delta += Constant::PI2;
    }
    // now delta lies in (0--Constant::PI2)
    if (delta < Constant::PIH)
        return true;
    else
        return false;
}

/* String Functions */

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


/* Great Circle related Functions */
/*
// calculate normal vector of great circle generated by two point
void calc_norm_vec(double *norm_vec, double beglat, double beglon, double endlat, double endlon)
{
    double begx = cos(beglat) * cos(beglon);
    double begy = cos(beglat) * sin(beglon);
    double begz = sin(beglat);
    double endx = cos(endlat) * cos(endlon);
    double endy = cos(endlat) * sin(endlon);
    double endz = sin(endlat);

    norm_vec[0] = begy * endz - endy * begz;
    norm_vec[1] = begz * endx - endz * begx;
    norm_vec[2] = begx * endy - endx * begy;
}

// normalize normal vector if required
void do_normalize_vec(double *vec, int dim)
{
    double norm2 = 0.0;
    for (int i = 0; i < dim; i++)
        norm2 += vec[i] * vec[i];
    norm2 = sqrt(norm2);
    for (int i = 0; i < dim; i++)
        vec[i] /= norm2;
}
void do_normalize_vec(double *vec)
{
    do_normalize_vec(vec, 3);
}

// calculate intersection of greate circle and latitude segment
void calc_gc_with_lat_segment(double &intrsct_lon, double *norm_vec, double lat, double lon_low, double lon_high)
{
    double z0 = sin(lat);       // latitude value of latitude segment
    double r = sqrt(1 - z0*z0); // radius of latitude plane

    // trianglar formula Acos(x) + Bsin(x) = sqrt(A^2+B^2)cos(x + phi)
    double R = sqrt(norm_vec[0] * norm_vec[0] + norm_vec[1] * norm_vec[1]); // sqrt(A^2+B^2)
    double theta0 = atan2(norm_vec[1], norm_vec[0]);    // phi

    check_longitude(theta0, Constant::ZERO, Constant::PI2);                 // adjust range from -Constant::PI--Constant::PI to Constant::ZERO--Constant::PI2
    intrsct_lon = acos(-norm_vec[2]*z0/r/R) + theta0;
}

// calculate intersection of great circle and longitude segment
void calc_gc_with_lon_segment(double &intrsct_lat, double *norm_vec, double lon, double lat_low, double lat_high)
{
}

// calculate intersection of great circle and great circle segment

void calc_gc_with_gc_segment(Point * intrsct, Segment * beg, Segment * end)
{
}
*/
