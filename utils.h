#ifndef _UTILS_H_
#define _UTILS_H_ 1

#include "constants.h"

#define EPS (1e-10)
/* Numerical Functions */

// max|min value of an array with length given
double maxval(const double *head, const int len);
double minval(const double *head, const int len);
void check_longitude(double &longitude, double bottom, double top);
void check_longitude(double &longitude);
void check_latitude(double &latitude);
void check_longitude_all(double *lon, int size, double bottom, double top);
void check_longitude_all(double *lon, int size);
void check_latitude_all(double *lat, int size);

// return minimal integer
inline int MIN(int a, int b) {  return ((a > b)?b:a); };
// return maximize integer
inline int MAX(int a, int b) {  return ((a > b)?a:b); };
// return absolute float
inline double ABS(double d) {   return ((d > 0)?(d):(-d)); };
// check whether a number is zero (nearly for float)
inline bool zero(double num) {  return (ABS(num) < 1e-30); };
// check whether a number is nonzero
inline bool nonzero(double num) {   return !zero(num); };

// compare two longitudes based on West-->East Increase Direction
bool lessequal(double lon1, double lon2);
bool greaterequal(double lon1, double lon2);
bool le(double lon1, double lon2);
bool ge(double lon1, double lon2);

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
