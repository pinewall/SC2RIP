#ifndef _NETCDF_ERROR_H_
#define _NETCDF_ERROR_H_ 1

#include "kinds.h"
#include "constants.h"

#include <netcdf.h>
#include <stdio.h>
#include <stdlib.h>

extern "C"
extern void netcdf_error_handler(int istat);

#endif
