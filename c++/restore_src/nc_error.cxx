#include "nc_error.h"
extern "C"
void netcdf_error_handler(int istat)
{
    if (istat != NC_NOERR)
    {
        printf("NETCDF Error: %s\n", nc_strerror(istat));
        exit(-1);
    }
}
