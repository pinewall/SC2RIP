#include "netcdf_utils.h"

NCUtils::NCUtils(char * file)
{
    istat = nc_open(file, NC_NOWRITE, &fileid);
    netcdf_error_handler(istat);
}

NCUtils::~NCUtils()
{
    istat = nc_close(fileid);
    netcdf_error_handler(istat);
}

void NCUtils::netcdf_error_handler()
{
    if (istat != NC_NOERR)
    {
        printf("NETCDF Error: %s\n", nc_strerror(istat));
        exit(-1);
    }
}

void NCUtils::fetch_att_text(char *dst, char * varname, char *attr_name)
{
    if (varname == 0)
    {
        istat = nc_get_att_text(fileid, NC_GLOBAL, attr_name, dst);
        netcdf_error_handler(istat);
    }
    else
    {
        istat = nc_inq_varid(fileid, varname, &varid);
        netcdf_error_handler(istat);
        istat = nc_get_att_text(fileid, varid, attr_name, dst);
        netcdf_error_handler(istat);
    }
}

void NCUtils::fetch_dim(int & dst, char * dim_name)
{
    istat = nc_inq_dimid(fileid, dim_name, &varid);
    netcdf_error_handler(istat);
    istat = nc_inq_dimlen(fileid, varid, dst);
    netcdf_error_handler(istat);
}

void NCUtils::fetch_var_int(int * dst, char * var_name)
{
    istat = nc_inq_varid(fileid, var_name, &varid);
    netcdf_error_handler(istat);
    istat = nc_get_var_int(fileid, varid, dst);
    netcdf_error_handler(istat);
}

void NCUtils::fetch_var_double(double * dst, char * var_name)
{
    istat = nc_inq_varid(fileid, var_name, &varid);
    netcdf_error_handler(istat);
    istat = nc_get_var_double(fileid, varid, dst);
    netcdf_error_handler(istat);
}
