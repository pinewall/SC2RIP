#ifndef _NETCDF_UTILS_H_
#define _NETCDF_UTILS_H_ 1

#include <netcdf.h>

class NCUtils
{
    private:
        int fileid;     // file id from nc_open
        int varid;      // variable id used in variable-related functions
        int istat;      // status variable storing return value

    public:
        NCUtils(char * file);
        ~NCUtils();
        void netcdf_error_handler(int istat);
        void fetch_global_attr_text(char * dst, char * attr_name);
        void fetch_dim(int & dst, char * dim_name);
        void fetch_var_int(int * dst, char * var_name);
        void fetch_var_double(double * dst, char * var_name);
};

#endif
