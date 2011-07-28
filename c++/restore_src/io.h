#ifndef _IO_H_
#define _IO_H_ 1

#define     NETCDF_TYPE         "netcdf"

#define     FILE_STATUS_NEW     "new"
#define     FILE_STATUS_OLD     "old"

#include <netcdf.h>

// print error string
void log(const char *error);

// data file IO class
class DATA_IO
{
    protected:
        char    IO_type [64];
        char    filename[256];
        char    status[64];                 // file status

    public:
        char * getIOType();
        char * getFilename();
        char * getFileStatus();
        bool   setIOType        (const char * io_type);
        bool   setFilename      (const char * file_name);
        bool   setFileStatus    (const char * file_status);

};

class NETCDF_IO : public DATA_IO
{
    private:
        int     fileID;                 // file ID returned by open/create
        int     variableID;
        int     dimensionID;
        int     dimensionID_array[64];
        int     dimensionLen_array[64];
        char    data_type[64];
        int     ncstat;                     // used to record last status
    public:
        NETCDF_IO (const char * io_type, const char * file_name, const char * file_status);


        bool    NETCDF_IO_Create       (int mode);
        bool    NETCDF_IO_Open         (int mode);
        bool    NETCDF_IO_Close        (void);

        bool    NETCDF_IO_Read_Dim     (unsigned int & dimension_len, const char * dimension_name);                                                     
        bool    NETCDF_IO_Read_Dim     (unsigned int & dimension_len, int dimension_id);                                                     
        bool    NETCDF_IO_Read_Var     (void ** read_buf, unsigned int ** dimension_len, int  & num_of_dimensions, const char * variable_name);         // read_buf allocated by NETCDF_IO
        bool    NETCDF_IO_Free         (void * read_buf, void * dimension_len);
        bool    NETCDF_IO_Write_Dim    (int   dimension_len, const char * dimension_name);
        bool    NETCDF_IO_Write_Var    (void * write_buf, int * dimension_len, int    num_of_dimensions, const char * variable_name);
};

#endif
