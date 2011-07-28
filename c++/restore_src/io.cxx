#include "io.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

// print error string
void log(const char *error)
{
    printf(error);
    printf("\n");
}


/* class DATA_IO */
char * DATA_IO::getIOType()
{
    return IO_type;
}

char * DATA_IO::getFilename()
{
    return filename;
}

char * DATA_IO::getFileStatus()
{
    return status;
}

bool   DATA_IO::setIOType (const char * io_type)
{
    return (strcpy(IO_type, io_type) != NULL);
}

bool   DATA_IO::setFilename (const char * file_name)
{
    return (strcpy(filename, file_name) != NULL);
}

bool   DATA_IO::setFileStatus (const char * file_status)
{
    return (strcpy(status, file_status) != NULL);
}

/* class NETCDF_IO */
NETCDF_IO::NETCDF_IO (const char * io_type, const char * file_name, const char * file_status)
{
    //assert(setIOType(io_type));
    //assert(setFilename(file_name));
    //assert(setFileStatus(file_status));
    strcpy(IO_type, io_type);
    strcpy(filename, file_name);
    strcpy(status, file_status);
}

bool NETCDF_IO::NETCDF_IO_Create (int  mode)
{
    if (strcmp(status, FILE_STATUS_NEW) == 0)
    {
        ncstat = nc_create(filename, mode, &fileID);
        if (ncstat == NC_NOERR)
        {
            //assert(setFileStatus(FILE_STATUS_OLD));
            strcpy(status, FILE_STATUS_OLD);
            return true;
        }
        else
        {
            fprintf(stderr, "%s\n", nc_strerror(ncstat));
            fprintf(stderr, "%s\n", "@ NETCDF_IO_Create()");
            return false;
        }
    }
    else
    {
        fprintf(stderr, "%s\n", "Cannot Create File that is Exist!");
        fprintf(stderr, "%s\n", "@ NETCDF_IO_Create()");
        return false;
    }

}

bool NETCDF_IO::NETCDF_IO_Open (int mode)
{
    if (strcmp(status, FILE_STATUS_OLD) == 0)
    {
        ncstat = nc_open(filename, mode, &fileID);
        if (ncstat == NC_NOERR)
            return true;
        else
        {
            fprintf(stderr, "%s\n", nc_strerror(ncstat));
            fprintf(stderr, "%s\n", "@ NETCDF_IO_Open()");
            return false;
        }
    }
    else
    {
        fprintf(stderr, "%s\n", "Cannot Open File that is not Exist!");
        fprintf(stderr, "%s\n", "@ NETCDF_IO_Open()");
        return false;
    }
}

bool NETCDF_IO::NETCDF_IO_Close (void)
{
    if (strcmp(status, FILE_STATUS_OLD) == 0)
    {
        ncstat = nc_close(fileID);
        if (ncstat == NC_NOERR)
            return true;
        else
        {
            fprintf(stderr, "%s\n", nc_strerror(ncstat));
            fprintf(stderr, "%s\n", "@ NETCDF_IO_Close()");
            return false;
        }
    }
    else
    {
        fprintf(stderr, "%s\n", "Cannot Close File that is not Exist!");
        fprintf(stderr, "%s\n", "@ NETCDF_IO_Close()");
        return false;
    }
}

bool NETCDF_IO::NETCDF_IO_Read_Dim (unsigned int & dimension_len, const char * dimension_name)
{
    ncstat = nc_inq_dimid (fileID, dimension_name, & dimensionID);
    if (ncstat != NC_NOERR)
    {
        fprintf(stderr, "%s\n", nc_strerror(ncstat));
        fprintf(stderr, "%s\n", "@ NETCDF_IO_Read_Dim()");
        return false;
    }
    ncstat = nc_inq_dimlen (fileID, dimensionID, &dimension_len);
    if (ncstat != NC_NOERR)
    {
        fprintf(stderr, "%s\n", nc_strerror(ncstat));
        fprintf(stderr, "%s\n", "@ NETCDF_IO_Read_Dim()");
        return false;
    }
    return true;
}

bool NETCDF_IO::NETCDF_IO_Read_Dim (unsigned int & dimension_len, int dimension_id)
{
    ncstat = nc_inq_dimlen (fileID, dimension_id, &dimension_len);
    if (ncstat != NC_NOERR)
    {
        fprintf(stderr, "%s\n", nc_strerror(ncstat));
        fprintf(stderr, "%s\n", "@ NETCDF_IO_Read_Dim()");
        return false;
    }
    return true;
}

bool NETCDF_IO::NETCDF_IO_Read_Var (void ** read_buf, unsigned int ** dimension_len, int & num_of_dimensions, const char * variable_name)
{
    int i, dimension_all;
    nc_type variable_type;

    // read variable dimension info and variable type info
    ncstat = nc_inq_varid (fileID, variable_name, & variableID);
    if (ncstat != NC_NOERR)
    {
        fprintf(stderr, "%s\n", nc_strerror(ncstat));
        fprintf(stderr, "%s\n", "@ NETCDF_IO_Read_Var()");
        return false;
    }

    ncstat = nc_inq_vartype (fileID, variableID, & variable_type);
    if (ncstat != NC_NOERR)
    {
        fprintf(stderr, "%s\n", nc_strerror(ncstat));
        fprintf(stderr, "%s\n", "@ NETCDF_IO_Read_Var()");
        return false;
    }

    ncstat = nc_inq_varndims (fileID, variableID, & num_of_dimensions);
    if (ncstat != NC_NOERR)
    {
        fprintf(stderr, "%s\n", nc_strerror(ncstat));
        fprintf(stderr, "%s\n", "@ NETCDF_IO_Read_Var()");
        return false;
    }
    
    ncstat = nc_inq_vardimid (fileID, variableID, dimensionID_array);
    if (ncstat != NC_NOERR)
    {
        fprintf(stderr, "%s\n", nc_strerror(ncstat));
        fprintf(stderr, "%s\n", "@ NETCDF_IO_Read_Var()");
        return false;
    }

    dimension_all = 1;
    *dimension_len = (unsigned int *) malloc (sizeof(unsigned int) * num_of_dimensions);
    assert(dimension_len);
    for (i = 0; i < num_of_dimensions; i++)
    {
        assert(NETCDF_IO_Read_Dim((*dimension_len)[i], dimensionID_array[i]));
        dimension_all *= (*dimension_len)[i];
    }

    // allocate memory space for read_buf; read data
    if (variable_type == NC_DOUBLE)
    {
        strcpy (data_type, "double");
        * read_buf = malloc (sizeof(double) * dimension_all);
        ncstat = nc_get_var_double (fileID, variableID, (double *)(*read_buf));
    }
    else if (variable_type == NC_INT)
    {
        strcpy (data_type, "int");
        * read_buf = malloc (sizeof(int   ) * dimension_all);
        ncstat = nc_get_var_int (fileID, variableID, (int *)(*read_buf));
    }
    else if (variable_type == NC_FLOAT)
    {
        strcpy (data_type, "float");
        * read_buf = malloc (sizeof(float ) * dimension_all);
        ncstat = nc_get_var_float (fileID, variableID, (float *)(*read_buf));
    }
    else if (variable_type == NC_SHORT)
    {
        strcpy (data_type, "short");
        * read_buf = malloc (sizeof(short ) * dimension_all);
        ncstat = nc_get_var_short (fileID, variableID, (short *)(*read_buf));
    }
    else if (variable_type == NC_CHAR)
    {
        strcpy (data_type, "char");
        * read_buf = malloc (sizeof(char  ) * dimension_all);
        ncstat = nc_get_var_schar (fileID, variableID, (signed char *)(*read_buf));
    }
    else                   // NC_BYTE
    {
        strcpy (data_type, "unsigned char");
        * read_buf = malloc (sizeof(unsigned char  ) * dimension_all); 
        ncstat = nc_get_var_uchar (fileID, variableID, (unsigned char *)(*read_buf));
    }
    if (ncstat != NC_NOERR)
    {
        fprintf(stderr, "%s\n", nc_strerror(ncstat));
        fprintf(stderr, "%s\n", "@ NETCDF_IO_Read_Var()");
        return false;
    }

    return true;
}

bool NETCDF_IO::NETCDF_IO_Free (void * read_buf, void * dimension_len)
{
    free(read_buf);
    free(dimension_len);
    return true;
}
