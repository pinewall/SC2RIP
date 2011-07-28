#include "kinds.h"
#include "constants.h"
#include "iounits.h"
#include "timers.h"
#include "grids.h"
#include "remap_vars.h"
#include "namelist.h"
#include "utils.h"
#include "remap_write.h"
#include "remap_algorithms.h"

#include "neighborhood.h"
#include "gradient.h"
#include "io.h"

#include "debug.h"

#include <unistd.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>

#define USE_NETCDF_IO       1

bool prepare_test_data (const double * source_center_lat_coords, const double * source_center_lon_coords, int source_grid_size, 
                    const double * destination_center_lat_coords, const double * destination_center_lon_coords, int destination_grid_size,
                    double ** source_array, double ** destination_array, 
                    double ** destination_first_order_remap, double ** destination_second_order_remap,
                    double ** first_order_remap_error, double ** second_order_remap_error,
                    char *  test_case);

bool calculate_remap (double * destination_remap, double destination_grid_size, 
                    const double * source_array, double source_grid_size, char * weight_ncfile);
                    
bool calculate_max_ave_error (const double * destination_array, double destination_grid_size, 
                    const double * destination_first_order_remap, const double * destination_second_order_remap,
                    double * first_order_remap_remap, double * second_order_remap_error);

bool destroy_test_data (double * source_array, double * destination_array,
                    double * destination_first_order_remap, double * destination_second_order_remap,
                    double * first_order_remap_error, double * second_order_remap_error);

int main(int argc, char ** argv)
{
    /* cmdline parsing */
    if (argc < 2)
    {
        printf("Usage: ./sc2rip_test test_case [test_order]\n");
        printf("\t[1] only first-order test\n");
        printf("\t[2] only second-order test\n");
        printf("\t[3] both first and second order tests\n");
        return -1;
    }
    char test_case  [64];
    char test_order [64];
    strcpy(test_case, argv[1]);
    if (argc == 3)
    {
        strcpy(test_order, argv[2]);
    }

    // initialize timers
    for (int i = 0; i < timers->get_max_timers(); i++)
    {
        timers->clear(i);
    }

    // information of remap grids
    grid_init(grid1_file, grid2_file);
    printf("Computing remappings between <<%s>> and <<%s>>\n", grid1_name, grid2_name);

#if CHECK_MASK
    for (int i = 0; i < grid1_size; i++)
    {
        if (grid1_mask[i])
            printf("%d ", 1);
        else
            printf("%d ", 0);
        if ((i+1) % 20 == 0)
            putchar(10);
    }

    for (int i = 0; i < grid2_size; i++)
    {
        if (grid2_mask[i])
            printf("%d ", 1);
        else
            printf("%d ", 0);
        if ((i+1) % 20 == 0)
            putchar(10);
    }
#endif
    /* prepare test data */
    double * source_array, * destination_array;
    double * destination_first_order_remap, * destination_second_order_remap;
    double * first_order_remap_error, * second_order_remap_error;

    printf("grid1_size: %d\n", grid1_size);
    printf("grid2_size: %d\n", grid2_size);

    assert(grid1_size != 0);
    assert(grid2_size != 0);

    prepare_test_data (grid1_center_lat, grid1_center_lon, grid1_size,
                       grid2_center_lat, grid2_center_lon, grid2_size,
                       &source_array, &destination_array,
                       &destination_first_order_remap, &destination_second_order_remap,
                       &first_order_remap_error, &second_order_remap_error,
                       test_case);


    /* calculate conservative remapping result */
    calculate_remap (destination_first_order_remap,     grid2_size, source_array, grid1_size, "SSQoutput1.nc");
    calculate_remap (destination_second_order_remap,    grid2_size, source_array, grid1_size, "SSQoutput2.nc");

    /* maximal and average error */
    calculate_max_ave_error(destination_array, grid2_size, 
            destination_first_order_remap, destination_second_order_remap,
            first_order_remap_error, second_order_remap_error);

    /* may be added netCDF IO */

    /* deallocate memory */
    destroy_test_data(source_array, destination_array, 
                destination_first_order_remap, destination_second_order_remap,
                first_order_remap_error, second_order_remap_error);
    //finalize_intersection();
    finalize_remap_conserv();
    finalize_remap_vars();
    finalize_grids();
    return 0;
}


bool prepare_test_data(const double * source_center_lat_coords, const double * source_center_lon_coords, int source_grid_size, 
                    const double * destination_center_lat_coords, const double * destination_center_lon_coords, int destination_grid_size,
                    double ** source_array, double ** destination_array, 
                    double ** destination_first_order_remap, double ** destination_second_order_remap,
                    double ** first_order_remap_error, double ** second_order_remap_error,
                    char *  test_case)
{
    /* prepare test data from analytical function */
    * source_array                    = (double *) malloc (source_grid_size       * sizeof(double));
    * destination_array               = (double *) malloc (destination_grid_size  * sizeof(double));
    * destination_first_order_remap   = (double *) malloc (destination_grid_size  * sizeof(double));
    * destination_second_order_remap  = (double *) malloc (destination_grid_size  * sizeof(double));
    * first_order_remap_error         = (double *) malloc (destination_grid_size  * sizeof(double));
    * second_order_remap_error        = (double *) malloc (destination_grid_size  * sizeof(double));

    // may be added by netCDF IO
    const bool * source_mask         = grid1_mask;
    const bool * destination_mask    = grid2_mask;
    
    memset(* source_array,                    0, source_grid_size         * sizeof(double));      // zero initalization
    memset(* destination_array,               0, destination_grid_size    * sizeof(double));
    memset(* destination_first_order_remap,   0, destination_grid_size    * sizeof(double));
    memset(* destination_second_order_remap,  0, destination_grid_size    * sizeof(double));
    memset(* first_order_remap_error,         0, destination_grid_size    * sizeof(double));
    memset(* second_order_remap_error,        0, destination_grid_size    * sizeof(double));

    // test case
    double   lat, lon;
    if (strcmp(test_case, "1") == 0)
    {
        double length = 0.1 * 2 * 3.14159265359;
        double grid_tmp;
        for (int cell = 0; cell < source_grid_size; cell ++)
        {
            lat = source_center_lat_coords[cell];
            lon = source_center_lon_coords[cell];
            (*source_array)[cell] = cos(lat) * cos(lon);
            grid_tmp = acos(-(*source_array)[cell]) / length;
            if (grid_tmp <= ONE)
                (*source_array)[cell] = TWO + cos(PI * grid_tmp);
            else
                (*source_array)[cell] = ONE;

            if (!source_mask[cell])
                (*source_array)[cell] = ZERO;
        }
        for (int cell = 0; cell < destination_grid_size; cell ++)
        {
            lat = destination_center_lat_coords[cell];
            lon = destination_center_lon_coords[cell];
            (*destination_array)[cell] = cos(lat) * cos(lon);
            grid_tmp = acos(-(*destination_array)[cell]) / length;
            if (grid_tmp <= ONE)
                (*destination_array)[cell] = TWO + cos(PI * grid_tmp);
            else
                (*destination_array)[cell] = ONE;

            if (!destination_mask[cell])
                (*destination_array)[cell] = ZERO;
        }
    }
    else if (strcmp(test_case, "2") == 0)
    {
        for (int cell = 0; cell < source_grid_size; cell ++)
        {
            lat = source_center_lat_coords[cell];
            lon = source_center_lon_coords[cell];
            if (source_mask[cell])
                (*source_array)[cell] = TWO + cos(lat) * cos(lat) * cos(TWO * lon);
            else
                (*source_array)[cell] = ZERO;
        }
        for (int cell = 0; cell < destination_grid_size; cell ++)
        {
            lat = destination_center_lat_coords[cell];
            lon = destination_center_lon_coords[cell];
            if (destination_mask[cell])
                (*destination_array)[cell] = TWO + cos(lat) * cos(lat) * cos(TWO * lon);
            else
                (*destination_array)[cell] = ZERO;
        }
    }
    else if (strcmp(test_case, "3") == 0)
    {
        for (int cell = 0; cell < source_grid_size; cell ++)
        {
            lat = source_center_lat_coords[cell];
            lon = source_center_lon_coords[cell];
            if (source_mask[cell])
                (*source_array)[cell] = TWO + pow(cos(lat), 16) * cos(16 * lon);
            else
                (*source_array)[cell] = ZERO;
        }
        for (int cell = 0; cell < destination_grid_size; cell ++)
        {
            lat = destination_center_lat_coords[cell];
            lon = destination_center_lon_coords[cell];
            if (destination_mask[cell])
                (*destination_array)[cell] = TWO + pow(cos(lat), 16) * cos(16 * lon);
            else
                (*destination_array)[cell] = ZERO;
        }
    }
    else
    {
        printf("Test Case:\n\t[1] Cosine Hill Function\n\t[2] Y(2,2) Spherical Function\n\t[3] Y(16, 32) Spherical Function\n");
        return false;
    }
    return true;
}

bool destroy_test_data(double * source_array, double * destination_array,
                    double * destination_first_order_remap, double * destination_second_order_remap,
                    double * first_order_remap_error, double * second_order_remap_error)
{
    free(source_array);
    free(destination_array);
    free(destination_first_order_remap);
    free(destination_second_order_remap);
    free(first_order_remap_error);
    free(second_order_remap_error);
    if (source_array == NULL && destination_array == NULL &&
            NULL == destination_first_order_remap && NULL == destination_second_order_remap &&
            NULL == first_order_remap_error && NULL == second_order_remap_error)
        return true;
    else
        return false;
}

bool calculate_remap (double * destination_remap, double destination_grid_size, 
                    const double * source_array, double source_grid_size, char * weight_ncfile)
{
    unsigned int        num_links           = 0;
    unsigned int        num_wgts            = 0;
    int                 src_index_rank      = 0;
    int                 dst_index_rank      = 0;
    int                 remap_matrix_rank   = 0;
    unsigned int    *   src_index_dimension_lens;
    unsigned int    *   dst_index_dimension_lens;
    unsigned int    *   remap_matrix_dimension_lens;
    unsigned int    *   src_index;
    unsigned int    *   dst_index;
    double          *   weights; 
#if USE_NETCDF_IO
    NETCDF_IO * netcdf = new NETCDF_IO ("netcdf", weight_ncfile, FILE_STATUS_OLD);
    assert(netcdf != NULL);

    assert(netcdf->NETCDF_IO_Open (NC_NOCLOBBER));

    // read dimensions
    assert(netcdf->NETCDF_IO_Read_Dim (num_links, "num_links"));
    assert(netcdf->NETCDF_IO_Read_Dim (num_wgts, "num_wgts"));

    // read variables
    assert(netcdf->NETCDF_IO_Read_Var ((void **)(&src_index), &src_index_dimension_lens, src_index_rank, "src_address"));
    assert(netcdf->NETCDF_IO_Read_Var ((void **)(&dst_index), &dst_index_dimension_lens, dst_index_rank, "dst_address"));
    assert(netcdf->NETCDF_IO_Read_Var ((void **)(&weights), &remap_matrix_dimension_lens, remap_matrix_rank, "remap_matrix"));
#else
    int ncid;
    nc_open (weight_ncfile, NC_NOCLOBBER, &ncid);
    int num_links_dim_id, num_wgts_dim_id;
    int src_index_id, dst_index_id, weights_id;
    nc_inq_dimid (ncid, "num_links", &num_links_dim_id);
    nc_inq_dimid (ncid, "num_wgts", &num_wgts_dim_id);
    nc_inq_varid (ncid, "src_address", &src_index_id);
    nc_inq_varid (ncid, "dst_address", &dst_index_id);
    nc_inq_varid (ncid, "remap_matrix", &weights_id);

    nc_inq_dimlen (ncid, num_links_dim_id, &num_links);
    nc_inq_dimlen (ncid, num_wgts_dim_id, &num_wgts);
    printf("num_links = %d\n", num_links);
    printf("num_wgts = %d\n", num_wgts);

    src_index = (int *) malloc (sizeof(int) * num_links);
    dst_index = (int *) malloc (sizeof(int) * num_links);
    weights = (double *) malloc (sizeof(double) * num_links * num_wgts);

    nc_get_var_int (ncid, src_index_id, src_index);
    nc_get_var_int (ncid, dst_index_id, dst_index);
    nc_get_var_double (ncid, weights_id, weights);

    nc_close(ncid);
#endif
    // matrix vector multiply in sparse pattern
    for (int link = 0; link < num_links; link ++)
    {
        //printf("(%d\t%d\t%f\n", src_index[link], dst_index[link], weights[link * num_wgts]);
        assert(src_index[link] < source_grid_size);
        assert(dst_index[link] < destination_grid_size);
        destination_remap[dst_index[link]] += source_array[src_index[link]] * weights[link * num_wgts];
    }

    // deallocate memory
#if USE_NETCDF_IO
    assert(netcdf->NETCDF_IO_Free (src_index, src_index_dimension_lens));
    assert(netcdf->NETCDF_IO_Free (dst_index, dst_index_dimension_lens));
    assert(netcdf->NETCDF_IO_Free (weights, remap_matrix_dimension_lens));
#else
    free(src_index);
    free(dst_index);
    free(weights);
#endif
    return true;
}

bool calculate_max_ave_error (const double * destination_array, double destination_grid_size, 
                    const double * destination_first_order_remap, const double * destination_second_order_remap,
                    double * first_order_remap_error, double * second_order_remap_error)
{
    double  max_error1 = -999;
    double  min_error1 =  999;
    int     max_index1 =    0;
    int     min_index1 =    0;  
    double  max_error2 = -999;
    double  min_error2 =  999;
    int     max_index2 =    0;
    int     min_index2 =    0;
    double  ave_error1 =  0.0;
    double  ave_error2 =  0.0;
    double  tmp_error1 =  0.0;
    double  tmp_error2 =  0.0;

    printf("Exact        First        Second\n");
    for (int cell = 0; cell < destination_grid_size; cell ++)
    {
        first_order_remap_error[cell]   = destination_first_order_remap[cell]   - destination_array[cell];
        second_order_remap_error[cell]  = destination_second_order_remap[cell]  - destination_array[cell];
        //printf("%3.6f    %3.6f    %3.6f\n", destination_array[cell], destination_first_order_remap[cell], destination_second_order_remap[cell]);
        if (nonzero(destination_array[cell]))
        {
            first_order_remap_error[cell]  /= destination_array[cell];
            tmp_error1 = first_order_remap_error[cell];
            second_order_remap_error[cell] /= destination_array[cell];
            tmp_error2 = second_order_remap_error[cell];
            if (tmp_error1 > max_error1)
            {
                max_error1 = tmp_error1;
                max_index1 = cell;
            }
            
            if (tmp_error1 < min_error1)
            {
                min_error1 = tmp_error1;
                min_index1 = cell;
            }

            if (tmp_error2 > max_error2)
            {
                max_error2 = tmp_error2;
                max_index2 = cell;
            }

            if (tmp_error2 < min_error2)
            {
                min_error2 = tmp_error2;
                min_index2 = cell;
            }
        }
        ave_error1 += first_order_remap_error[cell];
        ave_error2 += second_order_remap_error[cell];
    }
    ave_error1 /= destination_grid_size;
    ave_error2 /= destination_grid_size;
    printf("Max First-order Error %3.16f @ %d\n", max_error1, max_index1);
    printf("Min First-order Error %3.16f @ %d\n", min_error1, min_index1);
    printf("Ave First-order Error %3.16f\n", ave_error1);
    printf("Max Second-order Error %3.16f @ %d\n", max_error2, max_index2);
    printf("Min Second-order Error %3.16f @ %d\n", min_error2, min_index2);
    printf("Ave Second-order Error %3.16f\n", ave_error2);
    return true;
}
