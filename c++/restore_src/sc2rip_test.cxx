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

#include "debug.h"

#include <unistd.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <iostream>

#define _TEST_CONSERV_REMAP_            1
#define _TEST_CONSERV_REMAP_2ND_ORDER_  1
#define _TEST_CASE_                     3
#define _TEST_CASE_1_                   ((_TEST_CASE_ - 2)*(_TEST_CASE_ - 3))
#define _TEST_CASE_2_                   ((_TEST_CASE_ - 3)*(_TEST_CASE_ - 1))
#define _TEST_CASE_3_                   ((_TEST_CASE_ - 1)*(_TEST_CASE_ - 2))
using namespace std;
int main()
{

    // initialize timers
    for (int i = 0; i < timers->get_max_timers(); i++)
    {
        timers->clear(i);
    }

    // read input namelist; here we use namelist.cxx instead
    int unit = get_unit();
    release_unit(unit);

    // remap method settings; map_type declared in remap_vars.h
    if (streqls(map_method, "conservative"))
    {
        map_type = MAP_TYPE_CONSERV;
        luse_grid_centers = false;
    }
    else if (streqls(map_method, "bilinear"))
    {
        map_type = MAP_TYPE_BILINEAR;
        luse_grid_centers = true;
    }
    else if (streqls(map_method, "bicubic"))
    {
        map_type = MAP_TYPE_BICUBIC;
        luse_grid_centers = true;
    }
    else if (streqls(map_method, "distwgt"))
    {
        map_type = MAP_TYPE_BICUBIC;
        luse_grid_centers = true;
    }
    else
    {
        cerr << "Unknown mapping method" << endl;
        return 1;
    }

    // normalize option settings; norm_opt declared in remap_vars.h
    if (streqls(normalize_opt, "none"))
        norm_opt = NORM_OPT_NONE;
    else if (streqls(normalize_opt, "fracarea"))
        norm_opt = NORM_OPT_FRACAREA;
    else if (streqls(normalize_opt, "destarea"))
        norm_opt = NORM_OPT_DESTAREA;
    else
    {
        cerr << "Unknown normalization option" << endl;
        return 1;
    }
    printf("norm_opt: %d\n", norm_opt);
    
    // initialize grid information for both grids
    grid_init(grid1_file, grid2_file);
    cout << "Computing remappings between <<" << grid1_name << ">> and <<" << grid2_name << ">>" << endl;
    
#if _DEBUG_GRID_
    // test grid util
    grid_debug();
#endif
    // initialize some remapping variables
    init_remap_vars();
    printf("map_type: %d\n", map_type);
    printf("num_wts: %d\n", num_wts);
    int     * scrip_source_cell_index = new int [63112];
    int     * scrip_destination_cell_index = new int [63112];
    double  * scrip_remap_matrix = new double [63112 * 3];

    FILE    * scrip_wts;
    scrip_wts = fopen("src_address", "rb");
    fread(scrip_source_cell_index, sizeof(int), 63112, scrip_wts);
    fclose(scrip_wts);

    scrip_wts = fopen("dst_address", "rb");
    fread(scrip_destination_cell_index, sizeof(int), 63112, scrip_wts);
    fclose(scrip_wts);

    scrip_wts = fopen("remap_matrix", "rb");
    fread(scrip_remap_matrix, sizeof(double), 63112 * 3, scrip_wts);
    fclose(scrip_wts);

    for (int i = 0; i < 63112; i ++)
    {
        scrip_source_cell_index[i] ++;
        scrip_destination_cell_index[i] ++;
    }

    /* prepare test data from analytical function */
    double * array1 = new double [grid1_size];      // exact data as input
    double * array2 = new double [grid2_size];      // exact data for comparing with remapping results
    printf("grid1_size: %d\n", grid1_size);
    printf("grid2_size: %d\n", grid2_size);
    double   lat, lon;
    for (int cell = 0; cell < grid1_size; cell ++)
    {
        lat = grid1_center_lat[cell];
        lon = grid1_center_lon[cell];
    #if _TEST_CASE_1_
        array1[cell] = 1.0;
    #endif
    #if _TEST_CASE_2_
        array1[cell] = cos(lat) * cos(lon);
        //array1[cell] = cos(lat) * cos(lon);   // for conservative remapping, point value is not a good method, we should calculate area-average value
    #endif
    #if _TEST_CASE_3_
        array1[cell] = 2 + cos(lon);
    #endif
    }

    for (int cell = 0; cell < grid2_size; cell ++)
    {
        lat = grid2_center_lat[cell];
        lon = grid2_center_lon[cell];
    #if _TEST_CASE_1_
        array2[cell] = 1.0;
    #endif
    #if _TEST_CASE_2_
        array2[cell] = cos(lat) * cos(lon);
    #endif
    #if _TEST_CASE_3_
        array2[cell] = 2 + cos(lon);
    #endif
    }

    /* Calculate remapping result */
    double * array2_1st = new double [grid2_size];           // allocate memory for remapping results
    double * array2_2nd = new double [grid2_size];
    double * array2_all = new double [grid2_size];
    double *    error1  = new double [grid2_size];
    double *    error2  = new double [grid2_size];
    
    memset(array2_1st, 0, grid2_size * sizeof(double));      // zero initalization
    memset(array2_2nd, 0, grid2_size * sizeof(double));
    memset(array2_all, 0, grid2_size * sizeof(double));
    memset(error1,     0, grid2_size * sizeof(double));
    memset(error2,     0, grid2_size * sizeof(double));


    /* calculate first-order conservative remapping result */
    for (int link = 0; link < 63112; link ++)
    {
        array2_1st[scrip_destination_cell_index[link]] += scrip_remap_matrix[link * 3] * array1[scrip_source_cell_index[link]];
    }
#if 0
    for (int link = 0; link < num_links_map; link ++)
    {
        array2_1st[grid2_add_map[link]] += wts_map[link * num_wts] * array1[grid1_add_map[link]];
    }
#endif

    #if _TEST_CONSERV_REMAP_2ND_ORDER_
#if 0
    double * src_grad_lat  = new double [grid1_size];       // allocate memory for gradient
    double * src_grad_lon  = new double [grid1_size];
    memset(src_grad_lat, 0, sizeof(double) * grid1_size);   // zero initialization
    memset(src_grad_lon, 0, sizeof(double) * grid1_size);
#endif

    Neighborhood    * neighbor  = new Neighborhood("source", grid1_size, -1);
    Gradient        * gradient  = new Gradient("latlon", grid1_size, array1, grid1_center_lat, grid1_center_lon);
    for (int cell = 0; cell < grid1_size; cell ++)
    {
        gradient->calculate_gradient(cell, array1[cell], grid1_center_lat[cell], grid1_center_lon[cell], neighbor->get_index_of_neighbor_cells(cell), neighbor->get_num_of_neighbor_cells(cell));
    }

#if 0
    // test gradients
    for (int cell = 0; cell < grid1_size; cell ++)
    {
        printf("grad_lat(%3.6f) grad_lon(%3.6f)\n", src_grad_lat[cell], src_grad_lon[cell]);
    }
#endif
#if 0
    for (int link = 0; link < num_links_map; link ++)
    {
        array2_2nd[grid2_add_map[link]] += wts_map[link * num_wts + 1] * src_grad_lat[grid1_add_map[link]];
        array2_2nd[grid2_add_map[link]] += wts_map[link * num_wts + 2] * src_grad_lon[grid1_add_map[link]];
    }
#endif
    for (int link = 0; link < 63112; link ++)
    {
        array2_2nd[scrip_destination_cell_index[link]] += scrip_remap_matrix[link * 3 + 1] * gradient->get_gradient(scrip_source_cell_index[link])[0];
        array2_2nd[scrip_destination_cell_index[link]] += scrip_remap_matrix[link * 3 + 2] * gradient->get_gradient(scrip_source_cell_index[link])[1];
    }
    #endif

    /* add remapping results together */
    for (int cell = 0; cell < grid2_size; cell ++)
    {
        array2_all[cell] = array2_1st[cell] + array2_2nd[cell];
    }
#if 0
    /* test result */
    for (int cell = 0; cell < grid2_size; cell ++)
    {
        if (grid2_frac[cell] > 0.001)
        {
            array2_1st[cell] /= grid2_frac[cell];
            array2_all[cell] /= grid2_frac[cell];
        }
        else
        {
            array2_1st[cell] = 0.0;
            array2_all[cell] = 0.0;
        }
    }
#endif

    double max_error1 = -1;
    int    max_index1 = 0;
    double max_error2 = -1;
    int    max_index2 = 0;
    double ave_error1 = 0.0;
    double ave_error2 = 0.0;
    double tmp_error1 = 0.0;
    double tmp_error2 = 0.0;

    printf("Exact        First        Second\n");
    for (int cell = 0; cell < grid2_size; cell ++)
    {
        error1[cell] = ABS(array2_1st[cell] - array2[cell]);
        error2[cell] = ABS(array2_all[cell] - array2[cell]);
        //printf("%3.6f    %3.6f    %3.6f\n", array2[cell], array2_1st[cell], array2_all[cell]);
        if (nonzero(array2[cell]))
        {
            tmp_error1 = error1[cell] / array2[cell];
            tmp_error2 = error2[cell] / array2[cell];
            if (tmp_error1 > max_error1)
            {
                max_error1 = tmp_error1;
                max_index1 = cell;
            }

            if (tmp_error2 > max_error2)
            {
                max_error2 = tmp_error2;
                max_index2 = cell;
            }

            ave_error1 += error1[cell];
            ave_error2 += error2[cell];
        }
    }
    ave_error1 /= grid2_size;
    ave_error2 /= grid2_size;
    printf("Max First-order Error %3.6f @ %d\n", max_error1, max_index1);
    printf("Ave First-order Error %3.6f\n", ave_error1);
    printf("Max Second-order Error %3.6f @ %d\n", max_error2, max_index2);
    printf("Ave Second-order Error %3.6f\n", ave_error2);

    delete array1;
    delete array2;
    delete array2_1st;
    delete array2_2nd;
    delete array2_all;
    delete error1;
    delete error2;
    delete scrip_source_cell_index;
    delete scrip_destination_cell_index;
    delete scrip_remap_matrix;
    //finalize_intersection();
    finalize_remap_conserv();
    finalize_remap_vars();
    finalize_grids();
    return 0;
}
