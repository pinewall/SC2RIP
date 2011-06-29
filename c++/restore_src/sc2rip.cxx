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
    printf("map_type: %d\n", map_type);
    printf("num_wts: %d\n", num_wts);

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

    // call appropriate interpolation setup routine based on type of remapping requested
    switch(map_type)
    {
        case MAP_TYPE_CONSERV:
            remap_conserv();
            break;
        case MAP_TYPE_BILINEAR:
            remap_bilinear();
            break;
        case MAP_TYPE_BICUBIC:
            remap_bicubic();
            break;
        case MAP_TYPE_DISTWGT:
            remap_distwgt();
            break;
        default:
            cerr << "Invalid Map Type" << endl;
            return 1;

    }
#if _DEBUG_AREA_FRAC_
    //printf("num_links_map: %d\n", num_links_map);
    printf("grid1 area\tgrid1 frac\n");
    for (int i = 0; i < grid1_size; i++)
        printf("%6d\t%3.6f\t%3.6f\n", i, grid1_area[i], grid1_frac[i]);
    printf("grid2 area\tgrid2 frac\n");
    for (int i = 0; i < grid2_size; i++)
        printf("%6d\t%3.6f\t%3.6f\n", i, grid2_area[i], grid2_frac[i]);
#endif

    // reduce size of remapping arrays and then write reammping info to a file
    if (num_links_map != max_links_map)
        resize_remap_vars(num_links_map - max_links_map);
    
    printf("compare result between me and SCRIP before write_remap\n");
    FILE * wts_out;     // file for output weights
    wts_out = fopen("weights.out", "w");
    if (wts_out == (FILE *)0)
        printf("cannot open weights.out file\n");
    fwrite(wts_map, sizeof(double), num_links_map * num_wts, wts_out);
    fclose(wts_out);

    //printf("ready to write weights to NETCDF\n");
    write_remap(map_name, interp_file, output_opt);


#if _TEST_CONSERV_REMAP_
    /* prepare test data from analytical function */
    double * array1 = new double [grid1_size];      // exact data as input
    double * array2 = new double [grid2_size];      // exact data for comparing with remapping results
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
        array1[cell] = 2 + cos(lat);
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
        array2[cell] = 2 + cos(lat);
    #endif
    }

    /* Calculate remapping result */
    double * array2_1st = new double [grid2_size];           // allocate memory for remapping results
    double * array2_2nd = new double [grid2_size];
    double * array2_all = new double [grid2_size];
    
    memset(array2_1st, 0, grid2_size * sizeof(double));      // zero initalization
    memset(array2_2nd, 0, grid2_size * sizeof(double));
    memset(array2_all, 0, grid2_size * sizeof(double));


    /* calculate first-order conservative remapping result */
    num_wts = 3;
    for (int link = 0; link < num_links_map; link ++)
    {
        array2_1st[grid2_add_map[link]] += wts_map[link * num_wts] * array1[grid1_add_map[link]];
    }

    #if _TEST_CONSERV_REMAP_2ND_ORDER_
    /* calculate first-order conservative remapping result */
    double * src_grad_lat  = new double [grid1_size];       // allocate memory for gradient
    double * src_grad_lon  = new double [grid1_size];
    memset(src_grad_lat, 0, sizeof(double) * grid1_size);   // zero initialization
    memset(src_grad_lon, 0, sizeof(double) * grid1_size);
    
    Neighborhood * neighbor  = new Neighborhood("source", -1);

    //neighbor->calculate_gradient_lat (src_grad_lat, src_array, grid1_size);
    //neighbor->calculate_gradient_lon (src_grad_lon, src_array, grid1_size);
    neighbor->calculate_gradient_latlon (src_grad_lat, src_grad_lon, array1, grid1_size);

#if 0
    // test gradients
    for (int cell = 0; cell < grid1_size; cell ++)
    {
        printf("grad_lat(%3.6f) grad_lon(%3.6f)\n", src_grad_lat[cell], src_grad_lon[cell]);
    }
#endif

    for (int link = 0; link < num_links_map; link ++)
    {
        array2_2nd[grid2_add_map[link]] += wts_map[link * num_wts + 1] * src_grad_lat[grid1_add_map[link]];
        array2_2nd[grid2_add_map[link]] += wts_map[link * num_wts + 2] * src_grad_lon[grid1_add_map[link]];
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
    printf("Exact        First        Second\n");
    for (int cell = 0; cell < grid2_size; cell ++)
    {
        printf("%3.6f    %3.6f    %3.6f\n", array2[cell], array2_1st[cell], array2_all[cell]);
    }
#endif
    delete array1;
    delete array2;
    delete array2_1st;
    delete array2_2nd;
    delete array2_all;
    //finalize_intersection();
    finalize_remap_conserv();
    finalize_remap_vars();
    finalize_grids();
    return 0;
}
