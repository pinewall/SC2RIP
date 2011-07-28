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
#include <assert.h>

#define USE_SECOND_ORDER    1

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

#if BINARY_WEIGHT
    printf("compare result between me and SCRIP before write_remap\n");
    FILE * wts_out;     // file for output weights
    wts_out = fopen("weights.out", "w");
    if (wts_out == (FILE *)0)
        printf("cannot open weights.out file\n");
    fwrite(wts_map, sizeof(double), num_links_map * num_wts, wts_out);
    fclose(wts_out);
#endif

    // write SCRIP conservative remap matrix
    write_remap(map_name, interp_file, output_opt);
    
#if USE_SECOND_ORDER
    /* prepare neighborhood and gradient */
    Neighborhood    *   neighbor            = new Neighborhood("source", grid1_size, -1);
    double          *   delta_theta         = new double [128];    
    double          *   delta_varphi        = new double [128];
    double              sum_theta;
    double              sum_varphi;
    double              sum_theta_varphi;   
    double              sum_theta_square;  
    double              sum_varphi_square;
    double              determine;       
    
    int                 num_grad_links      = 0;
    for (int i = 0; i < neighbor->get_num_of_cells(); i++)          // consider neighbors
        num_grad_links += neighbor->get_num_of_neighbor_cells(i);
    num_grad_links += grid1_size;                                   // consider diaginal line

    int             *   row_index           = new int [num_grad_links];
    int             *   col_index           = new int [num_grad_links];
    double          *   grad_lat            = new double [num_grad_links];
    double          *   grad_lon            = new double [num_grad_links];

    int sparse_matrix_id = 0;
    bool diag_element = true;
    for (int cell = 0; cell < grid1_size; cell ++)
    {
        int     num_neighbor = neighbor->get_num_of_neighbor_cells(cell);
        int *   neighbor_id  = neighbor->get_index_of_neighbor_cells(cell);
        sum_theta = 0;
        sum_varphi = 0;
        for (int neighbor_cell = 0; neighbor_cell < num_neighbor; neighbor_cell ++)
        {
            delta_theta[neighbor_cell] = grid1_center_lat[neighbor_id[neighbor_cell]] - grid1_center_lat[cell];
            delta_varphi[neighbor_cell] = grid1_center_lon[neighbor_id[neighbor_cell]] - grid1_center_lon[cell];
            check_longitude(delta_varphi[neighbor_cell], -PI, PI);
            sum_theta               += delta_theta[neighbor_cell];
            sum_varphi              += delta_varphi[neighbor_cell];
            sum_theta_varphi  += delta_theta[neighbor_cell] * delta_varphi[neighbor_cell];
            sum_theta_square  += delta_theta[neighbor_cell] * delta_theta[neighbor_cell];
            sum_varphi_square += delta_varphi[neighbor_cell] * delta_varphi[neighbor_cell];
        }
        determine = sum_theta_square * sum_varphi_square 
                  - sum_theta_varphi * sum_theta_varphi;
        assert(determine > 1e-30);

        // write value to matrix
        diag_element = true;
        for (int neighbor_cell = 0; neighbor_cell < num_neighbor; neighbor_cell ++)
        {
            if (diag_element && neighbor_cell < num_neighbor - 1 && neighbor_id[neighbor_cell + 1] > cell)
            {
                diag_element = false;
                row_index[sparse_matrix_id] = cell;
                col_index[sparse_matrix_id] = cell;
                grad_lat [sparse_matrix_id] = - sum_varphi_square * sum_theta + sum_theta_varphi * sum_varphi;
                grad_lat [sparse_matrix_id] /= determine;

                grad_lon [sparse_matrix_id] = - sum_theta_square * sum_varphi + sum_theta_varphi * sum_theta;
                grad_lon [sparse_matrix_id] /= determine;
                sparse_matrix_id ++;
            }
            row_index[sparse_matrix_id] = cell;
            col_index[sparse_matrix_id] = neighbor_cell;
            grad_lat [sparse_matrix_id] = sum_varphi_square * delta_theta[neighbor_cell] - sum_theta_varphi * delta_varphi[neighbor_cell];
            grad_lat [sparse_matrix_id] /= determine;

            grad_lon [sparse_matrix_id] = sum_theta_square * delta_varphi[neighbor_cell] - sum_theta_varphi * delta_theta[neighbor_cell];
            grad_lon [sparse_matrix_id] /= determine;
            sparse_matrix_id ++;

            if (diag_element && neighbor_cell == num_neighbor - 1)
            {
                diag_element = false;
                row_index[sparse_matrix_id] = cell;
                col_index[sparse_matrix_id] = cell;
                grad_lat [sparse_matrix_id] = - sum_varphi_square * sum_theta + sum_theta_varphi * sum_varphi;
                grad_lat [sparse_matrix_id] /= determine;

                grad_lon [sparse_matrix_id] = - sum_theta_square * sum_varphi + sum_theta_varphi * sum_theta;
                grad_lon [sparse_matrix_id] /= determine;
                sparse_matrix_id ++;
            }
        }
    }
    assert(sparse_matrix_id == num_grad_links);

#if CHECK_GRADIENT
    printf("grat_lat\n");
    for (int i = 0; i < num_grad_links; i++)
    {
        if (grad_lat[i] > 0)
            printf(" ");
        printf("%5.10f\n", grad_lat[i]);
    }
    printf("grad_lon\n");
    for (int i = 0; i < num_grad_links; i ++)
    {
        if (grad_lon[i] > 0)
            printf(" ");
        printf("%5.10f\n", grad_lon[i]);
    }
#endif
    
    int old_num_links = num_links_map;
    int * src_address = new int [num_links_map];
    int * dst_address = new int [num_links_map];
    double * map_weights  = new double [num_links_map * num_wts];
    memcpy(src_address, grid1_add_map, sizeof(int) * num_links_map);
    memcpy(dst_address, grid2_add_map, sizeof(int) * num_links_map);
    memcpy(map_weights, wts_map, sizeof(double) * num_links_map * num_wts);

    finalize_remap_vars();
    delete [] src_link_add;
    delete [] dst_link_add;

    map_type = 999;                     // to make num_wgts = 1
    first_call_store = true;
    init_remap_vars();
    printf("max_links_map: %d\n", max_links_map);
    printf("num_links_map: %d\n", num_links_map);
    printf("old_num_links: %d\n", old_num_links);
    printf("num_grad_links: %d\n", num_grad_links);
    printf("map_type: %d\n", map_type);
    printf("num_wts: %d\n", num_wts);

    /*
       dst*src by src*src 
       ==> (dst_w, src_w, value_w) -- (row_g, col_g, grad_g) any two triple entries
       ==> if (src_w == row_g) then value_w * grad_g contribues to (dst_w, col_g) in result matrix
    */
    double sparse_tmp;
    double * weight_tmp = new double [1];
    for (int w = 0; w < old_num_links; w ++)
    {
        * weight_tmp = map_weights[w * 3];
        store_link_cnsrv(src_address[w], dst_address[w], weight_tmp, num_wts);
        for (int g = 0; g < num_grad_links; g ++)
        {
            if (src_address[w] == row_index[g])
            {
                sparse_tmp =    map_weights[w * 3 + 1] * grad_lat[g];
                sparse_tmp +=   map_weights[w * 3 + 2] * grad_lon[g];
                store_link_cnsrv(col_index[g], dst_address[w],  &sparse_tmp, num_wts);
            }
        }
    }
#endif
    printf("num_links_map: %d\n", num_links_map);
    printf("max_links_map: %d\n", max_links_map);
    // reduce size of remapping arrays and then write reammping info to a file
    if (num_links_map != max_links_map)
        resize_remap_vars(num_links_map - max_links_map);

    // write our conservative remap matrix
    write_remap(map_name, "SSQoutput2.nc", output_opt);

#if USE_SECOND_ORDER
    // free
    delete [] row_index;
    delete [] col_index;
    delete [] grad_lat;
    delete [] grad_lon;
    delete [] delta_theta;
    delete [] delta_varphi;
#endif
    //finalize_intersection();
    finalize_remap_conserv();
    finalize_remap_vars();
    finalize_grids();
    return 0;
}
