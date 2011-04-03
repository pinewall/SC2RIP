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

#include <unistd.h>
#include <stdio.h>
#include <iostream>

#define _DEBUG_ 1
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
    
    // initialize grid information for both grids
    grid_init(grid1_file, grid2_file);
    cout << "Computing remappings between <<" << grid1_name << ">> and <<" << grid2_name << ">>" << endl;
    
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

    // reduce size of remapping arrays and then write reammping info to a file
    if (num_links_map != max_links_map)
        resize_remap_vars(num_links_map - max_links_map);
    write_remap(map_name, interp_file, output_opt);

#if _DEBUG_
    // test grid util
    cout << "grid1_size = " << grid1_size << endl;
    cout << "grid1_rank = " << grid1_rank << endl;
    cout << "grid1_mask[0] = " << grid1_mask[0] << endl;
    //cout << "grid1_center_lat[0] = " << grid1_center_lat[0] << endl;
    //cout << "grid1_center_lon[0] = " << grid1_center_lon[0] << endl;
    cout << "grid1_corners_max = " << grid1_corners_max << endl;
    cout << "GRID1 corner lats" << endl;
    int index = 0;
    for (int i = 0; i < grid1_size; i++)
    {
        for (int j = 0; j < grid1_corners_max; j++)
        {
            //cout << grid1_center_lat[index + j] << "\t";
            printf("%1.5f  ", grid1_corner_lat[index + j]);
        }
        //cout << endl;
        printf("\n");
        index += grid1_corners_max;
    }
    cout << "GRID1 corner lons" << endl;
    index = 0;
    for (int i = 0; i < grid1_size; i++)
    {
        for (int j = 0; j < grid1_corners_max; j++)
        {
            //cout << grid1_corner_lat[index + j] << "\t";
            printf("%1.5f  ", grid1_corner_lon[index + j]);
        }
        //cout << endl;
        printf("\n");
        index += grid1_corners_max;
    }

    // test bounding box
    cout << "GRID1 bounding box" << endl;
    index = 0;
    for (int i = 0; i < grid1_size; i++)
    {
        //cout << grid1_bound_box[index++] << "\t";
        //cout << grid1_bound_box[index++] << "\t";
        //cout << grid1_bound_box[index++] << "\t";
        //cout << grid1_bound_box[index++] << "\t";
        //cout << grid1_bound_box[index++] << "\t";
        //cout << grid1_bound_box[index++] << endl;
        printf("%1.5f  ", grid1_bound_box[index++]);
        printf("%1.5f  ", grid1_bound_box[index++]);
        printf("%1.5f  ", grid1_bound_box[index++]);
        printf("%1.5f  ", grid1_bound_box[index++]);
        printf("%1.5f  ", grid1_bound_box[index++]);
        printf("%1.5f\n", grid1_bound_box[index++]);
    }
    cout << "grid2_size = " << grid2_size << endl;
    cout << "grid2_rank = " << grid2_rank << endl;
    cout << "grid2_mask[0] = " << grid2_mask[0] << endl;
    //cout << "grid2_center_lat[0] = " << grid2_center_lat[0] << endl;
    //cout << "grid2_center_lon[0] = " << grid2_center_lon[0] << endl;
    cout << "grid2_corners_max = " << grid2_corners_max << endl;
    cout << "GRID2 corner lats" << endl;
    index = 0;
    for (int i = 0; i < grid2_size; i++)
    {
        for (int j = 0; j < grid2_corners_max; j++)
        {
            //cout << grid2_center_lat[index + j] << "\t";
            printf("%1.5f  ", grid2_corner_lat[index + j]);
        }
        //cout << endl;
        printf("\n");
        index += grid2_corners_max;
    }
    cout << "GRID2 corner lons" << endl;
    index = 0;
    for (int i = 0; i < grid2_size; i++)
    {
        for (int j = 0; j < grid2_corners_max; j++)
        {
            //cout << grid2_corner_lat[index + j] << "\t";
            printf("%1.5f  ", grid2_corner_lon[index + j]);
        }
        //cout << endl;
        printf("\n");
        index += grid2_corners_max;
    }

    // test bounding box
    cout << "GRID2 bounding box" << endl;
    index = 0;
    for (int i = 0; i < grid2_size; i++)
    {
        //cout << grid2_bound_box[index++] << "\t";
        //cout << grid2_bound_box[index++] << "\t";
        //cout << grid2_bound_box[index++] << "\t";
        //cout << grid2_bound_box[index++] << "\t";
        //cout << grid2_bound_box[index++] << "\t";
        //cout << grid2_bound_box[index++] << endl;
        printf("%1.5f  ", grid2_bound_box[index++]);
        printf("%1.5f  ", grid2_bound_box[index++]);
        printf("%1.5f  ", grid2_bound_box[index++]);
        printf("%1.5f  ", grid2_bound_box[index++]);
        printf("%1.5f  ", grid2_bound_box[index++]);
        printf("%1.5f\n", grid2_bound_box[index++]);
    }
#endif
    return 0;
}
