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

#include "debug.h"

#include <unistd.h>
#include <stdio.h>
#include <iostream>

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
    //write_remap(map_name, interp_file, output_opt);
    //finalize_remap_vars();
    //finalize_intersection();

    return 0;
}
