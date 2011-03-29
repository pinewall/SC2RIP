#include "timers.h"
#include "iounits.h"
#include "grids.h"
#include <unistd.h>
#include <stdio.h>
#include <iostream>
using namespace std;
int main()
{
    Timers *timers = new Timers();
    timers->start(0);
    sleep(1);
    timers->stop(0);
    cout << timers->get(0) << endl;
    timers->print(0);

    int unit = get_unit();
    release_unit(unit);

    //netcdf_error_handler(NC_NOERR+1);

    // test grid util
    grid_init("T42.nc", "POP43.nc");
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
            printf("%1.5f  ", grid2_center_lat[index + j]);
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
            printf("%1.5f  ", grid2_corner_lat[index + j]);
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
    return 0;
}
