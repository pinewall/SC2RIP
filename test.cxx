#include "timers.h"
#include "iounits.h"
#include "grids.h"
#include <unistd.h>
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
    for (int i = 0; i < 10; i++)
    {
        int index = 0;
        for (int j = 0; j < grid1_corners_max; j++)
        {
            cout << grid1_center_lat[index + j] << "\t";
        }
        cout << endl;
        index += grid1_corners_max;
    }
    for (int i = 0; i < 10; i++)
    {
        int index = 0;
        for (int j = 0; j < grid1_corners_max; j++)
        {
            cout << grid1_corner_lat[index + j] << "\t";
        }
        cout << endl;
        index += grid1_corners_max;
    }

    // test bounding box
    cout << "bounding box of grid1" << endl;
    int index = 0;
    for (int i = 0; i < grid1_size; i++)
    {
        cout << grid1_bound_box[index++] << "\t";
        cout << grid1_bound_box[index++] << "\t";
        cout << grid1_bound_box[index++] << "\t";
        cout << grid1_bound_box[index++] << "\t";
        cout << grid1_bound_box[index++] << "\t";
        cout << grid1_bound_box[index++] << endl;
    }
    return 0;
}
