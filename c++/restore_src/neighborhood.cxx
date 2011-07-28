/*
 *  This Class is used to generate neighborhood info for gradient calculating (currently least-square method only)
 *  Modification Time: 2011-06-22
 *  Authoer: pinewall
 */

#include "neighborhood.h"
#include <string.h>

Neighborhood::Neighborhood ( char * grid_name, int grid_size, int num_of_neighbors )
{
    /* initialization of data */
    num_of_cells = grid_size;
    num_of_neighbor_cells   = new int [num_of_cells];

    // allocate memory
    index_of_neighbor_cells = new int * [num_of_cells];

    // just use level1 neighbors
    if (num_of_neighbors == -1)
    {
        num_of_different_lats = 64;
        num_of_different_lons = 128;
        double num_of_level1_neighbors;
        int cell = 0;
        for (int lat = 0; lat < num_of_different_lats; lat ++)
        {
            if (lat == 0 || lat == num_of_different_lats - 1)
                num_of_level1_neighbors = 5;
            else
                num_of_level1_neighbors = 8;
            for (int lon = 0; lon < num_of_different_lons; lon ++)
            {
                num_of_neighbor_cells[cell] = num_of_level1_neighbors;
                index_of_neighbor_cells[cell] = new int [num_of_neighbor_cells[cell]];
                cell ++;
            }
        }
    }
    else
    {
        num_of_neighbor_cells = new int [num_of_cells];
        for (int cell = 0; cell < num_of_cells; cell ++)
        {
            num_of_neighbor_cells[cell] = num_of_neighbors;
            index_of_neighbor_cells[cell] = new int [num_of_neighbors];
        }
    }

    /* calculate neighborhood */
    if (num_of_neighbors == -1)
    {
        find_latlon_neighborhood();
    }
    else
    {
        find_exact_num_neighborhood();
    }
}

/*
 *  eight points case: 
 *      6   5   4
 *      7   *   3
 *      0   1   2
 *  five points case:
 *      0   *   4       or      3   2   1
 *      1   2   3               4   *   0
 */
void Neighborhood::find_latlon_neighborhood()
{
    int left, right;
    int north, south, curr;
    int cell = 0;
    for (int lat = 0; lat < num_of_different_lats; lat ++)
    {
        north   = (lat + 1) * num_of_different_lons;
        south   = (lat - 1) * num_of_different_lons;
        curr    = lat * num_of_different_lons;
        for (int lon = 0; lon < num_of_different_lons; lon ++)
        {
            left    = (lon + num_of_different_lons - 1) % num_of_different_lons;
            right   = (lon + 1) % num_of_different_lons;
            if (lat == num_of_different_lats - 1)
            {
                index_of_neighbor_cells[cell][0] = curr + left;
                index_of_neighbor_cells[cell][1] = south + left;
                index_of_neighbor_cells[cell][2] = south + lon;
                index_of_neighbor_cells[cell][3] = south + right;
                index_of_neighbor_cells[cell][4] = curr + right;
            }
            else if (lat == 0)
            {
                index_of_neighbor_cells[cell][0] = curr + right;
                index_of_neighbor_cells[cell][1] = north + right;
                index_of_neighbor_cells[cell][2] = north + lon;
                index_of_neighbor_cells[cell][3] = north + left;
                index_of_neighbor_cells[cell][4] = curr + left;
            }
            else
            {
                index_of_neighbor_cells[cell][0] = south + left;
                index_of_neighbor_cells[cell][1] = south + lon;
                index_of_neighbor_cells[cell][2] = south + right;
                index_of_neighbor_cells[cell][3] = curr + right;
                index_of_neighbor_cells[cell][4] = north + right;
                index_of_neighbor_cells[cell][5] = north + lon;
                index_of_neighbor_cells[cell][6] = north + left;
                index_of_neighbor_cells[cell][7] = curr + left;
            }
            cell ++;
        }
    }
}

void Neighborhood::find_exact_num_neighborhood()
{
}

int Neighborhood::get_num_of_cells()
{
    return num_of_cells;
}

int Neighborhood::get_num_of_neighbor_cells(int cell)
{
    return num_of_neighbor_cells[cell];
}

int * Neighborhood::get_index_of_neighbor_cells(int cell)
{
    return index_of_neighbor_cells[cell];
}
