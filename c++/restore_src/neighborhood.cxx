/*
 *  This Class is used to generate neighborhood info for gradient calculating (currently least-square method only)
 *  Modification Time: 2011-06-22
 *  Authoer: pinewall
 */

#include "grids.h"
#include "gradient.h"
#include "remap_conserv.h"
#include "neighborhood.h"
#include <string.h>

Neighborhood::Neighborhood ( char * grid_name, int num_of_neighbors )
{
    /* initialization of data */
    if (strcmp(grid_name,"source") == 0)
    {
        num_of_cells = grid1_size;
    }
    else if (strcmp(grid_name, "destination") == 0)
    {
        num_of_cells = grid2_size;
    }
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

    // redirect const pointers to extern variables we need
    if (strcmp(grid_name,"source") == 0)
    {
        grid_center_lat     = grid1_center_lat;
        grid_center_lon     = grid1_center_lon;
        grid_centroid_lat   = grid1_centroid_lat;
        grid_centroid_lon   = grid1_centroid_lon;
        grid_corner_lat     = grid1_corner_lat;
        grid_corner_lat     = grid1_corner_lon;
        grid_mask           = grid1_mask;
        grid_corners        = grid1_corners;
    }
    else if (strcmp(grid_name, "destination") == 0)
    {
        grid_center_lat     = grid2_center_lat;
        grid_center_lon     = grid2_center_lon;
        grid_centroid_lat   = grid2_centroid_lat;
        grid_centroid_lon   = grid2_centroid_lon;
        grid_corner_lat     = grid2_corner_lat;
        grid_corner_lat     = grid2_corner_lon;
        grid_mask           = grid2_mask;
        grid_corners        = grid2_corners;
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
            if (lat == 0)
            {
                index_of_neighbor_cells[cell][0] = curr + left;
                index_of_neighbor_cells[cell][1] = south + left;
                index_of_neighbor_cells[cell][2] = south + lon;
                index_of_neighbor_cells[cell][3] = south + right;
                index_of_neighbor_cells[cell][4] = curr + right;
            }
            else if (lat == num_of_different_lats - 1)
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

void Neighborhood::calculate_gradient_latlon (double * result_lat, double * result_lon, const double * field_value, int grid_size)
{
    for (int cell = 0; cell < grid_size; cell ++)
    {
        calculate_derivative_latlon_least_square_method(result_lat[cell], result_lon[cell], cell, grid_centroid_lat, grid_centroid_lon, 
            num_of_neighbor_cells[cell], index_of_neighbor_cells[cell], grid_center_lat, grid_center_lon, 
            field_value);
    }
}

void Neighborhood::calculate_gradient_lat (double * result_lat, const double * field_value, int grid_size)
{
    for (int cell = 0; cell < grid_size; cell ++)
    {
        calculate_derivative_lat_least_square_method(result_lat[cell], cell, grid_centroid_lat, grid_centroid_lon, 
            num_of_neighbor_cells[cell], index_of_neighbor_cells[cell], grid_center_lat, grid_center_lon, 
            field_value);
    }
}

void Neighborhood::calculate_gradient_lon (double * result_lon, const double * field_value, int grid_size)
{
    for (int cell = 0; cell < grid_size; cell ++)
    {
        calculate_derivative_lon_least_square_method(result_lon[cell], cell, grid_centroid_lat, grid_centroid_lon, 
            num_of_neighbor_cells[cell], index_of_neighbor_cells[cell], grid_center_lat, grid_center_lon, 
            field_value);
    }
}
