#include "remap_vars.h"
#include "utils.h"
#include "grids.h"

#include <stdio.h>

// declear extern variables
int num_links_map;          // actual number of links for remapping
int max_links_map;          // current size of link arrays
int num_wts;                // num of weights used in remapping
int map_type;               // identifier for remapping method
int norm_opt;               // option for noramlization (conserv only)
int resize_increment;       // defalut amout to increase array size

int *grid1_add_map;         // grid1 address for each link in remapping
int *grid2_add_map;         // grdi2 address for each link in remapping

double *wts_map;            // map weights for each link (wts, links)

// init remapping variables
void init_remap_vars()
{
    // determine the number of weights
    switch (map_type)
    {
        case MAP_TYPE_CONSERV:
          num_wts = 3;
          break;
        case MAP_TYPE_BILINEAR:
          num_wts = 1;
          break;
        case MAP_TYPE_BICUBIC:
          num_wts = 4;
          break;
        case MAP_TYPE_DISTWGT:
          num_wts = 1;
          break;
    default:
          num_wts = 1;
    }

    /** initialize num_links and set max_links to four times the largest
     *  of the destination grid sizes initally (can be changed later).
     *  set a default resize increment to increase the size fo link
     *  arrays if the number of links exceeds the inital size
     **/
    int grid1_size = 8192;
    int grid2_size = 32768;
    num_links_map = 0;
    max_links_map = 8 * grid2_size;
    //max_links_map = 4 * grid2_size;
    resize_increment = 0.1 * MAX(grid1_size, grid2_size);

    // allocate address and weight arrays for mapping
    grid1_add_map = new int [max_links_map];
    grid2_add_map = new int [max_links_map];
    wts_map = new double [max_links_map * num_wts];
    memset(grid1_add_map, -1, sizeof(int) * max_links_map);
    memset(grid2_add_map, -1, sizeof(int) * max_links_map);
    memset(wts_map, 0.0, sizeof(double) * max_links_map * num_wts);
}

// resize remapping arrays by increasing(decreasing) the max_links by increment
void resize_remap_vars(int increment)
{
    // allocate temporaries to hold original values
    int mxlinks = max_links_map;
    int *add1_tmp = grid1_add_map;
    int *add2_tmp = grid2_add_map;
    double *wts_tmp = wts_map;
    
    // allocate new arrays
    max_links_map += increment;     // new size of links
    grid1_add_map = new int [max_links_map];
    grid2_add_map = new int [max_links_map];
    wts_map = new double [max_links_map * num_wts];

    // move back original values; using System Call to decrease overhead
    mxlinks = MIN(mxlinks, max_links_map);      // amout to move
    memcpy(grid1_add_map, add1_tmp, mxlinks * sizeof(int));
    memcpy(grid2_add_map, add2_tmp, mxlinks * sizeof(int));
    memcpy(wts_map, wts_tmp, mxlinks * sizeof(double) * num_wts);

    // delete original values
    delete [] add1_tmp;
    delete [] add2_tmp;
    delete [] wts_tmp;

}

void finalize_remap_vars()
{
    delete [] grid1_add_map;
    delete [] grid2_add_map;
    delete [] wts_map;
}
