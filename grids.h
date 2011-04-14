#ifndef _GRIDS_H_
#define _GRIDS_H_ 1

#include <string.h>
#include <netcdf.h>

#include "kinds.h"
#include "constants.h"
#include "iounits.h"
#include "nc_error.h"

#define _2D_BOUND_BOX_ 1
#define _3D_BOUND_BOX_ (1 - _2D_BOUND_BOX_)
#if _2D_BOUND_BOX_
    #define BOUNDBOX_SIZE 4
#else
    #define BOUNDBOX_SIZE 6
#endif
#define deg2rad (PI/180)
#define     ERR     netcdf_error_handler(ncstat);

extern unsigned int grid1_size, grid2_size;  // total points on each grid
extern unsigned int grid1_rank, grid2_rank;  // rank of each grid
/* NOTE: Using variable corner number */
extern unsigned int grid1_corners_max, grid2_corners_max;    // max number of corners for total grid
extern unsigned int *grid1_corners, *grid2_corners;  //number of corners for each cell
extern int *grid1_dims, *grid2_dims;    // size of each grid dimension
extern char *grid1_name, *grid2_name;   // name for each grid
extern char *grid1_units, *grid2_units;  // units for grid (degrees or radians)

/* grid coordinates and masks */
extern bool *grid1_mask, *grid2_mask;   // flag which cells participate
extern double *grid1_center_lat, *grid1_center_lon; //lat/lon coordinates for each grid center in radians
extern double *grid2_center_lat, *grid2_center_lon;
extern double *grid1_corner_lat, *grid1_corner_lon;
extern double *grid2_corner_lat, *grid2_corner_lon;
extern double *grid1_area, *grid2_area;     //total area of each grid cell
extern double *grid1_frac, *grid2_frac;     //fractional area of grid cells participating in remapping
extern double *grid1_area_in, *grid2_area_in;   // total area of each grid cell from netCDF file
/* NOTE: Using three-dimention bound_box to avoid pole problem */
extern double *grid1_bound_box, *grid2_bound_box;    //lat/lon bound box for use in restricting grid searches

/* bins for restricting searches */
extern int *bin_addr1, *bin_addr2;  // min|max index for each cell in bin
extern double *bin_lats, *bin_lons; // min|max lat/lon for each srch bin


/** functions **/
// grid init
extern void grid_init(char *grid1_file, char* grid2_file);
// grid init for source grid
extern void grid_init_src(char *grid_src_file);
// grid init for destination grid
extern void grid_init_dst(char *grid_dst_file);
// calculate bounding box for each grid
void grid_cal_boundbox(double *boundbox, bool *grid_mask, int grid_size, double *center_lat, double *center_lon, double *corner_lat, double *corner_lon, unsigned int *grid_corners, unsigned int grid_corners_max);
// init search bins
void grid_srch_bin_init();
// assign search bin address, that is min|max index of cell
void grid_assign_srch_bin(double *boundbox, int *addr, int grid_size);
// grid debug
void grid_debug();
#endif
