#ifndef _GRIDS_H_
#define _GRIDS_H_ 1

#include <string.h>
#include <netcdf.h>

#include "kinds.h"
#include "iounits.h"
#include "nc_error.h"
#include "great_circle_base.h"

#define _2D_BOUND_BOX_ 1
#define _3D_BOUND_BOX_ (1 - _2D_BOUND_BOX_)
#if _2D_BOUND_BOX_
    #define BOUNDBOX_SIZE 4
#else
    #define BOUNDBOX_SIZE 6
#endif
#define     ERR     netcdf_error_handler(ncstat);

/* external variables */
extern Grid * grid_src, * grid_dst;     // source and destination Grids 

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
