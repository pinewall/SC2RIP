#ifndef _NAMELIST_
#define _NAMELIST_ 1

extern bool luse_grid_centers;  //use centers for bounding box
extern bool luse_grid1_area;    //use area from grid file
extern bool luse_grid2_area;    //use area from grid file
extern char *restrict_type;     // type of bins to use
extern int num_srch_bins;       // num of bins for restricted srch

#endif
