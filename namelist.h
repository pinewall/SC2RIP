#ifndef _NAMELIST_
#define _NAMELIST_ 1

#define SCRIP_CONVENTION 1
#define CSM_CONVENTION 2

extern char *grid1_file;        // netCDF file name
extern char *grid2_file;        // netCDF file name
extern char *interp_file;       // remap wts file name
extern char *map_name;          // name of remapping
extern char *map_method;        // remapping algorithm name
extern char *normalize_opt;     // normalize option (conserv only)
extern int output_opt;          // remap result convention (SCRIP or NCAR-CSM)
extern char *restrict_type;     // type of bins to use
extern int num_srch_bins;       // num of bins for restricted srch
extern bool luse_grid1_area;    // use area from grid file
extern bool luse_grid2_area;    // use area from grid file

extern bool luse_grid_centers;  // use centers for bounding box
#endif
