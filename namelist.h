#ifndef _NAMELIST_
#define _NAMELIST_ 1

extern char *grid1_file, *grid2_file;   // netCDF file name
extern char *interp_file1, *interp_file2;   // remap wts file name
extern char *map_name;          // name of remapping
extern char *map_method;        // remapping algorithm name
extern char *normalize_opt;     // normalize option (conserv only)
extern char *output_opt;        // remap result file name
extern char *restrict_type;     // type of bins to use
extern int num_srch_bins;       // num of bins for restricted srch
extern bool luse_grid1_area;    // use area from grid file
extern bool luse_grid2_area;    // use area from grid file

extern bool luse_grid_centers;  // use centers for bounding box
#endif
