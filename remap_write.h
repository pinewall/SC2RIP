#ifndef _REAMP_WRITE_
#define _REMAP_WRITE_ 1

// call correct output routine based on output format choice (scrip or csm)
void write_remap(char *map_name, char *interp_file, int output_opt);

// writes remap data to a netCDF file using SCRIP conventions
void write_remap_scrip(char *map_name, char *interp_file, int direction);

// writes remap data to a netCDF file using NCAR-CSM conventions
void write_remap_csm(char *map_name, char *interp_file, int direction);

/** sorts address and weight arrays based on the 
  * destination address with the source address as a secondary 
  * sorting criterion. the method is a standard heap sort.
  **/
void sort_add(int *add1, int *add2, double *weights);

#endif
