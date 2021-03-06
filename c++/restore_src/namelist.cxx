/* This file define extern variables */

#include "namelist.h"

char * grid1_file = "T42.nc";
char * grid2_file = "G128.nc";
//char * grid2_file = "POP43.nc";
char * interp_file = "SSQoutput1.nc";
char * map_name = "T42-to-G128";
//char * map_name = "T42-to-POP43";
char * map_method = "conservative";
char * normalize_opt = "destarea";
int output_opt = SCRIP_CONVENTION;
char *restrict_type = "latitude";
int num_srch_bins = 300;
bool luse_grid1_area = false;
bool luse_grid2_area = false;

bool luse_grid_centers;         // defined in sc2rip.cxx
