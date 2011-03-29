#ifndef _REMAP_VARS_
#define _REMAP_VARS_ 1

#include "kinds.h"
#include "constants.h"
#include "grids.h"

#define NORM_OPT_NONE 1
#define NORM_OPT_DSTAREA 2
#define NORM_OPT_FRACAREA 3

#define MAP_TYPE_CONSERV 1
#define MAP_TYPE_BILINEAR 2
#define MAP_TYPE_BICUBIC 3
#define MAP_TYPE_DISTWGT 4

extern int max_links_map;      // current size of link arrays
extern int num_links_map;      // actual number of links for remapping
extern int num_wts;             // num of weights used in remapping
extern int map_type;            // identifier for remapping method
extern int norm_type;           // option for normalization (conserv only)
extern int resize_increment;    // default amout to increase array size

extern int *grid1_add_map;     // grid1 address for each link in mapping 1
extern int *grid2_add_map;     // grid2 address for each link in mapping 1

extern double *wts_map;        // map weights for each link (wts, links)

void init_remap_vars();         // init remapping variables
void resize_remap_vars(int increment);      // resize remapping variable scale

#endif
