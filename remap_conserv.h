#ifndef _REMAP_CONSERV_
#define _REMAP_CONSERV_ 1

#include "kinds.h"
#include "constants.h"
#include "timers.h"
#include "grids.h"
#include "remap_vars.h"

#include "intersection.h"
#include "line_integral.h"
#include "store_link_cnsrv.h"

// threshhold precompile constants
#define NORTH_THRESH (2.00)
#define SOUTH_THRESH (-2.00)
#define INTEGRATE_AROUND_SRC_GRID 1
#define INTEGRATE_AROUND_DST_GRID 2

// extern variables
extern int num_srch_cells;      // num cells in restricted search arrays
extern int *srch_add;           // global address of cells in srch arrays
extern double *srch_corner_lat;    // lat of each corner of srch cells
extern double *srch_corner_lon;    // lon of each corner of srch cells

/** this routine traces the perimeters of every grid cell on each
 *  grid checking for intersections with the other grid and computing
 *  line integrals for each subsegment
 **/
void remap_conserv();

/* integrate around each cell on one grid */
int conserv_sweep(int choice, double *grid1_controid_lat, double *grid1_centroid_lon, double *grid2_centroid_lat, double *grid2_centroid_lon);

#endif
