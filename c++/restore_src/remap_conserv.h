#ifndef _REMAP_CONSERV_H_
#define _REMAP_CONSERV_H_ 1

#include "kinds.h"
#include "constants.h"
#include "timers.h"
#include "grids.h"
#include "remap_vars.h"

#include "intersection.h"
#include "line_integral.h"
#include "store_link_cnsrv.h"

#include "namelist.h"
#include "utils.h"
#include <memory.h>

#define INTEGRATE_AROUND_SRC_GRID 1
#define INTEGRATE_AROUND_DST_GRID 2

#define SRCH_SIZE 1024
#define SRCH_CORNER_MAX 8

extern double *grid2_centroid_lat; // centroid coords on each grid
extern double *grid2_centroid_lon;
extern double *grid1_centroid_lat;
extern double *grid1_centroid_lon;

/** this routine traces the perimeters of every grid cell on each
 *  grid checking for intersections with the other grid and computing
 *  line integrals for each subsegment
 **/
void remap_conserv();

/* integrate around each cell on one grid */
int conserv_sweep(int choice, double *grid1_controid_lat, double *grid1_centroid_lon, double *grid2_centroid_lat, double *grid2_centroid_lon);

void finalize_remap_conserv();

// finalize intersection
void finalize_intersection();
#endif
