#ifndef _INTERSECTION_H
#define _INTERSECTION_H 1

#include <stdio.h>
#include "utils.h"
#include "constants.h"

// threshhold precompile constants
#define NORTH_THRESH (2.00)
#define SOUTH_THRESH (-2.00)
#define tiny (1e-14)

// declare extern variables
extern unsigned int num_srch_cells;      // num cells in restricted search arrays
extern int *srch_add;           // global address of cells in srch arrays
extern double *srch_corner_lat;    // lat of each corner of srch cells
extern double *srch_corner_lon;    // lon of each corner of srch cells
extern int last_loc;        // save location when crossing threshold
extern bool lthresh;        // flags segments crossing threshold bndy
extern double intrsct_lat_off;     // latitude coords offset for next search
extern double intrsct_lon_off;     // longitude coords offset for next search

// finalize intersection
void finalize_intersection();

/** this routine finds the next intersection of a destination grid
  * line with the line segment given by beglon, endlon, etc.
  * a coincidence flag is returnd if the segment is entirely
  * coincident with an ocean grid line. the cells in which to search
  * for an intersection must have already been restricted in the
  * calling routine
  **/
void intersection(int &location, 
        double &intrsct_lat, double &intrsct_lon, bool &lcoinc, 
        double beglat, double beglon, double endlat, double endlon, 
        double *begseg, bool lbegin, bool lrevers);

void pole_intersection(int &location, 
        double &intrsct_lat, double &intrsct_lon, bool &lcoinc, 
        double beglat, double beglon, double endlat, double endlon, 
        double *begseg, bool lbegin, bool lrevers);

/** this routine determin whether a point --given by lat/lon pair-- 
  * in a cell with several corners --corner_num-- and lat/lon coords 
  * @param[corner_lat, corner_lon]: pointer to header of corner info
  * @param[corner_num]: number of corners
  * @param[objlat, objlon]: point to judge whether in district
  * @param[assistlat, assistlon]: assistant when obj on edge of district
  **/
bool is_point_in(double *corner_lat, double *corner_lon, 
        int corner_num, double objlat, double objlon, 
        double assistlat, double assistlon, bool &lcoinc, bool lrevers);

#endif
