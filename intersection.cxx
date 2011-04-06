#include "intersection.h"

// define extern variables
int last_loc;               // save location when crossing threshold
bool lthresh = false;       // flags segments crossing threshold bndy
double intrsct_lat_off;     // latitude coords offset for next search
double intrsct_lon_off;     // longitude coords offset for next search

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
        double *begseg, bool lbegin, bool lrevers)
{
    // local variables
    int n, next_n, cell, srch_corners, pole_loc;
    bool loutside;              // flags points outside grid
    double lon1, lon2;          // local longitude variables for segment
    double lat1, lat2;          // local latitude variables for segment
    double grdlon1, grdlon2;    // local longitude variables for grid cell
    double grdlat1, grdlat2;    // local latitude variables for grid cell
    double vec1_lat, vec1_lon;  // vectors and cross product used during grid search
    double cross_product;
    double eps, offset;         // small offset away from intersect
    double s1, s2, determ;      // variables used for linear solve to find intersection
    double mat1, mat2, mat3, mat4, rhs1, rhs2;

    // initialize defaults, flags, etc
}
