#include "intersection.h"


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
    location = -1;
    intrsct_lat = 0.0;
    intrsct_lon = 0.0;
    lcoinc = true;
    return;
}
