#include "intersection.h"
#include "debug.h"

#include <math.h>

// define extern variables
int             last_loc = -1;              // save location when crossing threshold
bool            lthresh = false;            // flags segments crossing threshold bndy
bool            luse_last = false;
double          intrsct_lat_off;            // latitude coords offset for next search
double          intrsct_lon_off;            // longitude coords offset for next search
double          intrsct_x;
double          intrsct_y;

int             avoid_pole_count = 0;
double          avoid_pole_offset = tiny;

unsigned int    num_srch_cells;
int     *       srch_add = NULL;            // global address of cells in srch arrays
double  *       srch_corner_lat = NULL;     // lat of each corner of srch cells
double  *       srch_corner_lon = NULL;     // lon of each corner of srch cells
double  *       srch_corner_x = NULL;
double  *       srch_corner_y = NULL;

/** this routine finds the next intersection of a destination grid
  * line with the line segment given by beglon, endlon, etc.
  * a coincidence flag is returnd if the segment is entirely
  * coincident with an ocean grid line. the cells in which to search
  * for an intersection must have already been restricted in the
  * calling routine
  **/
void ll_ll_intersection(int &location, 
        double &intrsct_lat, double &intrsct_lon, bool &lcoinc, 
        double beglat, double beglon, double endlat, double endlon, 
        double *begseg, bool lbegin, bool lrevers)
{
    // local variables
    double vec1_lat, vec1_lon;      // first vector
    double vec2_lat, vec2_lon;      // second vector
    double cp;                      // cross product value
    int n, next_n, cell, srch_corners, pole_loc;
    bool loutside;              // flags points outside grid
    double lon1, lon2;          // local longitude variables for segment
    double lat1, lat2;          // local latitude variables for segment
    double grdlon1, grdlon2;    // local longitude variables for grid cell
    double grdlat1, grdlat2;    // local latitude variables for grid cell
    double eps, offset;         // small offset away from intersect
    double s1, s2, determ;      // variables used for linear solve to find intersection
    double mat1, mat2, mat3, mat4, rhs1, rhs2;

    // initialize defaults, flags, etc
    location = -1;
    lcoinc = false;
    intrsct_lat = endlat;
    intrsct_lon = endlon;
    if (num_srch_cells == 0)    // no work to do
        return;

    // pole transformation control
    if (beglat > NORTH_THRESH || beglat < SOUTH_THRESH)
    {
        if (lthresh)
            location = last_loc;

        pole_intersection(location, intrsct_lat, intrsct_lon, lcoinc, lthresh, beglat, beglon, endlat, endlon, begseg, lbegin, lrevers);

        if (lthresh)
        {
            last_loc = location;
            intrsct_lat_off = intrsct_lat;
            intrsct_lon_off = intrsct_lon;
            printf("cross threshold @ location = %d\n", location);
            printf("\tintersection @ (%3.6f, %3.6f)\n", intrsct_lat, intrsct_lon);
        }
        return;
    }
    
    // normal case of intersections
    loutside = false;
    if (lbegin)                 // using begin of full segment
    {
        lat1 = beglat;
        lon1 = beglon;
    }
    else
    {
        lat1 = intrsct_lat_off; // using shift away from intersection
        lon1 = intrsct_lon_off;
    }
    lat2 = endlat;
    lon2 = endlon;

    if ((lon2 - lon1) > THREE * PIH)
        lon2 -= PI2;
    else if ((lon2 - lon1) < - THREE * PIH)
        lon2 += PI2;

    //printf("init: %f %f %f %f\n", lat1, lon1, lat2, lon2);
    s1 = ZERO;                  // parameter controling endpoints shift

    // search for location of this segment in ocean grid using cross product method to determine whether a point is enclosed by a cell
    bool quit = false;
    srch_corners = 4;

    // otherwise normal search algorithm
    while (!quit)
    {

        // if last segment crossed threshold, use that location
        if (lthresh)
        {
            for (cell = 0; cell < num_srch_cells; cell ++)
            {
                if (srch_add[cell] == last_loc)
                {
                    location = last_loc;
                    eps = tiny;
                    quit = true;
                    break;
                }
            }
        }

        if (quit)
            break;

        for (cell = 0; cell < num_srch_cells; cell ++)
        {
            for (n = 0; n < srch_corners; n++)
            {
                next_n = (n + 1) % srch_corners;
                vec1_lat = srch_corner_lat[cell*srch_corners + next_n] - srch_corner_lat[cell*srch_corners + n];
                vec1_lon = srch_corner_lon[cell*srch_corners + next_n] - srch_corner_lon[cell*srch_corners + n];
                vec2_lat = lat1 - srch_corner_lat[cell*srch_corners + n];
                vec2_lon = lon1 - srch_corner_lon[cell*srch_corners + n];

                // if endpoint coincident with vertex, offset the endpoint
                if (zero(vec2_lat) && zero(vec2_lon))
                {
                    lat1 = lat1 +  1e-10 * (lat2 - lat1);
                    lon1 = lon1 + 1e-10 * (lon2 - lon1);
                    vec2_lat = lat1 - srch_corner_lat[cell*srch_corners + n];
                    vec2_lon = lon1 - srch_corner_lon[cell*srch_corners + n];
                }   

                // check for 0.PI2 crossing
                check_longitude(vec1_lon, - PI, PI);
                check_longitude(vec2_lon, - PI, PI);
                //printf("vec lon: %f %f\n", vec1_lon, vec2_lon);

                /* cross product method
                    using (x1,y1,0)X(x2,y2,0) = (0,0,x1y2-y1x2)
                    to decide direction */
                cp = vec1_lon * vec2_lat - vec2_lon * vec1_lat;
                //printf("cp 1: %f\n", cp);

            /** if the cross product for a side is zero, the point
             *     lies exactly on the side or the side is degenerate
             *     (zero length). if degenerate, set the cross 
             *     product to a positive number. otherwise perform
             *     another cross product between the side and the
             *     segment itself.
             * if this cross product is also zero, the line is
             *     coincident with the cell boundary - perform the
             *     dot product and only choose the cell if the dot
             *     product is positive (parallel vs anti-parallel).
            **/
                if (zero(cp))
                {
                    if (nonzero(vec1_lat) || nonzero(vec1_lon))
                    {
                        vec2_lat = lat2 - lat1;
                        vec2_lon = lon2 - lon1;

                        check_longitude(vec2_lon, - PI, PI);

                        cp = vec1_lon * vec2_lat - vec2_lon * vec1_lat;
                    }
                    else
                        cp = ONE;

                    if (zero(cp))
                    {
                        lcoinc = true;
                        cp = vec1_lon * vec2_lon + vec1_lat * vec2_lat;
                        if (lrevers)
                            cp = -cp;
                    }
                }   

                //printf("cp 2: %f\n", cp);
                // if cp is less than zero, this cell doesn't work
                if (cp < ZERO)
                    break;
            }
            
            if (n == srch_corners)
            {
                location = srch_add[cell];
                // if the beginning of this segment was outside the grid, ivert the segment so the intersection found will be the first intersection with the grid
                if (loutside)
                {
                    lat2 = beglat;
                    lon2 = beglon;
                    location = -1;
                    eps = - tiny;
                }
                else
                    eps = tiny;
                quit = true;
#if _DEBUG_SRCH_
                printf("found in %d\n", location);
#endif
                break;
            }
            // otherwise move on to the next cell
        }
        if (quit)
            break;
        
        // if still no cell found, the point lies outside the grid. take some baby steps along the segment to see if any part of the segment lies inside the grid
        loutside = true;
        //printf("loutside\n");
        s1 += 0.5;
        lat1 = beglat + s1 * (endlat - beglat);
        lon1 = beglon + s1 * (lon2 - beglon);

        // reached the end of the segment and still outside the grid, return no intersection
        if (s1 >= ONE)
            return;
    }


    /** now that a cell found, search for the next intersection.
      * loop over sides of the cell to find intersection with side
      * must check all sides for coincidences or intersections
     **/
    for (n = 0; n < srch_corners; n++)
    {
        next_n = (n + 1) % srch_corners;

        grdlon1 = srch_corner_lon[cell * srch_corners + n];
        grdlon2 = srch_corner_lon[cell * srch_corners + next_n];
        grdlat1 = srch_corner_lat[cell * srch_corners + n];
        grdlat2 = srch_corner_lat[cell * srch_corners + next_n];

        // set up linear system to solve for intersection

        mat1 = lat2 - lat1;
        mat2 = grdlat1 - grdlat2;
        mat3 = lon2 - lon1;
        mat4 = grdlon1 - grdlon2;
        rhs1 = grdlat1 - lat1;
        rhs2 = grdlon1 - lon1;

        check_longitude(mat3, -PI, PI);
        check_longitude(mat4, -PI, PI);
        check_longitude(rhs2, -PI, PI);

        determ = mat1*mat4 - mat2*mat3;

        /** if the determinant is zero, the segments are either
          *     parallel or coincident. Coincidences were detected
          *     above so do nothing
          * if the determinant is non-zero, solve for the linear
          *     parameters s for the intersection point on each line
          *     segment
          * if 0 < s1, s2 < 1 then the segment intersects with this side
          *     return the point of intersection (adding a small 
          *     number so the intersection is off the grid (line)
          */
        if (ABS(determ) > 1e-30)
        {
            s1 = (rhs1*mat4 - mat2*rhs2) / determ;
            s2 = (mat1*rhs2 - rhs1*mat3) / determ;

            if (s2 >= ZERO && s2 <= ONE &&
                s1 >  ZERO && s1 <= ONE)
            {
                // recompute intersection based on full segment so the intersections are consistent for both sweeps
                if (! loutside)
                {
                    mat1 = lat2 - begseg[0];
                    mat3 = lon2 - begseg[1];
                    rhs1 = grdlat1 - begseg[0];
                    rhs2 = grdlon1 - begseg[1];
                }
                else
                {
                    mat1 = begseg[0] - endlat;
                    mat3 = begseg[1] - endlon;
                    rhs1 = grdlat1 - endlat;
                    rhs2 = grdlon1 - endlon;
                }

                check_longitude(mat3);
                check_longitude(rhs2);

                determ = mat1*mat4 - mat2*mat3;

                // re-check non-zero of determ
                if (nonzero(determ))
                {
                    s1 = (rhs1*mat4 - mat2*rhs2) / determ;
                    s2 = (mat1*rhs2 - rhs1*mat3) / determ;

                    offset = s1 + eps / determ;
                    if (offset > ONE)
                        offset = ONE;

                    if (! loutside)
                    {
                        intrsct_lat = begseg[0] + mat1*s1;
                        intrsct_lon = begseg[1] + mat3*s1;
                        intrsct_lat_off = begseg[0] + mat1*offset;
                        intrsct_lon_off = begseg[1] + mat3*offset;
                    }
                    else
                    {
                        intrsct_lat = endlat + mat1*s1;
                        intrsct_lon = endlon + mat3*s1;
                        intrsct_lat_off = endlat + mat1*offset;
                        intrsct_lon_off = endlon + mat3*offset;
                    }
                    break;
                }           /* determ nonzero */
            }               /* s1 s2 */
        }                   /* determ nonzero */
        // no intersection this side, move on to next side
    }                       /* intrsct loop */

    /*  if the segment crosses a pole threshold, reset the intersection
     *  to be the threhold latitude. only check if this was not a
     *  threshold segment since sometimes coordinate transform can end
     *  up on other side of threshold again
     */
    if (lthresh)
    {
        if (intrsct_lat < NORTH_THRESH && intrsct_lat > SOUTH_THRESH)
        {
            lthresh = false;
        }
    }
    else if (lat1 > ZERO && intrsct_lat > NORTH_THRESH)
    {
        intrsct_lat = NORTH_THRESH + tiny;
        intrsct_lat_off = NORTH_THRESH + eps * mat1;
        s1 = (intrsct_lat - begseg[0]) / mat1;
        intrsct_lon =       begseg[1] + s1 * mat3;
        intrsct_lon_off =   begseg[1] + (s1 + eps) * mat3;
        last_loc = location;
        lthresh = true;
    }
    else if (lat1 < ZERO && intrsct_lat < SOUTH_THRESH)
    {
        intrsct_lat = SOUTH_THRESH - tiny;
        intrsct_lat_off = SOUTH_THRESH + eps * mat1;
        s1 = (intrsct_lat - begseg[0]) / mat1;
        intrsct_lon =       begseg[1] + s1 * mat3;
        intrsct_lon_off =   begseg[1] + (s1 + eps) * mat3;
        last_loc = location;
        lthresh = true;
    }
}

void pole_intersection(int &location, 
        double &intrsct_lat, double &intrsct_lon, bool &lcoinc, bool &lthresh,
        double beglat, double beglon, double endlat, double endlon, 
        double *begseg, bool lbegin, bool lrevers)
{
    // local variables
    int n, next_n, cell, corner, srch_corners, pole_loc;

    bool loutside;                  // flags points outside grid

    double PI4, rns;                // for pole conversion
    double vec1_x, vec1_y;      // first vector
    double vec2_x, vec2_y;      // second vector
    double cp;                      // cross product value
    double y1, y2;              // local longitude variables for segment
    double x1, x2;              // local latitude variables for segment
    double grdy1, grdy2;        // local longitude variables for grid cell
    double grdx1, grdx2;        // local latitude variables for grid cell
    double eps, offset;             // small offset away from intersect
    double s1, s2, determ;          // variables used for linear solve to find intersection
    double mat1, mat2, mat3, mat4, rhs1, rhs2;
    double begsegx, begsegy;
    double begx, begy, endx, endy;

    // initialize defaults, flags, etc
    if (!lthresh)
        location = -1;
    lcoinc = false;
    intrsct_lat = endlat;
    intrsct_lon = endlon;
    
    loutside = false;
    s1 = ZERO;                  // parameter controling endpoints shift

    /* convert coordinates */
    if (beglat > ZERO)
    {
        PI4 = PI / 4;
        rns = ONE;
    }
    else
    {
        PI4 = - PI / 4;
        rns = - ONE;
    }

    if (luse_last)
    {
        x1 = intrsct_x;
        y1 = intrsct_y;
    }
    else
    {
        x1 = rns * TWO * sin(PI4 - HALF * beglat) * cos(beglon);
        y1 =       TWO * sin(PI4 - HALF * beglat) * sin(beglon);
        luse_last = true;
    }
    x2 = rns * TWO * sin(PI4 - HALF * endlat) * cos(endlon);
    y2 =       TWO * sin(PI4 - HALF * endlat) * sin(endlon);

    for (cell = 0; cell < num_srch_cells; cell ++)
    {
        for (corner = 0; corner < 4; corner ++)
        {
            srch_corner_x[cell*4 + corner] = rns * TWO * sin(PI4 - HALF * srch_corner_lat[cell*4 + corner]) * cos(srch_corner_lon[cell*4 + corner]);
            srch_corner_y[cell*4 + corner] =       TWO * sin(PI4 - HALF * srch_corner_lat[cell*4 + corner]) * sin(srch_corner_lon[cell*4 + corner]);
        }
    }

    begx = x1;
    begy = y1;
    endx = x2;
    endy = y2;
    begsegx = rns * TWO * sin(PI4 - HALF * begseg[0]) * cos(begseg[1]);
    begsegy =       TWO * sin(PI4 - HALF * begseg[0]) * sin(begseg[1]);
    intrsct_x = endx;
    intrsct_y = endy;

    // search for location of this segment in ocean grid using cross product method to determine whether a point is enclosed by a cell
    bool quit = false;
    srch_corners = 4;

    while (!quit)
    {
        // if last segment crossed threshold, use that location
        if (lthresh)
        {
            for (cell = 0; cell < num_srch_cells; cell ++)
            {
                if (srch_add[cell] == last_loc)
                {
                    eps = tiny;
                    quit = true;
                    break;
                }
            }
        }

        if(quit)
            break;

        // otherwise normal search algorithm
        for (cell = 0; cell < num_srch_cells; cell ++)
        {
            for (n = 0; n < srch_corners; n++)
            {
                next_n = (n + 1) % srch_corners;
                vec1_x = srch_corner_x[cell*srch_corners + next_n] - srch_corner_x[cell*srch_corners + n];
                vec1_y = srch_corner_y[cell*srch_corners + next_n] - srch_corner_y[cell*srch_corners + n];
                vec2_x = x1 - srch_corner_x[cell*srch_corners + n];
                vec2_y = y1 - srch_corner_y[cell*srch_corners + n];

                // if endpoint coincident with vertex, offset the endpoint
                if (zero(vec2_x) && zero(vec2_y))
                {
                    x1 = x1 + 1e-10 * (x2 - x1);
                    y1 = y1 + 1e-10 * (y2 - y1);
                    vec2_x = x1 - srch_corner_x[cell*srch_corners + n];
                    vec2_y = y1 - srch_corner_y[cell*srch_corners + n];
                }   

                cp = vec1_x * vec2_y - vec2_x * vec1_y;

            /** if the cross product for a side is zero, the point
             *     lies exactly on the side or the side is degenerate
             *     (zero length). if degenerate, set the cross 
             *     product to a positive number. otherwise perform
             *     another cross product between the side and the
             *     segment itself.
             * if this cross product is also zero, the line is
             *     coincident with the cell boundary - perform the
             *     dot product and only choose the cell if the dot
             *     product is positive (parallel vs anti-parallel).
            **/
                if (zero(cp))
                {
                    if (nonzero(vec1_x) || nonzero(vec1_y))
                    {
                        vec2_x = x2 - x1;
                        vec2_y = y2 - y1;
                        cp = vec1_x * vec2_y - vec2_x * vec1_y;
                    }
                    else
                        cp = ONE;

                    if (zero(cp))
                    {
                        lcoinc = true;
                        cp = vec1_x * vec2_x + vec1_y * vec2_y;
                        if (lrevers)
                            cp = -cp;
                    }
                }   

                // if cp is less than zero, this cell doesn't work
                if (cp < ZERO)
                    break;
            }
            
            if (n == srch_corners)
            {
                location = srch_add[cell];
                // if the beginning of this segment was outside the grid, ivert the segment so the intersection found will be the first intersection with the grid
                if (loutside)
                {
                    x2 = begx;
                    y2 = begy;
                    location = -1;
                    eps = - tiny;
                }
                else
                    eps = tiny;
                quit = true;
#if _DEBUG_SRCH_
                printf("found in %d\n", location);
#endif
                break;
            }
            // otherwise move on to the next cell
        }
        if (quit)
            break;
        
        // if still no cell found, the point lies outside the grid. take some baby steps along the segment to see if any part of the segment lies inside the grid
        loutside = true;
        //printf("loutside\n");
        s1 += 0.001;
        x1 = begx + s1 * (x2 - begx);
        y1 = begy + s1 * (y2 - begy);

        // reached the end of the segment and still outside the grid, return no intersection
        if (s1 >= ONE)
        {
            luse_last = false;
            return;
        }
    }


    /** now that a cell found, search for the next intersection.
      * loop over sides of the cell to find intersection with side
      * must check all sides for coincidences or intersections
     **/
    for (n = 0; n < srch_corners; n++)
    {
        next_n = (n + 1) % srch_corners;

        grdy1 = srch_corner_y[cell * srch_corners + n];
        grdy2 = srch_corner_y[cell * srch_corners + next_n];
        grdx1 = srch_corner_x[cell * srch_corners + n];
        grdx2 = srch_corner_x[cell * srch_corners + next_n];

        // set up linear system to solve for intersection

        mat1 = x2 - x1;
        mat2 = grdx1 - grdx2;
        mat3 = y2 - y1;
        mat4 = grdy1 - grdy2;
        rhs1 = grdx1 - x1;
        rhs2 = grdy1 - y1;

        determ = mat1*mat4 - mat2*mat3;

        /** if the determinant is zero, the segments are either
          *     parallel or coincident. Coincidences were detected
          *     above so do nothing
          * if the determinant is non-zero, solve for the linear
          *     parameters s for the intersection point on each line
          *     segment
          * if 0 < s1, s2 < 1 then the segment intersects with this side
          *     return the point of intersection (adding a small 
          *     number so the intersection is off the grid (line)
          */
        if (ABS(determ) > 1e-30)
        {
            s1 = (rhs1*mat4 - mat2*rhs2) / determ;
            s2 = (mat1*rhs2 - rhs1*mat3) / determ;

            if (s2 >= ZERO && s2 <= ONE &&
                s1 >  ZERO && s1 <= ONE)
            {
                // recompute intersection based on full segment so the intersections are consistent for both sweeps
                if (! loutside)
                {
                    mat1 = x2 - begsegx;
                    mat3 = y2 - begsegy;
                    rhs1 = grdx1 - begsegx;
                    rhs2 = grdy1 - begsegy;
                }
                else
                {
                    mat1 = x2 - endx;
                    mat3 = y2 - endy;
                    rhs1 = grdx1 - endx;
                    rhs2 = grdy1 - endy;
                }

                determ = mat1*mat4 - mat2*mat3;

                // re-check non-zero of determ
                if (nonzero(determ))
                {
                    s1 = (rhs1*mat4 - mat2*rhs2) / determ;
                    s2 = (mat1*rhs2 - rhs1*mat3) / determ;

                    if (! loutside)
                    {
                        intrsct_x = begsegx + mat1*s1;
                        intrsct_y = begsegy + mat3*s1;
                    }
                    else
                    {
                        intrsct_x = endx + mat1*s1;
                        intrsct_y = endy + mat3*s1;
                    }

                    // convert back to lat/lon coordinates
                    intrsct_lon = rns * atan2(intrsct_y, intrsct_x);
                    if (intrsct_lon < ZERO)
                        intrsct_lon = intrsct_lon + PI2;

                    if (ABS(intrsct_x) > 1e-10)
                        intrsct_lat = (PI4 - asin(rns * HALF * intrsct_x / cos(intrsct_lon))) * TWO;
                    else if (ABS(intrsct_y) > 1e-10)
                        intrsct_lat = (PI4 - asin(      HALF * intrsct_y / sin(intrsct_lon))) * TWO;
                    else
                        intrsct_lat = TWO * PI4;

                    // add offset in transformed space for next pass
                    if (s1 - eps / determ < ONE)
                    {
                        intrsct_x = intrsct_x - mat1 * (eps / determ);
                        intrsct_y = intrsct_y - mat3 * (eps / determ);
                    }
                    else
                    {
                        if (! loutside)
                        {
                            intrsct_x = endx;
                            intrsct_y = endy;
                            intrsct_lat = endlat;
                            intrsct_lon = endlon;
                        }
                        else
                        {
                            intrsct_x = begsegx;
                            intrsct_y = begsegy;
                            intrsct_lat = begseg[0];
                            intrsct_lon = begseg[1];
                        }
                    }
                }           /* determ nonzero */
            }               /* s1 s2 */
        }                   /* determ nonzero */
        // no intersection this side, move on to next side
    }                       /* intrsct loop */

    /*
     *  if segment manages to cross over pole, shift the beginning
     *  endpoint in order to avoid hitting pole directly
     *  (it is ok for endpoint to be pole point)
     */
    if (ABS(intrsct_x) < 1e-10 && ABS(intrsct_y) < 1e-10 &&
            nonzero(endx) && nonzero(endy))
    {
        if (avoid_pole_count > 2)
        {
            avoid_pole_count = 0;
            avoid_pole_offset = 10 * avoid_pole_offset;
        }

        cp = begsegx * (endy - begsegy) - begsegy * (endx - begsegx);
        intrsct_lat = begseg[0];
        if (cp * intrsct_lat > ZERO)
        {
            intrsct_lon = beglon    + avoid_pole_offset;
            begseg[1]   = begseg[1] + avoid_pole_offset;
        }
        else
        {
            intrsct_lon = beglon    - avoid_pole_offset;
            begseg[1]   = begseg[1] - avoid_pole_offset;
        }

        avoid_pole_count ++;
        luse_last = false;
    }
    else
    {
        avoid_pole_count    = 0;
        avoid_pole_offset   = tiny;
    }

    /*  if the segment crosses a pole threshold, reset the intersection
     *  to be the threhold latitude. only check if this was not a
     *  threshold segment since sometimes coordinate transform can end
     *  up on other side of threshold again
     */
    if (lthresh)
    {
        if (intrsct_lat < NORTH_THRESH && intrsct_lat > SOUTH_THRESH)
        {
            lthresh = false;
        }
    }
    else if (beglat > ZERO && intrsct_lat < NORTH_THRESH)
    {
        mat4 = endlat - begseg[0];
        mat3 = endlon - begseg[0];
        check_longitude(mat3, - PI, PI);
        intrsct_lat = NORTH_THRESH - tiny;
        s1          = (NORTH_THRESH - begseg[0]) / mat4;
        intrsct_lon = begseg[1] + s1 * mat3;
        luse_last   = false;
        lthresh     = true;
    }
    else if (beglat < ZERO && intrsct_lat > SOUTH_THRESH)
    {
        mat4 = endlat - begseg[0];
        mat3 = endlon - begseg[0];
        check_longitude(mat3, - PI, PI);
        intrsct_lat = SOUTH_THRESH + tiny;
        s1          = (SOUTH_THRESH - begseg[0]) / mat4;
        intrsct_lon = begseg[1] + s1 * mat3;
        luse_last   = false;
        lthresh     = true;
    }

    /* if reach end of segment, do not use x,y intersect on next entry */
    
    if (zero(intrsct_lat - endlat) || zero(intrsct_lon - endlon))
        luse_last = false;
}

