#ifndef _LINE_INTEGRAL_
#define _LINE_INTEGRAL_ 1

/** this routine computes the line integral of the flux function
  * that results in the interpolation weights. the line is defined
  * by the input lat/lon of the endpoints
  **/
void line_integral(double *weights, int num_weights, 
        double in_phi1, double in_phi2, double theta1, double theta2,
        double grid1_lat, double grid1_lon, double grid2_lat, double grid2_lon);

#endif
