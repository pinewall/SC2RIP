#include "line_integral.h"
#include "constants.h"
#include "utils.h"
#include <math.h>

/** this routine computes the line integral of the flux function
  * that results in the interpolation weights. the line is defined
  * by the input lat/lon of the endpoints
  **/
void line_integral(double *weights, int num_weights, 
        double in_phi1, double in_phi2, double theta1, double theta2,
        double grid1_lat, double grid1_lon, double grid2_lat, double grid2_lon)
{
    // local variables
    double dphi, sinth1, sinth2, costh1, costh2, fac;
    double phi1, phi2, phidiff1, phidiff2, sinint;
    double f1, f2, fint;

    // weights for the general case based on a trapezoidal approx to the integrals
    sinth1 = sin(theta1);
    sinth2 = sin(theta2);
    costh1 = cos(theta1);
    costh2 = cos(theta2);

    dphi = in_phi1 - in_phi2;
    check_longitude(dphi);
    dphi = dphi * HALF;

    // the first weight is the area overlap integral. the second and fourth are second-order latitude gradient weights.
    weights[0] = dphi * (sinth1 + sinth2);
    weights[num_weights] = dphi * (sinth1 + sinth2);
    weights[1] = dphi * (costh1 + costh2 + (theta1*sinth1 + theta2*sinth2));
    weights[num_weights+1] = dphi * (costh1 + costh2 + (theta1*sinth1 + theta2*sinth2));

    // the third and fifth weights are for the second-order phi gradient component. must be careful of longitude range
    f1 = HALF * (costh1*sinth1 + theta1);
    f2 = HALF * (costh2*sinth2 + theta2);

    phi1 = in_phi1 - grid1_lon;
    check_longitude(phi1);

    phi2 = in_phi2 - grid1_lon;
    check_longitude(phi2);

    if ((phi2-phi1) < PI && (phi2-phi1) > - PI)     // normal case
    {
        weights[2] = dphi * (phi1*f1 + phi2*f2);
    }
    else            // cross longitude bi-value line case
    {
        if (phi1 > ZERO)
            fac = PI;
        else
            fac = - PI;
        fint = f1 + (f2-f1) * (fac-phi1) / ABS(dphi);
        weights[2] = HALF * phi1 * (phi1-fac) * f1 -
                     HALF * phi2 * (phi2+fac) * f2 +
                     HALF * fac * (phi1+phi2) * fint;
    }

    phi1 = in_phi1 - grid2_lon;
    check_longitude(phi1);

    phi2 = in_phi2 - grid2_lon;
    check_longitude(phi2);

    if ((phi2-phi1) < PI && (phi2-phi1) > - PI)     // normal case
    {
        weights[num_weights+2] = dphi * (phi1*f1 + phi2*f2);
    }
    else            // cross longitude bi-value line case
    {
        if (phi1 > ZERO)
            fac = PI;
        else
            fac = - PI;
        fint = f1 + (f2-f1) * (fac-phi1) / ABS(dphi);
        weights[num_weights+2] = HALF * phi1 * (phi1-fac) * f1 -
                     HALF * phi2 * (phi2+fac) * f2 +
                     HALF * fac * (phi1+phi2) * fint;
    }
}
