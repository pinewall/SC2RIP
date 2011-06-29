#include "gradient.h"
#include "utils.h"
#include "constants.h"

#include <math.h>
#include <stdio.h>

#define GRAD_TEST 0
/* Line integral Function using function pointer */
double line_integral_along_x(double x1, double y1, double x2, double y2, double (*f)(double, double))
{
    double result = 0.0;
    double tmp_point_x, tmp_point_y, tmp_function;
    int point, point_max = 1000;
    double dx = (x2 - x1) / point_max;
    double dy = (y2 - y1) / point_max;

    for (point = 0; point < point_max; point ++)
    {
        tmp_point_x = x1 + point * dx;
        tmp_point_y = y1 + point * dy;
        tmp_function = f(tmp_point_x, tmp_point_y) + f(tmp_point_x + dx, tmp_point_y + dy);
        tmp_function /= 2;
        result += tmp_function;
    }
    result *= dx;
    return result;
}

double line_integral_along_y(double x1, double y1, double x2, double y2, double (*f)(double, double))
{
    double result = 0.0;
    double tmp_point_x, tmp_point_y, tmp_function;
    int point, point_max = 1000;
    double dx = (x2 - x1) / point_max;
    double dy = (y2 - y1) / point_max;

    for (point = 0; point < point_max; point ++)
    {
        tmp_point_x = x1 + point * dx;
        tmp_point_y = y1 + point * dy;
        tmp_function = f(tmp_point_x, tmp_point_y) + f(tmp_point_x + dx, tmp_point_y + dy);
        tmp_function /= 2;
        result += tmp_function;
    }
    result *= dy;
    return result;
}

/* Bilinear Simulation ==> to calculate four parameters from z = a0 + a1*x + a2*y +  a3*x*y */
void calculate_bilinear_simulation(double * result, const double x1, const double y1, const double f1, const double x2, const double y2, const double f2, const double x3, const double y3, const double f3, const double x4, const double y4, const double f4)
{
    int row, col;
    // degenerate matrix problem By=g
    double B[3][3] = { 
                    x2 - x1,    y2 - y1,    x2*y2 - x1*y1,
                    x3 - x1,    y3 - y1,    x3*y3 - x1*y1,
                    x4 - x1,    y4 - y1,    x4*y4 - x1*y1,
    };

#if GRAD_TEST
    // test matrix B
    print_matrix(B);
    printf("\n");
#endif
    
    // Classical Company Matrix of B
    double B_classical_company[3][3] = {0.0};
    calculate_classical_company_matrix(B, B_classical_company);

#if GRAD_TEST
    // check reverse matrix
    int reverse_correct = test_reverse_matrix(B, B_classical_company);
    if (reverse_correct)
        printf("Reverse matrix ==> right!\n");
    else
        printf("Reverse matrix ==> wrong!\n");
#endif

    double delta_function[3] = {f2 - f1, f3 - f1, f4 - f1};
    for (row = 0; row < 3; row ++)
    {
        for (col = 0; col < 3; col ++)
        {
            result[row+1] += (B_classical_company[row][col] * delta_function[col]);
        }
    }

    result[0] = f1 - result[1]*x1 - result[2]*y1 - result[3]*x1*y1;
}

/* Matrix Related Calculation */
void calculate_classical_company_matrix(double origin[3][3], double company[3][3])
{
    int row, col;
    for (row = 0; row < 3; row ++)
    {
        for (col = 0; col < 3; col ++)
        {
            company[row][col] = origin[(col+1)%3][(row+1)%3] * origin[(col+2)%3][(row+2)%3]  -  origin[(col+2)%3][(row+1)%3] * origin[(col+1)%3][(row+2)%3];
        }
    }
    double determine_origin = 0.0;
    for (col = 0; col < 3; col ++)
    {
        determine_origin += origin[0][col] * company[col][0];
    }

    // store Reverse Matrix of B in B_classical_company
    for (row = 0; row < 3; row ++)
    {
        for (col = 0; col < 3; col ++)
        {
            company[row][col] /= determine_origin;
        }
    }
}

void calculate_matrix_multiply_vector(double * result, double matrix[3][3], double vector[3])
{
    int row, col;
    for (row = 0; row < 3; row ++)
    {
        result[row] = 0.0;
        for (col = 0; col < 3; col ++)
        {
            result[row] += matrix[row][col] * vector[col];
        }
    }
}


/* Gradient Calculation Interface */

/* Least Square Method */

// --for lat/lon-- //
void calculate_derivative_latlon_least_square_method(double & destination_gradient_lat, double & destination_gradient_lon, 
        int index_of_destination,  const double * destination_lat, const double * destination_lon, 
        int num_of_sample_points, int * index_of_samples, const double * sample_lat, const double * sample_lon, 
        const double * field_value)
{
    double sum_of_square_delta_lat = 0.0;
    double sum_of_square_delta_lon = 0.0;
    double sum_of_delta_lat_by_delta_lon = 0.0;
    double sum_of_delta_lat_by_delta_field = 0.0;
    double sum_of_delta_lon_by_delta_field = 0.0;

    double delta_lat_tmp, delta_lon_tmp, delta_field_tmp;
    for (int index = 0; index < num_of_sample_points; index ++)
    {
        delta_lat_tmp = sample_lat[index_of_samples[index]] - destination_lat[index_of_destination];
        delta_lon_tmp = sample_lon[index_of_samples[index]] - destination_lon[index_of_destination];
        delta_field_tmp = field_value[index_of_samples[index]] - field_value[index_of_destination];

        // check longitude before operations
        check_longitude(delta_lon_tmp, -PI, PI);

        sum_of_square_delta_lat += delta_lat_tmp * delta_lat_tmp;
        sum_of_square_delta_lon += delta_lon_tmp * delta_lon_tmp;
        sum_of_delta_lat_by_delta_lon += delta_lat_tmp * delta_lon_tmp;
        sum_of_delta_lat_by_delta_field += delta_lat_tmp * delta_field_tmp;
        sum_of_delta_lon_by_delta_field += delta_lon_tmp * delta_field_tmp;
    }
    
    // do matrix multifly vector
    destination_gradient_lat = sum_of_square_delta_lat * sum_of_delta_lat_by_delta_field - sum_of_delta_lat_by_delta_lon * sum_of_delta_lon_by_delta_field;
    destination_gradient_lon = - sum_of_delta_lat_by_delta_lon * sum_of_delta_lat_by_delta_field + sum_of_square_delta_lon * sum_of_delta_lon_by_delta_field;

    double determine_gradient = sum_of_square_delta_lat * sum_of_delta_lon_by_delta_field - sum_of_delta_lat_by_delta_lon * sum_of_delta_lat_by_delta_lon;
    if (zero(determine_gradient))
    {
        printf("Singular Matrix, No Gradient needed!\n");
        destination_gradient_lat = 0.0;
        destination_gradient_lon = 0.0;
    }
    else
    {
        destination_gradient_lat /= determine_gradient;
        destination_gradient_lon /= determine_gradient;
    }

    // for longitude gradient, we need to multiply abs(1/sin(theta))
    double sin_dst_lat = sin(destination_lat[index_of_destination]);
    destination_gradient_lon /= sin_dst_lat;
    if (sin_dst_lat < 0)
        destination_gradient_lon = - destination_gradient_lon;
}

void calculate_derivative_lat_least_square_method(double & destination_gradient_lat,  
        int index_of_destination,  const double * destination_lat, const double * destination_lon, 
        int num_of_sample_points, int * index_of_samples, const double * sample_lat, const double * sample_lon, 
        const double * field_value)
{
    double used_for_padding;
    calculate_derivative_latlon_least_square_method(destination_gradient_lat, used_for_padding, 
            index_of_destination, destination_lat, destination_lon,
            num_of_sample_points, index_of_samples, sample_lat, sample_lon,
            field_value);
}

void calculate_derivative_lon_least_square_method(double & destination_gradient_lon, 
        int index_of_destination,  const double * destination_lat, const double * destination_lon, 
        int num_of_sample_points, int * index_of_samples, const double * sample_lat, const double * sample_lon, 
        const double * field_value)
{
    double used_for_padding;
    calculate_derivative_latlon_least_square_method(used_for_padding, destination_gradient_lon, 
            index_of_destination, destination_lat, destination_lon,
            num_of_sample_points, index_of_samples, sample_lat, sample_lon,
            field_value);
}

// --fox x/y-- //
double calculate_derivative_x_least_square_method(double destination_coordination_x, double destination_coordination_y, double destination_function_value, int num_of_sample_points, double * sample_coordination_x, double * sample_coordination_y, double * sample_function_value)
{
    int i;
    double square_xx = 0.0;
    double square_yy = 0.0;
    double square_xy = 0.0;
    double square_fx = 0.0;
    double square_fy = 0.0;

    for (i = 0; i < num_of_sample_points; i++)
    {
        square_xx += ((sample_coordination_x[i] - destination_coordination_x) * (sample_coordination_x[i] - destination_coordination_x));
        square_yy += ((sample_coordination_y[i] - destination_coordination_y) * (sample_coordination_y[i] - destination_coordination_y));
        square_xy += ((sample_coordination_x[i] - destination_coordination_x) * (sample_coordination_y[i] - destination_coordination_y));
        square_fx += ((sample_function_value[i] - destination_function_value) * (sample_coordination_x[i] - destination_coordination_x));
        square_fy += ((sample_function_value[i] - destination_function_value) * (sample_coordination_y[i] - destination_coordination_y));
    }

    // Now just to solve linear equations (xx xy;xy yy) * sol = (fx;fy)
    double determine_lsm = square_xx * square_yy - square_xy * square_xy;
    double derivative_x_lsm, derivative_y_lsm;
    derivative_x_lsm =  square_yy * square_fx - square_xy * square_fy;
    derivative_y_lsm = -square_xy * square_fx + square_xx * square_fy;
    derivative_x_lsm /= determine_lsm;
    derivative_y_lsm /= determine_lsm;

    return derivative_x_lsm;
}

double calculate_derivative_y_least_square_method(double destination_coordination_x, double destination_coordination_y, double destination_function_value, int num_of_sample_points, double * sample_coordination_x, double * sample_coordination_y, double * sample_function_value)
{
    int i;
    double square_xx = 0.0;
    double square_yy = 0.0;
    double square_xy = 0.0;
    double square_fx = 0.0;
    double square_fy = 0.0;

    for (i = 0; i < num_of_sample_points; i++)
    {
        square_xx += ((sample_coordination_x[i] - destination_coordination_x) * (sample_coordination_x[i] - destination_coordination_x));
        square_yy += ((sample_coordination_y[i] - destination_coordination_y) * (sample_coordination_y[i] - destination_coordination_y));
        square_xy += ((sample_coordination_x[i] - destination_coordination_x) * (sample_coordination_y[i] - destination_coordination_y));
        square_fx += ((sample_function_value[i] - destination_function_value) * (sample_coordination_x[i] - destination_coordination_x));
        square_fy += ((sample_function_value[i] - destination_function_value) * (sample_coordination_y[i] - destination_coordination_y));
    }

    // Now just to solve linear equations (xx xy;xy yy) * sol = (fx;fy)
    double determine_lsm = square_xx * square_yy - square_xy * square_xy;
    double derivative_x_lsm, derivative_y_lsm;
    derivative_x_lsm =  square_yy * square_fx - square_xy * square_fy;
    derivative_y_lsm = -square_xy * square_fx + square_xx * square_fy;
    derivative_x_lsm /= determine_lsm;
    derivative_y_lsm /= determine_lsm;

    return derivative_y_lsm;
}

double calculate_derivative_x_area_average_method(double destination_coordination_x, double destination_coordination_y, double destination_function_value, int num_of_sample_points, double * sample_coordination_x, double * sample_coordination_y, double (*f)(double, double), int num_of_segments)
{
    int i, next_i;
    double dx = 0.3;
    double dy = 0.3;
    double derivative_x_aam = 0.0;
    double line_integral_segment_along_y = 0.0;
    for (i = 0; i < num_of_sample_points; i++)
    {
        next_i = (i + 1) % num_of_sample_points;        // next point index in counter-clockwise loop
        line_integral_segment_along_y = line_integral_along_y(sample_coordination_x[i], sample_coordination_y[i], sample_coordination_x[next_i], sample_coordination_y[next_i], f);
        derivative_x_aam += line_integral_segment_along_y;
    }
    double area = dx * dy * 4;      // area of considered district
    derivative_x_aam /= area;

    return derivative_x_aam;
}

double calculate_derivative_y_area_average_method(double destination_coordination_x, double destination_coordination_y, double destination_function_value, int num_of_sample_points, double * sample_coordination_x, double * sample_coordination_y, double (*f)(double, double), int num_of_segments)
{
    int i, next_i;
    double dx = 0.3;
    double dy = 0.3;
    double derivative_y_aam = 0.0;
    double line_integral_segment_along_x = 0.0;
    for (i = 0; i < num_of_sample_points; i++)
    {
        next_i = (i + 1) % num_of_sample_points;        // next point index in counter-clockwise loop
        line_integral_segment_along_x = line_integral_along_x(sample_coordination_x[i], sample_coordination_y[i], sample_coordination_x[next_i], sample_coordination_y[next_i], f);
        derivative_y_aam += line_integral_segment_along_x;
    }
    double area = dx * dy * 4;      // area of considered district
    derivative_y_aam /= area;

    derivative_y_aam *= -1;         // counter-clockwise ==> clockwise, just for integral for derivative_y

    return derivative_y_aam;
}

double calculate_derivative_x_bilinear_simulation_method(double destination_coordination_x, double destination_coordination_y, double destination_function_value, int num_of_sample_points, double * sample_coordination_x, double * sample_coordination_y, double * sample_function_value);

double calculate_derivative_y_bilinear_simulation_method(double destination_coordination_x, double destination_coordination_y, double destination_function_value, int num_of_sample_points, double * sample_coordination_x, double * sample_coordination_y, double * sample_function_value);
