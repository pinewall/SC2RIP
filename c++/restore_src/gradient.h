#ifndef _GRADIENT_H_
#define _GRADIENT_H_ 1

/* Matrix Operations */
void calculate_classical_company_matrix(double origin[3][3], double company[3][3]);
void calculate_matrix_multiply_vector(double * result, double matrix[3][3], double vector[3]);

/* Line integral Function using function pointer */
double line_integral_along_x(double x1, double y1, double x2, double y2, double (*f)(double, double));
double line_integral_along_y(double x1, double y1, double x2, double y2, double (*f)(double, double));

/* Bilinear Simulation ==> to calculate four parameters from z = a0 + a1*x + a2*y +  a3*x*y */
void calculate_bilinear_simulation(double * result, const double x1, const double y1, const double f1, const double x2, const double y2, const double f2, const double x3, const double y3, const double f3, const double x4, const double y4, const double f4);

/* Gradient Calculation Interface */
void calculate_derivative_latlon_least_square_method(double & destination_gradient_lat, double & destination_gradient_lon, 
        int index_of_destination,  const double * destination_lat, const double * destination_lon, 
        int num_of_sample_points, int * index_of_samples, const double * sample_lat, const double * sample_lon, 
        const double * field_value);

void calculate_derivative_lat_least_square_method(double & destination_gradient_lat,  
        int index_of_destination,  const double * destination_lat, const double * destination_lon, 
        int num_of_sample_points, int * index_of_samples, const double * sample_lat, const double * sample_lon, 
        const double * field_value);

void calculate_derivative_lon_least_square_method(double & destination_gradient_lon, 
        int index_of_destination,  const double * destination_lat, const double * destination_lon, 
        int num_of_sample_points, int * index_of_samples, const double * sample_lat, const double * sample_lon, 
        const double * field_value);

double calculate_derivative_x_least_square_method(double destination_coordination_x, double destination_coordination_y, double destination_function_value, int num_of_sample_points, double * sample_coordination_x, double * sample_coordination_y, double * sample_function_value);


double calculate_derivative_y_least_square_method(double destination_coordination_x, double destination_coordination_y, double destination_function_value, int num_of_sample_points, double * sample_coordination_x, double * sample_coordination_y, double * sample_function_value);

double calculate_derivative_x_area_average_method(double destination_coordination_x, double destination_coordination_y, double destination_function_value, int num_of_sample_points, double * sample_coordination_x, double * sample_coordination_y, double (*f)(double, double), int num_of_segments);

double calculate_derivative_y_area_average_method(double destination_coordination_x, double destination_coordination_y, double destination_function_value, int num_of_sample_points, double * sample_coordination_x, double * sample_coordination_y, double (*f)(double, double), int num_of_segments);

double calculate_derivative_x_bilinear_simulation_method(double destination_coordination_x, double destination_coordination_y, double destination_function_value, int num_of_sample_points, double * sample_coordination_x, double * sample_coordination_y, double * sample_function_value);

double calculate_derivative_y_bilinear_simulation_method(double destination_coordination_x, double destination_coordination_y, double destination_function_value, int num_of_sample_points, double * sample_coordination_x, double * sample_coordination_y, double * sample_function_value);


#endif
