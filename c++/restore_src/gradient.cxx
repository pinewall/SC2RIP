#include "gradient.h"
#include "utils.h"
#include "constants.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

Gradient::Gradient(char * grad_type, int arr_size, double * arr_value, double * coord_lat, double * coord_lon)
{
    gradient_type = new char [256];
    strcpy(gradient_type, grad_type);
    if ((strcmp(gradient_type, "lat") == 0) &&
        (strcmp(gradient_type, "lon") == 0))
        gradient_num = 1;
    else if (strcmp(gradient_type, "latlon") == 0)
        gradient_num = 2;
    else
    {
        printf("Unknown gradient type\n");
        return;
    }
    
    array_size      = arr_size;
    array_value     = new double [array_size];
    coordinate_lat  = new double [array_size];
    coordinate_lon  = new double [array_size];
    memcpy(array_value, arr_value, sizeof(double) * array_size);
    memcpy(coordinate_lat, coord_lat, sizeof(double) * array_size);
    memcpy(coordinate_lon, coord_lon, sizeof(double) * array_size);

    delta_lat   = new double [100];
    delta_lon   = new double [100];
    delta_value = new double [100];
    memset(delta_lat,   0, sizeof(double) * 100);
    memset(delta_lon,   0, sizeof(double) * 100);
    memset(delta_value, 0, sizeof(double) * 100);

    gradient_result = new double [array_size * gradient_num];
}

void Gradient::calculate_gradient(int cell, double destination_value, double destination_lat, double destination_lon, int * index_of_neighbor, int num_of_neighbor)
{
    // check and save delta info
    for (int i = 0; i < num_of_neighbor; i++)
    {
        if (index_of_neighbor[i] >= array_size)
        {
            printf("Error: Index of neighbors exceeds Max value!\n");
            return;
        }

        delta_lat[i] = coordinate_lat[index_of_neighbor[i]] - destination_lat;
        delta_lon[i] = coordinate_lon[index_of_neighbor[i]] - destination_lon;
        check_longitude(delta_lon[i], - PI, PI);
        delta_value[i] = array_value[index_of_neighbor[i]] - destination_value;
    }

    double sum_square_delta_lat = 0.0;
    double sum_square_delta_lon = 0.0;
    double sum_delta_lat_delta_value = 0.0;
    double sum_delta_lon_delta_value = 0.0;
    double sum_delta_lat_delta_lon = 0.0;
    for (int i = 0; i < num_of_neighbor; i++)
    {
        sum_square_delta_lat += delta_lat[i] * delta_lat[i];
        sum_square_delta_lon += delta_lon[i] * delta_lon[i];
        sum_delta_lat_delta_value += delta_lat[i] * delta_value[i];
        sum_delta_lon_delta_value += delta_lon[i] * delta_value[i];
        sum_delta_lat_delta_lon += delta_lat[i] * delta_lon[i];
    }
    
    double gradient_lat = 0.0;
    double gradient_lon = 0.0;
    gradient_lat =  sum_square_delta_lat * sum_delta_lat_delta_value - 
                    sum_delta_lat_delta_lon * sum_delta_lon_delta_value;
    gradient_lon = -sum_delta_lat_delta_lon * sum_delta_lat_delta_value +                     sum_square_delta_lon * sum_delta_lon_delta_value;
    gradient_lat /= (sum_square_delta_lat * sum_square_delta_lon -
                     sum_delta_lat_delta_lon * sum_delta_lat_delta_lon);
    gradient_lon /= (sum_square_delta_lat * sum_square_delta_lon -
                     sum_delta_lat_delta_lon * sum_delta_lat_delta_lon);

    if (strcmp(gradient_type, "lat") == 0)
    {
        gradient_result[cell] = gradient_lat;
    }
    else if (strcmp(gradient_type, "lon") == 0)
    {
        gradient_result[cell] = gradient_lon;
    }
    else if (strcmp(gradient_type, "latlon") == 0)
    {
        gradient_result[cell * gradient_num]        = gradient_lat;
        gradient_result[cell * gradient_num + 1]    = gradient_lon;
    }
}
double Gradient::get_gradient_num()
{
    return gradient_num;
}

double * Gradient::get_gradient(int index)
{
    return (gradient_result + index * gradient_num);
}

Gradient::~Gradient()
{
    delete [] gradient_type;
    delete [] array_value;
    delete [] coordinate_lat;
    delete [] coordinate_lon;
    delete [] delta_lat;
    delete [] delta_lon;
    delete [] delta_value;
    delete [] gradient_result;
}
