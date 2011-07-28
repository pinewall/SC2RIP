#ifndef _GRADIENT_H_
#define _GRADIENT_H_ 1

class Gradient
{
    private:
        char    *   gradient_type;
        int         gradient_num;
        int         array_size;
        double  *   array_value;
        double  *   coordinate_lat;
        double  *   coordinate_lon;

        double  *   delta_lat;
        double  *   delta_lon;
        double  *   delta_value;

        double  *   gradient_result;
    public:
                    Gradient(char * grad_type, int arr_size, double * arr_value, double * coord_lat, double * coord_lon);
                    ~Gradient();
        void        calculate_gradient(int cell, double destination_value, double destination_lat, double destination_lon, int * index_of_neighbor, int num_of_neighbor);
        double  *   get_gradient(int cell);
        double      get_gradient_num();
};
#endif
