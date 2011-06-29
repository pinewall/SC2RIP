#ifndef _NEIGHBORHOOD_H_
#define _NEIGHBORHOOD_H_ 1

/* neighborhood of each cell in same grid
 *  --Only for Gradient calculating now  
 */
class Neighborhood
{
    private:
        int     num_of_cells;               // number of cells in grid
        int     num_of_different_lats;      // number of different lats
        int     num_of_different_lons;      // number of different lons
        int *   num_of_neighbor_cells;      // number of neighbor cells considered for each cell

        int **  index_of_neighbor_cells;    // cell index array of neighbor cells for each cell
        void    find_latlon_neighborhood(); // toy case, just for testing T42 --> POP43 
        void    find_exact_num_neighborhood();

        // define some const pointer to use extern variables without modifying
        const double        *   grid_center_lat;
        const double        *   grid_center_lon;
        const double        *   grid_centroid_lat;
        const double        *   grid_centroid_lon;
        const double        *   grid_corner_lat;
        const double        *   grid_corner_lon;
        const bool          *   grid_mask;
        const unsigned int  *   grid_corners;

    public:
        Neighborhood(char * grid_name, int num_of_neighbors);
        void calculate_gradient_latlon(double * result_lat, double * result_lon, const double * field_value, int grid_size);
        void calculate_gradient_lat(double * result_lat, const double * field_value, int grid_size);
        void calculate_gradient_lon(double * result_lon, const double * field_value, int grid_size);
};

#endif
