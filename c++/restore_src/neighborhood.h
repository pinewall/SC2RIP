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

    public:
        Neighborhood(char * grid_name, int grid_size, int num_of_neighbors);

        int     get_num_of_cells();
        int     get_num_of_neighbor_cells(int cell);
        int *   get_index_of_neighbor_cells(int cell);   
};

#endif
