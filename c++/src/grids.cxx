#include "grids.h"
#include "utils.h"
#include "namelist.h"
#include <iostream>
using namespace std;

/* extern variable define */
extern Grid * grid_src;
extern Grid * grid_dst;

// grid variables init
void grid_init(char *grid_file, char *grid2_file)
{
    grid_init_src(grid_file);
    grid_init_dst(grid2_file);

    /* compute bounding boxes for restricting future grid searches */
    grid_cal_boundbox(grid_bound_box, grid_mask, grid_size, grid_center_lat, grid_center_lon, grid_corner_lat, grid_corner_lon, grid_corners, grid_corners_max);    // bounding box for src grid
    grid_cal_boundbox(grid2_bound_box, grid2_mask, grid2_size, grid2_center_lat, grid2_center_lon, grid2_corner_lat, grid2_corner_lon, grid2_corners, grid2_corners_max);    // bounding box for dst grid

    /* assign address ranges to search bins in order to further restrict later searches */
    grid_srch_bin_init();   // init bin_lats|bin_lons

    grid_assign_srch_bin(grid_bound_box, bin_addr1, grid_size);   // assign srch bin address
    grid_assign_srch_bin(grid2_bound_box, bin_addr2, grid2_size);
    
}

// src grid varibles init
void grid_init(char *grid_file, Grid * grid)
{
    /* local variables */
    unsigned int grid_size;        // size var define
    unsigned int grid_rank;        // rank var define
    int * grid_dims;               // dims var define
    unsigned int grid_corners_max;  // max num corners var define
    unsigned int * grid_corners;    //num corners var define
    char * grid_name;              // grid name var define
    char * grid_units;            // grid unit var define
    bool * grid_mask;              // grid mask var define
    double * grid_center_lat;   // grid center lat var define
    double * grid_center_lon;   // grid center lon var define
    double * grid_corner_lat;   // grid corner lat var define
    double * grid_corner_lon;   // grid corner lon var define
    double *grid_area;            // grid area var define
    double *grid_area_in;      // grid area in var define
    double *grid_frac;            // grid frac var define
    int *imask;     // integer mask read from file

    /* open grid files and read size/name data */
    // open netCDF file
    NCUtils * nc = new NCUtils(grid_file);

    // read grid_size
    nc->fetch_dim(grid_size, "grid_size");

    // read grid_rank
    nc->fetch_dim(grid_rank, "grid_rank");

    // read grid_corners
    nc->fetch_dim(grid_corners_max, "grid_corners");

    // init coresponding grid corner number
    grid_corners = new unsigned int[grid_size];
    for (int i = 0; i < grid_size; i++)
        grid_corners[i] = grid_corners_max;

    // allocate and read grid_dims
    grid_dims = new int[grid_rank];
    nc->fetch_var_int(grid_dims, "grid_dims");

    // read grid_name
    grid_name = new char[Constant::CHARLEN];       // grid name var define
    grid_units = new char[Constant::CHARLEN];      // grid units var definea
    nc->fetch_att_text(grid_name, 0, "title");
    
    // allocate cell mask, area; read mask
    grid_mask = new bool[grid_size];
    grid_area = new double[grid_size];
    grid_area_in = new double[grid_size];
    grid_frac = new double[grid_size];
    imask = new int[grid_size];

    nc->fetch_var_int(imask, "grid_imask");
    for (int i = 0; i < grid_size; i++)
    {
        grid_mask[i] = (imask[i] == 1);
        grid_area[i] = 0.0;
        grid_frac[i] = 0.0;
    }

    // allocate and read center/corner lat/lon
    grid_center_lat = new double[grid_size];
    grid_center_lon = new double[grid_size];
    grid_corner_lat = new double[grid_size * grid_corners_max];
    grid_corner_lon = new double[grid_size * grid_corners_max];

    nc->fetch_var_double(grid_center_lat, "grid_center_lat");
    nc->fetch_var_double(grid_center_lon, "grid_center_lon");
    nc->fetch_var_double(grid_corner_lat, "grid_corner_lat");
    nc->fetch_var_double(grid_corner_lon, "grid_corner_lon");

    // read var units
    nc->fetch_att_text(grid_units, "grid_center_lat", "units");

    // convert lat/lon units if required
    if (strcmp(grid_units, "radians") == 0)
        cout << "No Convertion Required" << endl;// no convertion required
    else if (strcmp(grid_units, "degrees") == 0)
    {
        cout << "Converting from Degrees to Radians" << endl;
        int index = 0;
        for (int i = 0; i < grid_size; i++)
        {
            grid_center_lat[i] *= Constant::deg2rad;
            grid_center_lon[i] *= Constant::deg2rad;
            for (int j = 0; j < grid_corners_max; j++)
            {
                grid_corner_lat[index + j] *= Constant::deg2rad;
                grid_corner_lon[index + j] *= Constant::deg2rad;
            }
            index += grid_corners_max;
        }
    }
    else
        cout << "No Units Info, Set Units as Radians" << endl;
        // no units info; used as radians
    delete nc;

    /* make sure input latitude/longitude in given range */
    // latitude -Constant::PIH -- Constant::PIH
    check_latitude_all(grid_center_lat, grid_size);
    check_latitude_all(grid_corner_lat, grid_size * grid_corners_max);

    // longitude Constant::ZERO -- Constant::PI2
    check_longitude_all(grid_center_lon, grid_size, Constant::ZERO, Constant::PI2);
    check_longitude_all(grid_corner_lon, grid_size * grid_corners_max
                                                    , Constant::ZERO, Constant::PI2);

    // initialize Grid
    grid = new Grid(grid_rank, grid_size, grid_dims, 
            grid_corners, grid_corners_max, 
            grid_center_lat, grid_center_lon,
            grid_corner_lat, grid_corner_lon,
            grid_masks);
}


// bounding box calculation for grid
void grid_cal_boundbox(double *boundbox, bool *grid_mask, int grid_size, double *center_lat, double *center_lon, double *corner_lat, double *corner_lon, unsigned int *grid_corners, unsigned int grid_corners_max)
{
    double dlat, dlon;  // lat/lon intevals for srch bins
    double tmp_lats[6], tmp_lons[6];    // temps for computing bounding boxes
    if (! luse_grid_centers)
    {
#if _3D_BOUND_BOX_
        // 3D bounding box
        double XYZ[3][grid_corners_max]; 
        int index = 0;
        int bbdex = 0;
        for (int i = 0; i < grid_size; i++)
        {
            for (int j = 0; j < grid_corners[i]; j++)
            {
                XYZ[0][j] = cos(corner_lat[index + j]) * sin(corner_lon[index + j]);
                XYZ[1][j] = cos(corner_lat[index + j]) * sin(corner_lon[index + j]);
                XYZ[2][j] = sin(corner_lat[index + j]);
            }
            index += grid_corners_max;
            boundbox[bbdex++] = minval(XYZ[0], grid_corners[i]);
            boundbox[bbdex++] = maxval(XYZ[0], grid_corners[i]);
            boundbox[bbdex++] = minval(XYZ[1], grid_corners[i]);
            boundbox[bbdex++] = maxval(XYZ[1], grid_corners[i]);
            boundbox[bbdex++] = minval(XYZ[2], grid_corners[i]);
            boundbox[bbdex++] = maxval(XYZ[2], grid_corners[i]);
        }
#endif

#if _2D_BOUND_BOX_
        // 2D bounding box
        int bbdex = 0;  // index for boundbox
        int index = 0;  // index for center/corner
        double * lat_head, * lon_head;  // head for corners for each cell
        int num_corner;
        for (int i = 0; i < grid_size; i++)
        {
            lat_head = &(corner_lat[index]);
            lon_head = &(corner_lon[index]);
            num_corner = grid_corners[i];
            boundbox[bbdex++] = minval(lat_head, num_corner);
            boundbox[bbdex++] = maxval(lat_head, num_corner);
            boundbox[bbdex++] = minval(lon_head, num_corner);
            boundbox[bbdex++] = maxval(lon_head, num_corner);
            index += grid_corners_max;
        }

        // consider two-value longitude around 0/Constant::PI2
        /*  using le() and ge() to deal two-value problem*/
        bbdex = 0;
        for (int i = 0; i < grid_size; i++)
        {
            if (boundbox[bbdex+3] - boundbox[bbdex+2] > Constant::PI)
            {
                boundbox[bbdex+2] = Constant::ZERO;
                boundbox[bbdex+3] = Constant::PI2;
            }
            bbdex += BOUNDBOX_SIZE;
        }
        /**/
        // try to check for cells that overlap poles
        bbdex = 0;
        for (int i = 0; i < grid_size; i++)
        {
            if ((center_lat[i]) > (boundbox[bbdex+1]))
            {
                printf("overlap @ %d: %3.6f > %3.6f\n", i, center_lat[i], boundbox[bbdex+1]);
                boundbox[bbdex+1] = Constant::PIH;
            }
            else if ((center_lat[i]) < (boundbox[bbdex]))
            {
                printf("overlap @ %d: %3.6f < %3.6f\n", i, center_lat[i], boundbox[bbdex]);
                boundbox[bbdex] = -Constant::PIH;
            }
            bbdex += BOUNDBOX_SIZE;
        }
#endif
    }

}

// init search bins
void grid_srch_bin_init()
{
    if (strcmp(restrict_type, "latitude") == 0)
    {
        cout << "Using latitude bins to restrict search" << endl;
        double dlat = Constant::PI / num_srch_bins;       // srch bin interval
        bin_lats = new double [2 * num_srch_bins];
        bin_lons = new double [2 * num_srch_bins];
        bin_addr1 = new int [2 * num_srch_bins];
        bin_addr2 = new int [2 * num_srch_bins];

        int latdex = 0;
        int londex = 0;
        int adex1 = 0;
        int adex2 = 0;
        for (int i = 0; i < num_srch_bins; i++)
        {
            bin_lats[latdex++] = i * dlat - Constant::PIH;
            bin_lats[latdex++] = (i+1) * dlat - Constant::PIH;
            bin_lons[londex++] = Constant::ZERO;
            bin_lons[londex++] = Constant::PI2;
            bin_addr1[adex1++] = grid_size;
            bin_addr1[adex1++] = -1;
            bin_addr2[adex2++] = grid2_size;
            bin_addr2[adex2++] = -1;
        }
    }
    else if (strcmp(restrict_type, "latlon") == 0)
    {
        cout << "Using latlon bins to restrict search" << endl;
        double dlat = Constant::PI / num_srch_bins;
        double dlon = Constant::PI2 / num_srch_bins;
        int num_srch_bins_2 = num_srch_bins * num_srch_bins;
        bin_lats = new double [2 * num_srch_bins_2];
        bin_lons = new double [2 * num_srch_bins_2];
        bin_addr1 = new int [2 * num_srch_bins_2];
        bin_addr2 = new int [2 * num_srch_bins_2];
        int latdex = 0;
        int londex = 0;
        int adex1 = 0;
        int adex2 = 0;
        for (int i = 0; i < num_srch_bins; i++)
        {
            for (int j = 0; j < num_srch_bins; j++)
            {
                bin_lats[latdex++] = j * dlat - Constant::PIH;
                bin_lats[latdex++] = (j+1) * dlat - Constant::PIH;
                bin_lons[londex++] = i * dlon;
                bin_lons[londex++] = (i+1) * dlon;
                bin_addr1[adex1++] = grid_size;
                bin_addr1[adex1++] = -1;
                bin_addr2[adex2++] = grid2_size;
                bin_addr2[adex2++] = -1;
            }
        }
    }
    else
    {
        cout << "Unkown search restiction method" << endl;
        exit(1);
    }
}

// assign search bin address, that is min|max index of cell
void grid_assign_srch_bin(double *boundbox, int *addr, int grid_size)
{
    int num_srch_bins_all = 0;
    if (strcmp(restrict_type, "latitude") == 0)
        num_srch_bins_all = num_srch_bins;
    else if (strcmp(restrict_type, "latlon") == 0)
        num_srch_bins_all = num_srch_bins * num_srch_bins;
    else
    {
        cout << "Unknown search restriction method" << endl;
        exit(1);
    }
    if (streqls(restrict_type, "latitude"))
    {
        int bbdex = 0;
        for (int cell = 0; cell < grid_size; cell ++)
        {
            for (int n = 0; n < num_srch_bins_all; n++)
            {
                if (boundbox[bbdex+1] >= bin_lats[n*2] &&
                    boundbox[bbdex] <= bin_lats[n*2+1])
                {
                    addr[n*2] = MIN(cell, addr[n*2]);
                    addr[n*2+1] = MAX(cell, addr[n*2+1]);
                }
            }
            bbdex += BOUNDBOX_SIZE;
        }
    }
    else if (streqls(restrict_type, "latlon"))
    {
        int bbdex = 0;
        for (int cell = 0; cell < grid_size; cell++)
        {
            for (int n = 0; n < num_srch_bins_all; n++)
            {
                if (boundbox[bbdex+1] >= bin_lats[n*2] && boundbox[bbdex] <= bin_lats[n*2+1])
                    if (le(boundbox[bbdex+2], bin_lons[n*2+1]) && 
                        ge(boundbox[bbdex+3], bin_lons[n*2]))
                    {
                        addr[n*2] = MIN(cell, addr[n*2]);
                        addr[n*2+1] = MAX(cell, addr[n*2+1]);
                    }
            }
            bbdex += BOUNDBOX_SIZE;
        }
    }
}

