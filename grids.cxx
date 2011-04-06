#include "grids.h"
#include "utils.h"
#include "namelist.h"
#include <iostream>
using namespace std;

/* extern variable define */
unsigned int grid1_size, grid2_size;        // size var define
unsigned int grid1_rank, grid2_rank;        // rank var define
int *grid1_dims, *grid2_dims;               // dims var define
unsigned int grid1_corners_max, grid2_corners_max;  // max num corners var define
unsigned int *grid1_corners, *grid2_corners;    //num corners var define
char *grid1_name, *grid2_name;              // grid name var define
char *grid1_units, *grid2_units;            // grid unit var define
bool *grid1_mask, *grid2_mask;              // grid mask var define
double *grid1_center_lat, *grid2_center_lat;// grid center lat var define
double *grid1_center_lon, *grid2_center_lon;// grid center lon var define
double *grid1_corner_lat, *grid2_corner_lat;// grid corner lat var define
double *grid1_corner_lon, *grid2_corner_lon;// grid corner lon var define
double *grid1_area, *grid2_area;            // grid area var define
double *grid1_area_in, *grid2_area_in;      // grid area in var define
double *grid1_frac, *grid2_frac;            // grid frac var define
double *grid1_bound_box, *grid2_bound_box;  // grid bound box var define
int *bin_addr1, *bin_addr2;                 // min|max index for each cell in bin
double *bin_lats, *bin_lons;                // min|max lat/lon for each srch bin
// grid variables init
void grid_init(char *grid1_file, char *grid2_file)
{
    grid_init_src(grid1_file);
    grid_init_dst(grid2_file);

    /* compute bounding boxes for restricting future grid searches */
    grid_cal_boundbox(grid1_bound_box, grid1_mask, grid1_size, grid1_center_lat, grid1_center_lon, grid1_corner_lat, grid1_corner_lon, grid1_corners, grid1_corners_max);    // bounding box for src grid
    grid_cal_boundbox(grid2_bound_box, grid2_mask, grid2_size, grid2_center_lat, grid2_center_lon, grid2_corner_lat, grid2_corner_lon, grid2_corners, grid2_corners_max);    // bounding box for dst grid

    /* assign address ranges to search bins in order to further restrict later searches */
    grid_srch_bin_init();   // init bin_lats|bin_lons
    grid_assign_srch_bin(grid1_bound_box, bin_addr1, grid1_size);   // assign srch bin address
    grid_assign_srch_bin(grid2_bound_box, bin_addr2, grid2_size);
}

// src grid varibles init
void grid_init_src(char *grid_src_file)
{
    cout << "SRC_GRID INIT" << endl;
    /* netCDF id for IO */
    int ncstat;             // netCDF status variable
    int nc_grid_id;         // netCDF grid file id
    int nc_gridsize_id;     // netCDF grid size dim id
    int nc_gridcorn_id;     // netCDF grid corner dim id
    int nc_gridrank_id;     // netCDF grid rank dim id
    int nc_gridarea_id;     // netCDF grid area dim id
    int nc_griddims_id;     // netCDF grid dimension size id
    int nc_grdimask_id;     // netCDF grid imask var id
    int nc_grdcrnrlat_id;   // netCDF grid corner/center lat/lon var id
    int nc_grdcrnrlon_id;
    int nc_grdcntrlat_id;
    int nc_grdcntrlon_id;

    /* local variables */
    int *imask;     // integer mask read from file

    /* open grid files and read size/name data */
    // open netCDF file
    ncstat = nc_open(grid_src_file, NC_NOWRITE, &nc_grid_id);
    ERR
    // read grid_size
    ncstat = nc_inq_dimid(nc_grid_id, "grid_size", &nc_gridsize_id);
    ERR
    ncstat = nc_inq_dimlen(nc_grid_id, nc_gridsize_id, &grid1_size);
    ERR
    // read grid_rank
    ncstat = nc_inq_dimid(nc_grid_id, "grid_rank", &nc_gridrank_id);
    ERR
    ncstat = nc_inq_dimlen(nc_grid_id, nc_gridrank_id, &grid1_rank);
    ERR
    // read grid_corners
    ncstat = nc_inq_dimid(nc_grid_id, "grid_corners", &nc_gridcorn_id);
    ERR
    ncstat = nc_inq_dimlen(nc_grid_id, nc_gridcorn_id, &grid1_corners_max);
    ERR
    // init coresponding grid corner number
    grid1_corners = new unsigned int[grid1_size];
    for (int i = 0; i < grid1_size; i++)
        grid1_corners[i] = grid1_corners_max;

    // allocate and read grid1_dims
    grid1_dims = new int[grid1_rank];
    ncstat = nc_inq_varid(nc_grid_id, "grid_dims", &nc_griddims_id);
    ERR
    ncstat = nc_get_var_int(nc_grid_id, nc_griddims_id, grid1_dims);
    ERR

    // read grid1_name
    grid1_name = new char[CHARLEN];       // grid name var define
    grid1_units = new char[CHARLEN];      // grid units var define
    ncstat = nc_get_att_text(nc_grid_id, NC_GLOBAL, "title", grid1_name);
    
    // allocate cell mask, area; read mask
    grid1_mask = new bool[grid1_size];
    grid1_area = new double[grid1_size];
    grid1_area_in = new double[grid1_size];
    grid1_frac = new double[grid1_size];
    imask = new int[grid1_size];
    ncstat = nc_inq_varid(nc_grid_id, "grid_imask", &nc_grdimask_id);
    ERR
    ncstat = nc_get_var_int(nc_grid_id, nc_grdimask_id, imask);
    for (int i = 0; i < grid1_size; i++)
    {
        grid1_mask[i] = (imask[i] == 1);
        grid1_area[i] = 0.0;
        grid1_frac[i] = 0.0;
    }

    // allocate and read center/corner lat/lon
    grid1_center_lat = new double[grid1_size];
    grid1_center_lon = new double[grid1_size];
    grid1_corner_lat = new double[grid1_size * grid1_corners_max];
    grid1_corner_lon = new double[grid1_size * grid1_corners_max];
    ncstat = nc_inq_varid(nc_grid_id, "grid_center_lat", &nc_grdcntrlat_id);
    ERR
    ncstat = nc_get_var_double(nc_grid_id, nc_grdcntrlat_id, grid1_center_lat);
    ERR
    ncstat = nc_inq_varid(nc_grid_id, "grid_center_lon", &nc_grdcntrlon_id);
    ERR
    ncstat = nc_get_var_double(nc_grid_id, nc_grdcntrlon_id, grid1_center_lon);
    ERR
    ncstat = nc_inq_varid(nc_grid_id, "grid_corner_lat", &nc_grdcrnrlat_id);
    ERR
    ncstat = nc_get_var_double(nc_grid_id, nc_grdcrnrlat_id, grid1_corner_lat);
    ERR
    ncstat = nc_inq_varid(nc_grid_id, "grid_corner_lon", &nc_grdcrnrlon_id);
    ERR
    ncstat = nc_get_var_double(nc_grid_id, nc_grdcrnrlon_id, grid1_corner_lon);
    ERR

    // read var units
    ncstat = nc_get_att_text(nc_grid_id, nc_grdcntrlat_id, "units", grid1_units);
    ERR

    // convert lat/lon units if required
    if (strcmp(grid1_units, "radians") == 0)
        cout << "No Convertion Required" << endl;// no convertion required
    else if (strcmp(grid1_units, "degrees") == 0)
    {
        cout << "Converting from Degrees to Radians" << endl;
        int index = 0;
        for (int i = 0; i < grid1_size; i++)
        {
            grid1_center_lat[i] *= deg2rad;
            grid1_center_lon[i] *= deg2rad;
            for (int j = 0; j < grid1_corners_max; j++)
            {
                grid1_corner_lat[index + j] *= deg2rad;
                grid1_corner_lon[index + j] *= deg2rad;
            }
            index += grid1_corners_max;
        }
    }
    else
        cout << "No Units Info, Set Units as Radians" << endl;
        // no units info; used as radians
    ncstat = nc_close(nc_grid_id);
    ERR

    /* make sure input latitude/longitude in given range */
    // latitude -PIH -- PIH
    grid_lat_range(grid1_center_lat, grid1_size);
    grid_lat_range(grid1_corner_lat, grid1_size * grid1_corners_max);

    // longitude ZERO -- PI2
    grid_lon_range(grid1_center_lon, grid1_size);
    grid_lon_range(grid1_corner_lon, grid1_size * grid1_corners_max);

    /* allocate bounding box array */
    grid1_bound_box = new double[grid1_size * BOUNDBOX_SIZE];
    for (int i = 0; i < grid1_size * BOUNDBOX_SIZE; i++)
        grid1_bound_box[i] = 0.0;
}

// dst grid varibles init
void grid_init_dst(char *grid_dst_file)
{
    cout << "DST_GRID INIT" << endl;
    /* netCDF id for IO */
    int ncstat;             // netCDF status variable
    int nc_grid_id;         // netCDF grid file id
    int nc_gridsize_id;     // netCDF grid size dim id
    int nc_gridcorn_id;     // netCDF grid corner dim id
    int nc_gridrank_id;     // netCDF grid rank dim id
    int nc_gridarea_id;     // netCDF grid area dim id
    int nc_griddims_id;     // netCDF grid dimension size id
    int nc_grdimask_id;     // netCDF grid imask var id
    int nc_grdcrnrlat_id;   // netCDF grid corner/center lat/lon var id
    int nc_grdcrnrlon_id;
    int nc_grdcntrlat_id;
    int nc_grdcntrlon_id;

    /* local variables */
    int *imask;     // integer mask read from file

    /* open grid files and read size/name data */
    // open netCDF file
    ncstat = nc_open(grid_dst_file, NC_NOWRITE, &nc_grid_id);
    ERR
    // read grid_size
    ncstat = nc_inq_dimid(nc_grid_id, "grid_size", &nc_gridsize_id);
    ERR
    ncstat = nc_inq_dimlen(nc_grid_id, nc_gridsize_id, &grid2_size);
    ERR
    // read grid_rank
    ncstat = nc_inq_dimid(nc_grid_id, "grid_rank", &nc_gridrank_id);
    ERR
    ncstat = nc_inq_dimlen(nc_grid_id, nc_gridrank_id, &grid2_rank);
    ERR
    // read grid_corners
    ncstat = nc_inq_dimid(nc_grid_id, "grid_corners", &nc_gridcorn_id);
    ERR
    ncstat = nc_inq_dimlen(nc_grid_id, nc_gridcorn_id, &grid2_corners_max);
    ERR
    // init coresponding grid corner number
    grid2_corners = new unsigned int[grid2_size];
    for (int i = 0; i < grid2_size; i++)
        grid2_corners[i] = grid2_corners_max;

    // allocate and read grid2_dims
    grid2_dims = new int[grid2_rank];
    ncstat = nc_inq_varid(nc_grid_id, "grid_dims", &nc_griddims_id);
    ERR
    ncstat = nc_get_var_int(nc_grid_id, nc_griddims_id, grid2_dims);
    ERR

    // read grid2_name
    grid2_name = new char[CHARLEN];       // grid name var define
    grid2_units = new char[CHARLEN];      // grid units var define
    ncstat = nc_get_att_text(nc_grid_id, NC_GLOBAL, "title", grid2_name);
    
    // allocate cell mask, area; read mask
    grid2_mask = new bool[grid2_size];
    grid2_area = new double[grid2_size];
    grid2_area_in = new double[grid2_size];
    grid2_frac = new double[grid2_size];
    imask = new int[grid2_size];
    ncstat = nc_inq_varid(nc_grid_id, "grid_imask", &nc_grdimask_id);
    ERR
    ncstat = nc_get_var_int(nc_grid_id, nc_grdimask_id, imask);
    for (int i = 0; i < grid2_size; i++)
    {
        grid2_mask[i] = (imask[i] == 2);
        grid2_area[i] = 0.0;
        grid2_frac[i] = 0.0;
    }

    // allocate and read center/corner lat/lon
    grid2_center_lat = new double[grid2_size];
    grid2_center_lon = new double[grid2_size];
    grid2_corner_lat = new double[grid2_size * grid2_corners_max];
    grid2_corner_lon = new double[grid2_size * grid2_corners_max];
    ncstat = nc_inq_varid(nc_grid_id, "grid_center_lat", &nc_grdcntrlat_id);
    ERR
    ncstat = nc_get_var_double(nc_grid_id, nc_grdcntrlat_id, grid2_center_lat);
    ERR
    ncstat = nc_inq_varid(nc_grid_id, "grid_center_lon", &nc_grdcntrlon_id);
    ERR
    ncstat = nc_get_var_double(nc_grid_id, nc_grdcntrlon_id, grid2_center_lon);
    ERR
    ncstat = nc_inq_varid(nc_grid_id, "grid_corner_lat", &nc_grdcrnrlat_id);
    ERR
    ncstat = nc_get_var_double(nc_grid_id, nc_grdcrnrlat_id, grid2_corner_lat);
    ERR
    ncstat = nc_inq_varid(nc_grid_id, "grid_corner_lon", &nc_grdcrnrlon_id);
    ERR
    ncstat = nc_get_var_double(nc_grid_id, nc_grdcrnrlon_id, grid2_corner_lon);
    ERR

    // read var units
    ncstat = nc_get_att_text(nc_grid_id, nc_grdcntrlat_id, "units", grid2_units);
    ERR

    // convert lat/lon units if required
    if (strcmp(grid2_units, "radians") == 0)
        cout << "No Convertion Required" << endl;// no convertion required
    else if (strcmp(grid2_units, "degrees") == 0)
    {
        int index = 0;
        for (int i = 0; i < grid2_size; i++)
        {
            grid2_corner_lat[i] *= deg2rad;
            grid2_corner_lon[i] *= deg2rad;
            for (int j = 0; j < grid2_corners_max; j++)
            {
                grid2_center_lat[index + j] *= deg2rad;
                grid2_center_lon[index + j] *= deg2rad;
            }
            index += grid2_corners_max;
        }
    }
    else
        cout << "No Units Info, Set Units as Radians" << endl;
        // no units info; used as radians
    ncstat = nc_close(nc_grid_id);
    ERR

    /* make sure input latitude/longitude in given range */
    // latitude -PIH -- PIH
    grid_lat_range(grid2_center_lat, grid2_size);
    grid_lat_range(grid2_corner_lat, grid2_size * grid2_corners_max);

    // longitude ZERO -- PI2
    grid_lon_range(grid2_center_lon, grid2_size);
    grid_lon_range(grid2_corner_lon, grid2_size * grid2_corners_max);

    /* allocate bounding box array */
    grid2_bound_box = new double[grid2_size * BOUNDBOX_SIZE];
    for (int i = 0; i < grid2_size * BOUNDBOX_SIZE; i++)
        grid2_bound_box[i] = 0.0;
}

// latitude/longitude range control
void grid_lat_range(double * lat, int len)
{
    for (int i = 0; i < len; i++)
    {
        if (lat[i] > PIH)
            lat[i] = PIH;
        if (lat[i] < -PIH)
            lat[i] = -PIH;
    }
}

void grid_lon_range(double * lon, int len)
{
    for (int i = 0; i < len; i++)
    {
        if (lon[i] > PI2)
            lon[i] -= PI2;
        if (lon[i] < ZERO)
            lon[i] += PI2;
    }
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
            bbdex += 2;
            index += grid_corners_max;
        }

        // consider two-value longitude around 0/PI2
        bbdex = 0;
        for (int i = 0; i < grid_size; i++)
        {
            if (boundbox[bbdex+3] - boundbox[bbdex+2] > PI)
            {
                boundbox[bbdex+2] = ZERO;
                boundbox[bbdex+3] = PI2;
            }
            bbdex += 6;
        }

        // try to check for cells that overlap poles
        bbdex = 0;
        for (int i = 0; i < grid_size; i++)
        {
            if (center_lat[i] > boundbox[bbdex+1])
                boundbox[bbdex+1] = PIH;
            if (center_lat[i] < boundbox[bbdex])
                boundbox[bbdex] = -PIH;
            bbdex += 6;
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
        double dlat = PI / num_srch_bins;       // srch bin interval
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
            bin_lats[latdex++] = i * dlat - PIH;
            bin_lats[latdex++] = (i+1) * dlat - PIH;
            bin_lons[londex++] = ZERO;
            bin_lons[londex++] = PI2;
            bin_addr1[adex1++] = grid1_size;
            bin_addr1[adex1++] = - ONE;
            bin_addr2[adex2++] = grid2_size;
            bin_addr2[adex2++] = - ONE;
        }
    }
    else if (strcmp(restrict_type, "latlon") == 0)
    {
        cout << "Using latlon bins to restrict search" << endl;
        double dlat = PI / num_srch_bins;
        double dlon = PI2 / num_srch_bins;
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
                bin_lats[latdex++] = j * dlat - PIH;
                bin_lats[latdex++] = (j+1) * dlat - PIH;
                bin_lons[londex++] = i * dlon;
                bin_lons[londex++] = (i+1) * dlon;
                bin_addr1[adex1++] = grid1_size;
                bin_addr1[adex1++] = - ONE;
                bin_addr2[adex2++] = grid2_size;
                bin_addr2[adex2++] = - ONE;
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
    int bbdex = 0;
    for (int cell = 0; cell < grid_size; cell++)
    {
        for (int n = 0; n < num_srch_bins_all; n++)
        {
            if (boundbox[bbdex] <= bin_lats[n*2] && boundbox[bbdex+1] >= bin_lats[n*2+1])
                if (boundbox[bbdex+2] <= bin_lons[n*2] && boundbox[bbdex+3] >= bin_lons[n*2+1])
                {
                    addr[n*2] = min(cell, addr[n*2]);
                    addr[n*2+1] = max(cell, addr[n*2+1]);
                }
        }
        bbdex += 6;
    }
}
