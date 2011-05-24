#include "great_circle_base.h"
#include "utils.h"
#include <math.h>
#include <stdio.h>


/** For Basic   class   <GCBase>    **/

GCBase::GCBase()
{
}

GCBase::~GCBase()
{
}

void GCBase::log()
{
    printf("log from class GCBase\n");
}

/*----------------------------------------*/


/** For Public  class   <Point>     **/

Point::Point(double latitude, double longitude)
{
    lat = latitude;
    lon = longitude;
    x = cos(lat) * cos(lon) * Point::Radius;
    y = cos(lat) * sin(lon) * Point::Radius;
    z = sin(lat)            * Point::Radius;
}

Point::Point(double x_coord, double y_coord, double z_coord)
{
    double r = sqrt(x_coord*x_coord + y_coord*y_coord + z_coord*z_coord);
    x = x_coord * Point::Radius / r;
    y = y_coord * Point::Radius / r;
    z = z_coord * Point::Radius / r;
    lat = asin(z/Point::Radius);
    double coslon = y / Point::Radius / cos(lat);
    double sinlon = x / Point::Radius / cos(lat);
    if (coslon > 0)
    {
        if (sinlon > 0)
            lon = asin(sinlon);
        else
            lon = asin(sinlon) + Constant::PI2;
    }
    else
    {
        if (sinlon > 0)
            lon = asin(sinlon) + Constant::PI;
        else
            lon = asin(sinlon) + Constant::PI * 3 / 2;
    }
}

Point::~Point()
{

}

Point * Point::cross_product(Point * right)
{
    double cp_x = y * right->z - right->y * z;
    double cp_y = z * right->x - right->z * x;
    double cp_z = x * right->y - right->x * y;
    return new Point(cp_x, cp_y, cp_z);
}

void Point::log()
{
    printf("log from class Point\n");
}

/*----------------------------------------*/


/** For Public  class   <Edge>      **/

Edge::Edge(double beglat, double beglon, double endlat, double endlon)
{
    start   = new Point (beglat, beglon);
    end     = new Point (endlat, endlon);
    norm_vec= start->cross_product(end);
}

Edge::~Edge()
{
}

void Edge::log()
{
    printf("log from class Edge\n");
}

/*----------------------------------------*/


/** For Public  class   <Cell>      **/

Cell::Cell(int cell_size, double cntr_lat, double cntr_lon, double * crnr_lat, double * crnr_lon)
{
    int next_corner;
    size = cell_size;
    corners = new Point * [cell_size];
    center = new Point(cntr_lat, cntr_lon);
    edges = new Edge * [cell_size];
    for (int corner = 0; corner < cell_size; corner ++)
    {
        next_corner = (corner + 1) % cell_size;
        corners[corner] = new Point(crnr_lat[corner], crnr_lon[corner]);
        edges[corner] = new Edge(crnr_lat[corner], crnr_lon[corner], crnr_lat[next_corner], crnr_lon[next_corner]);
    }
    area = 0.0;
    frac = 0.0;

    // calc bounding box
    calc_bbox();
}

Cell::~Cell()
{
    delete [] corners;
    delete center;
    delete [] edges;
    delete [] bbox;
    delete [] conns;
}

void Cell::calc_bbox()
{
    // initalize bounding box for this cell
    bbox = new double [4];

    // simply use plane method
    double lat_min = - Constant::PI2;
    double lat_max =   Constant::PI2;
    double lon_min =   corners[0]->lon;
    double lon_max =   corners[0]->lon;

    for (int corn = 1; corn < size; corn ++)
    {
        if (corners[corn]->lat < lat_min)
            lat_min = corners[corn]->lat;
        if (corners[corn]->lat > lat_max)
            lat_max = corners[corn]->lat;
        if (le(corners[corn]->lon, lon_min))
            lon_min = corners[corn]->lon;
        if (ge(corners[corn]->lon, lon_max))
            lon_max = corners[corn]->lon;
    }
    // judge whether this cell overlap pole
    int n, next_n;
    int same_direct_count = 1;  // count number of same direction from cp
    int standard, direction;  // set first cp result as standard direction

    Point * cp = GCUtils::cross_product(corners[1], corners[0]);
    standard = (cp->z > 0) ? 1 : -1;

    for (n = 1; n < size; n ++)
    {
        next_n = (n + 1) % size;
        delete cp;
        cp = GCUtils::cross_product(corners[next_n], corners[n]);
        direction = (cp->z > 0) ? 1 : -1;
        //direction *= standard;
        if (direction == standard)
            same_direct_count ++;
    }
    delete cp;

    // normal bounding box
    bbox[LAT_MIN] = lat_min;
    bbox[LAT_MAX] = lat_max;
    bbox[LON_MIN] = lon_min;
    bbox[LON_MAX] = lon_max;

    // Points in a circle, thus overlap case happens
    if (same_direct_count == size)
    {
        if (center->z > 0)
            bbox[LAT_MAX] =   Constant::PIH;
        else
            bbox[LAT_MIN] = - Constant::PIH;
        bbox[LON_MIN] = LON_MIN_VALUE;
        bbox[LON_MAX] = LON_MAX_VALUE;
    }
}

void Cell::log()
{
    printf("log from class Cell\n");
}

/*----------------------------------------*/


/** For Public  class   <Grid>      **/

Grid::Grid()
{
    rank = 0;       // zero integer
    size = 0;       // zero integer
    dims = 0;       // null integer pointer
    units = 0;      // zero integer
    cells = 0;      // null double pointer
    areas_in = 0;   // null double pointer
    fracs_in = 0;   // null double pointer
    imask = 0;      // null integer pointer
    bin_lats = 0;   // null double pointer
    bin_lons = 0;   // null double pointer
    bin_addr = 0;   // null integer pointer
    conns = 0;      // null integer pointer
}

Grid::Grid(int grid_rank, int grid_size, int * grid_dims, 
        int * grid_corners, int grid_corners_max,  
        double * grid_center_lat, double * grid_center_lon, 
        double * grid_corner_lat, double * grid_corner_lon, 
        int * cell_mask)
{
    double * crnr_lat = grid_corner_lat;
    double * crnr_lon = grid_corner_lon;
    
    // rank and dims
    rank = grid_rank;
    dims = new int [grid_rank];
    for (int d = 0; d < grid_rank; d ++)
        dims[d] = grid_dims[d];

    // size and related
    size = grid_size;
    cells = new Cell * [grid_size];
    mask = new bool [grid_size];

    for (int cell = 0; cell < grid_size; cell ++)
    {
        // for cells
        cells[cell] = new Cell (grid_corners[cell], 
                grid_center_lat[cell], grid_center_lon[cell], 
                crnr_lat, crnr_lon);
        crnr_lat += grid_corners[cell];
        crnr_lon += grid_corners[cell];
        
        // for masks
        mask[cell] = (cell_mask[cell] == 1);
    }

    // initialize search
    init_search();
}

Cell ** Grid::getCell()
{
    return cells;
}

int * Grid::getMask()
{
    return imask;
}

void Grid::init_search()
{
    // initialize latitude range for search bins
    bin_lats = new double* [Constant::num_srch_bins];
    for (int i = 0; i < Constant::num_srch_bins; i++)
        bin_lats[i] = new double[2];
    double dlat = Constant::PI / Constant::num_srch_bins;
    for (int i = 0; i < Constant::num_srch_bins; i++)
    {
        bin_lats[i][0] = - Constant::PI2 + i * dlat;
        bin_lats[i][1] = - Constant::PI2 + (i + 1) * dlat;
    }

    // initialize index range for each cell
    bin_addr    = new int * [size];
    for (int cell = 0; cell < size; cell ++)
    {
        bin_addr[cell] = new int[2];
        bin_addr[cell][0] = -1;
        bin_addr[cell][1] = -1;
    }
}

void Grid::log()
{
    printf("log from class Grid\n");
}

/*----------------------------------------*/


/** For Public   class   <GCUtils>    **/

GCUtils::GCUtils()
{
}

GCUtils::~GCUtils()
{
}

PointPosition GCUtils::getPointPosition(Point *p, Point *q)
{
    return NorthEast;
}

bool GCUtils::isCellContainsPoint(Cell * cell, Point * point)
{
    return false;
}

/*----------------------------------------*/

