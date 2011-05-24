#ifndef _GREAT_CIRCLE_BASE_H_
#define _GREAT_CIRCLE_BASE_H_ 1

#include <stdio.h>

/** Define some macros for units **/
#define     DEGREES     0
#define     RADIANS     1
/*----------------------------------------*/

/** Define some macros for readablity **/

#define     LAT_MIN     0
#define     LAT_MAX     1
#define     LON_MIN     2
#define     LON_MAX     3
/*----------------------------------------*/

/** Define some macros for Point relative postions **/

typedef int PointPosition;
#define     NorthEast   0
#define     NorthWest   1
#define     SouthEast   2
#define     SouthWest   3
/*----------------------------------------*/

/** Define some macros for Point connection way **/

typedef int PointConnection;
#define     CONN_GREAT_CIRCLE       0
#define     CONN_LATLON_LINEAR      1
#define     CONN_DEFAULT_WAY        1
/*----------------------------------------*/

/** Define some macros for typedef **/

/*----------------------------------------*/

/*
 *  Static class for constants
 */
class Constant
{
    public:
        static const double ZERO    = 0.0;
        static const double HALF    = 0.5;
        static const double ONE     = 1.0;
        static const double TWO     = 2.0;
        static const double THREE   = 3.0;
        static const double FOUR    = 4.0;
        static const double FIVE    = 5.0;
        static const double PI      = 3.14159265359;
        static const double PI2     = PI * 2;
        static const double PIH     = PI / 2;
        static const double TINY    = 1.E-14;
        static const double BIGNUM  = 1.E+20;
        static const int    CHARLEN = 80;
        static const double deg2rad = PI / 180;
        static const int    num_srch_bins = 300;
};

/*
 *  Basic class for Great Circle
 */
class GCBase
{
    public:
        // construction function
        GCBase();

        //deconstruction function
        ~GCBase();

        // log system for details checking
        virtual void log();
};
/*----------------------------------------*/

/* 
 *  Public class for Point on surface
 *  Two coordination representation: lat/lon-coord and xyz-coord
 */
class Point : public GCBase
{
    friend class GCUtils;
    friend class Cell;
    friend class Edge;
    
    private:
        double lat;
        double lon;
        double x;
        double y;
        double z;
    public:
        static const double Radius = 100.0;
        Point(double latitude, double longitude);
        Point(double x_coord, double y_coord, double z_coord);
        ~Point();
        void log();
};
/*----------------------------------------*/

/*
 *  Public class for Segment on surface
 *  including start Point, end Point and normal vector
 */
class Edge : public GCBase
{
    friend class GCUtils;
    friend class Cell;
    private:
        Point * start;              // start point of edge
        Point * end;                // end point of edge
        Point * norm_vec;           // normal vector of edge
        PointConnection conn;       // connection way of two Point (great circle or lat/lon)
    public:
        // construction function
        Edge(double beglat, double beglon, double endlat, double endlon, PointConnection connection);
        
        // deconstruction function
        ~Edge();

        // print Edge information
        void log();                 
};
/*----------------------------------------*/

/*
 *  Public class for Cell in Grid
 */
class Cell : public GCBase
{
    friend class GCUtils;
    friend class Grid;
    private:
        int     size;           // cell size, number of cell corners
        Point   ** corners;     // array of cell corners
        Point   * center;       // cell center
        Edge    ** edges;       // array of cell edges
        double  *  bbox;        // bounding box for restricting search
        double  area;           // cell area
        double  frac;           // cell frac
        PointConnection * conns;// Point connection way of for each Edge; if null, use default

        // calculate bound box for Grid
        void calc_bbox();       // calculate bounding box

    public:
        // construction function
        Cell(int cell_size, double cntr_lat, double cntr_lon, double * crnr_lat, double * crnr_lon, PointConnection connect);

        // deconstruction function
        ~Cell();
        
        // print Cell information
        void log();                     
};
/*----------------------------------------*/

/*
 *  Public class for Grid
 */
class Grid : public GCBase
{
    friend class GCUtils;
    private:
        char        name[Constant::CHARLEN];           
                                    // grid name string
        int         rank;           // grid rank (number of dimensions)
        int         size;           // grid size
        int     *   dims;           // grid size for each dimension
        int         units;          // units of grid (degrees or radians)
        Cell    **  cells;          // array of grid points in 1D struct
        Cell    *** cells_2D;       // array of grid points in 2D struct
        Cell    ****cells_3D;       // array of grid points in 3D struct
        double  *   areas_in;       // array of grid areas
        double  *   fracs_in;       // array of grid fracs
        bool    *   mask;           // array of grid masks
        double  **  bin_lats;       // latitude range of search bins
        double  **  bin_lons;       // longitude range of search bins
        int     **  bin_addr;       // cell index bound for search bins
        PointConnection ** conns;   // Point connection way for each Cell; if null, use default


    public:
        // default null construction function
        Grid();
        
        // construction function without area/frac input
        Grid(int grid_rank, int grid_size, int * grid_dims, int * grid_corners, int grid_corners_max, double * grid_center_lat, double * grid_center_lon, double * grid_corner_lat, double * grid_corner_lon, int * grid_masks);
        
        // construction function with area/frac input
        Grid(int grid_rank, int * grid_size, int * cell_size, double * center_lat, double * center_lon, double * corner_lat, double * corner_lon, int * cell_mask, double * cell_area_in, double * cell_frac_in);

        // deconstruction function
        ~Grid();

        /* Member functions */
        Cell ** getCell();
        int *   getMask();
        void init_search();

        // print Grid information
        void log();
};
/*----------------------------------------*/

/*
 *  Public class for Great Circle Utils 
 */
class GCUtils
{
    public:
        // construction function
        GCUtils();

        // deconstruction function
        ~GCUtils();
    
        // calculate cross procduct of two Points
        static Point * cross_product(Point *left, Point * right);

        // check relative postion of Point p and Point q
        static PointPosition getPointPosition(Point * p, Point * q);

        // judge whether Cell contains Point
        static bool isCellContainsPoint(Cell * cell, Point * point);

        // judge whether Cells intersect based on bounding box
        static bool isCellIntersectCell(Cell * cell, Cell * llec);

        // calculate search bins for source and destination Grid
        static void CalcSearchBins(Grid * src, Grid * dst);

};
/*----------------------------------------*/

#endif
