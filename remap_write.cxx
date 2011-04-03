#include "remap_write.h"

#include "kinds.h"
#include "constants.h"
#include "grids.h"
#include "remap_vars.h"
#include "nc_error.h"
#include "utils.h"
#include "namelist.h"
#include "io.h"

// variables used only in remap_write.o
char map_method_write[CHARLEN]; // character string for map_type
char norm_opt_write[CHARLEN];   // character string for normalization option
char history[CHARLEN];          // character string for history infomation
char convention[CHARLEN];       // character string for output convention

char cdate[12];                 // character data string

int *src_mask_int;              // integer masks to determine cells that participate in map
int *dst_mask_int;

// various netCDF identifiers used by output routines
int ncstat;                     // error flag for netCDF calls
int nc_file_id;                 // id for netCDF file
int nc_srcgrdsize_id;           // id for source grid size
int nc_dstgrdsize_id;           // id for destination grid size
int nc_srcgrdcorner_id;         // id for max number of source grid corners
int nc_dstgrdcorner_id;         // id for max number of destination grid corners
int nc_srcgrdrank_id;           // id for source grid rank
int nc_dstgrdrank_id;           // id for destination grid rank
int nc_numlinks_id;             // id for number of links in remapping
int nc_numwgts_id;              // id for number of weights for mapping
int nc_srcgrddims_id;           // id for source grid dimensions
int nc_dstgrddims_id;           // id for destination grid dimensions
int nc_srcgrdcntrlat_id;        // id for source grid center latitude
int nc_dstgrdcntrlat_id;        // id for destination grid center latitude
int nc_srcgrdcntrlon_id;        // id for source grid center longitude
int nc_dstgrdcntrlon_id;        // id for destination grid center longitude
int nc_srcgrdimask_id;          // id for source grid mask
int nc_dstgrdimask_id;          // id for destination grid mask
int nc_srcgrdcrnrlat_id;        // id for source grid corner latitude
int nc_dstgrdcrnrlat_id;        // id for destination grid corner latitude
int nc_srcgrdcrnrlon_id;        // id for source grid corner longitude
int nc_dstgrdcrnrlon_id;        // id for destination grid corner longitude
int nc_srcgrdarea_id;           // id for area of source grid cells
int nc_dstgrdarea_id;           // id for area of destination grid cells
int nc_srcgrdfrac_id;           // id for frac of source grid cells
int nc_dstgrdfrac_id;           // id for frac of destination grid cells
int nc_srcadd_id;               // id for map source address
int nc_dstadd_id;               // id for map destination address
int nc_rmpmatrix_id;            // id for remapping matrix

int *nc_dims2_id;               // netCDF ids for 2d array dims


// call correct output routine based on output format choice (scrip or csm)
void write_remap(char *map_name, char *interp_file_str, int output_opt)
{
    log("Write_Remap");

    // define some common variables to be used in all routines
    switch (norm_opt)           // normalization option
    {
        case NORM_OPT_NONE:
            strcpy(norm_opt_write, "none");
            break;
        case NORM_OPT_FRACAREA:
            strcpy(norm_opt_write, "fracarea");
            break;
        case  NORM_OPT_DESTAREA:
            strcpy(norm_opt_write, "destarea");
            break;
        default:
            log("Invalid Normalization Option Type");
            return;
    }
    
    switch (map_type)           // map method
    {
        case MAP_TYPE_CONSERV:
            strcpy(map_method_write, "Conservative remapping");
            break;
        case MAP_TYPE_BILINEAR:
            strcpy(map_method_write, "Bilinear remapping");
            break;
        case MAP_TYPE_BICUBIC:
            strcpy(map_method_write, "Bicubic remapping");
            break;
        case MAP_TYPE_DISTWGT:
            strcpy(map_method_write, "Distance weighted avg of nearest neighbors");
            break;
        default:
            log("Invaild Map Type");
            return;
    }

                                // date & history
    sysdate(cdate);
    strcpy(history, "Created: ");
    strcat(history, cdate);


    // sort address and weight arrays
    sort_add(grid2_add_map, grid1_add_map, wts_map, num_links_map, num_wts);

    // call appropriate output routine
    switch (output_opt)
    {
        case SCRIP_CONVENTION:
            write_remap_scrip(map_name, interp_file_str, 1);
            break;
        case CSM_CONVENTION:
            write_remap_csm  (map_name, interp_file_str, 1);
            break;
        default:
            log("Unknown output file convention");
            return;
    }

}

// writes remap data to a netCDF file using SCRIP conventions
void write_remap_scrip(char *map_name, char *interp_file_str, int direction)
{
    log("--Using SCRIP convention");
    
    // temp variables
    char grid1_ctmp[CHARLEN];       // character temp for grid1 names
    char grid2_ctmp[CHARLEN];       // character temp for grid2 names

    int itmp1;                      // integer temp
    int itmp2;                      // integer temp
    int itmp3;                      // integer temp
    int itmp4;                      // integer temp

    /* create netCDF file for mapping and define some global attributes */
    
    // create netCDF file named interp_file
    ncstat = nc_create (interp_file_str, NC_CLOBBER, &nc_file_id);
    ERR

    // map name
    ncstat = nc_put_att_text (nc_file_id, NC_GLOBAL, "title", len_trim(map_name), map_name);
    ERR

    // normalization option
    ncstat = nc_put_att_text (nc_file_id, NC_GLOBAL, "normalization", len_trim(norm_opt_write), norm_opt_write);
    ERR

    // map method
    ncstat = nc_put_att_text (nc_file_id, NC_GLOBAL, "map_method", len_trim(map_method_write), map_method_write);
    ERR

    // history
    ncstat = nc_put_att_text (nc_file_id, NC_GLOBAL, "history", len_trim(history), history);
    ERR

    // file convention
    strcpy(convention, "SCRIP CONVENTION");
    ncstat = nc_put_att_text (nc_file_id, NC_GLOBAL, "conventions", len_trim(convention), convention);
    ERR

    // source and destination grid names
    if (direction == 1)
    {
        strcpy(grid1_ctmp, "source_grid");
        strcpy(grid2_ctmp, "dest_grid");
    }
    else
    {
        strcpy(grid2_ctmp, "source_grid");
        strcpy(grid1_ctmp, "dest_grid");
    }
    ncstat = nc_put_att_text (nc_file_id, NC_GLOBAL, trim(grid1_ctmp), len_trim(grid1_name), grid1_name);
    ERR
    ncstat = nc_put_att_text (nc_file_id, NC_GLOBAL, trim(grid2_ctmp), len_trim(grid2_name), grid2_name);
    ERR

    /* prepare netCDF dimention info */

    // define grid size dimentions
    if (direction == 1)
    {
        itmp1 = grid1_size;
        itmp2 = grid2_size;
    }
    else
    {
        itmp1 = grid2_size;
        itmp2 = grid1_size;
    }
    ncstat = nc_def_dim (nc_file_id, "src_grid_size", itmp1, &nc_srcgrdsize_id);
    ERR
    ncstat = nc_def_dim (nc_file_id, "dst_grid_size", itmp2, &nc_dstgrdsize_id);
    ERR
    
    // define grid corner dimension
    if (direction == 1)
    {
        itmp1 = grid1_corners_max;
        itmp2 = grid2_corners_max;
    }
    else
    {
        itmp2 = grid1_corners_max;
        itmp1 = grid2_corners_max;
    }
    ncstat = nc_def_dim (nc_file_id, "src_grid_corners", itmp1, &nc_srcgrdcorner_id);
    ERR
    ncstat = nc_def_dim (nc_file_id, "dst_grid_corners", itmp2, &nc_dstgrdcorner_id);
    ERR

    // define grid rank dimensions
    if (direction == 1)
    {
        itmp1 = grid1_rank;
        itmp2 = grid2_rank;
    }
    else
    {
        itmp1 = grid2_rank;
        itmp2 = grid1_rank;
    }
    ncstat = nc_def_dim (nc_file_id, "src_grid_rank", itmp1, &nc_srcgrdrank_id);
    ERR
    ncstat = nc_def_dim (nc_file_id, "dst_grid_rank", itmp2, &nc_dstgrdrank_id);
    ERR

    // define map size dimensions
    ncstat = nc_def_dim (nc_file_id, "num_links", num_links_map, &nc_numlinks_id);
    ERR

    ncstat = nc_def_dim (nc_file_id, "num_wgts", num_wts, &nc_numwgts_id);
    ERR

    // define grid dimensions
    ncstat = nc_def_var (nc_file_id, "src_grid_dims", NC_INT, 1, &nc_srcgrdrank_id, &nc_srcgrddims_id);
    ERR
    ncstat = nc_def_var (nc_file_id, "dst_grid_dims", NC_INT, 1, &nc_dstgrdrank_id, &nc_dstgrddims_id);
    ERR

    /* define all arrays for netCDF descriptors */
    
    // define grid center latitude array
    ncstat = nc_def_var (nc_file_id, "src_grid_center_lat", NC_DOUBLE, 1, &nc_srcgrdsize_id, &nc_srcgrdcntrlat_id);
    ERR
    ncstat = nc_def_var (nc_file_id, "dst_grid_center_lat", NC_DOUBLE, 1, &nc_dstgrdsize_id, &nc_dstgrdcntrlat_id);
    ERR

    // define grid center longitude array
    ncstat = nc_def_var (nc_file_id, "src_grid_center_lon", NC_DOUBLE, 1, &nc_srcgrdsize_id, &nc_srcgrdcntrlon_id);
    ERR
    ncstat = nc_def_var (nc_file_id, "dst_grid_center_lon", NC_DOUBLE, 1, &nc_dstgrdsize_id, &nc_dstgrdcntrlon_id);
    ERR

    // define grid corner lat/lon arrays
    nc_dims2_id = new int [2];

    // define source grid corner latitude|longitude array
    nc_dims2_id[1] = nc_srcgrdcorner_id;
    nc_dims2_id[0] = nc_srcgrdsize_id;
    ncstat = nc_def_var (nc_file_id, "src_grid_corner_lat", NC_DOUBLE, 2, nc_dims2_id, &nc_srcgrdcrnrlat_id);
    ERR
    ncstat = nc_def_var (nc_file_id, "src_grid_corner_lon", NC_DOUBLE, 2, nc_dims2_id, &nc_srcgrdcrnrlon_id);
    ERR

    // define destination grid corner latitude|longitude array
    nc_dims2_id[1] = nc_dstgrdcorner_id;
    nc_dims2_id[0] = nc_dstgrdsize_id;
    ncstat = nc_def_var (nc_file_id, "dst_grid_corner_lat", NC_DOUBLE, 2, nc_dims2_id, &nc_dstgrdcrnrlat_id);
    ERR
    ncstat = nc_def_var (nc_file_id, "dst_grid_corner_lon", NC_DOUBLE, 2, nc_dims2_id, &nc_dstgrdcrnrlon_id);
    ERR

    // define units for all coordinate arrays
    if (direction == 1)
    {
        strcpy(grid1_ctmp, grid1_units);
        strcpy(grid2_ctmp, grid2_units);
    }
    else
    {
        strcpy(grid1_ctmp, grid2_units);
        strcpy(grid2_ctmp, grid1_units);
    }

    ncstat = nc_put_att_text (nc_file_id, nc_srcgrdcntrlat_id, "units", 7, grid1_ctmp);
    ERR
    ncstat = nc_put_att_text (nc_file_id, nc_dstgrdcntrlat_id, "units", 7, grid2_ctmp);
    ERR
    ncstat = nc_put_att_text (nc_file_id, nc_srcgrdcntrlon_id, "units", 7, grid1_ctmp);
    ERR
    ncstat = nc_put_att_text (nc_file_id, nc_dstgrdcntrlon_id, "units", 7, grid2_ctmp);
    ERR
    ncstat = nc_put_att_text (nc_file_id, nc_srcgrdcrnrlat_id, "units", 7, grid1_ctmp);
    ERR
    ncstat = nc_put_att_text (nc_file_id, nc_dstgrdcrnrlat_id, "units", 7, grid2_ctmp);
    ERR
    ncstat = nc_put_att_text (nc_file_id, nc_srcgrdcrnrlon_id, "units", 7, grid1_ctmp);
    ERR
    ncstat = nc_put_att_text (nc_file_id, nc_dstgrdcrnrlon_id, "units", 7, grid2_ctmp);
    ERR
    
    // define grid mask
    ncstat = nc_def_var (nc_file_id, "src_grid_mask", NC_INT, 1, &nc_srcgrdsize_id, &nc_srcgrdimask_id);
    ERR
    ncstat = nc_put_att_text (nc_file_id, nc_srcgrdimask_id, "units", 8, "unitless");
    ERR
    ncstat = nc_def_var (nc_file_id, "dst_grid_mask", NC_INT, 1, &nc_dstgrdsize_id, &nc_dstgrdimask_id);
    ERR
    ncstat = nc_put_att_text (nc_file_id, nc_dstgrdimask_id, "units", 8, "unitless");

    // define grid area arrays
    ncstat = nc_def_var (nc_file_id, "src_grid_area", NC_DOUBLE, 1, &nc_srcgrdsize_id, &nc_srcgrdarea_id);
    ERR
    ncstat = nc_put_att_text (nc_file_id, nc_srcgrdarea_id, "units", 14, "square radians");
    ERR
    ncstat = nc_def_var (nc_file_id, "dst_grid_area", NC_DOUBLE, 1, &nc_dstgrdsize_id, &nc_dstgrdarea_id);
    ERR
    ncstat = nc_put_att_text (nc_file_id, nc_dstgrdarea_id, "units", 14, "square radians");
    ERR

    // define grid fraction arrays
    ncstat = nc_def_var (nc_file_id, "src_grid_frac", NC_DOUBLE, 1, &nc_srcgrdsize_id, &nc_srcgrdfrac_id);
    ERR
    ncstat = nc_put_att_text (nc_file_id, nc_srcgrdfrac_id, "units", 8, "unitless");
    ERR
    ncstat = nc_def_var (nc_file_id, "dst_grid_frac", NC_DOUBLE, 1, &nc_dstgrdsize_id, &nc_dstgrdfrac_id);
    ERR
    ncstat = nc_put_att_text (nc_file_id, nc_dstgrdfrac_id, "units", 8, "unitless");
    ERR

    // define mapping arrays
    ncstat = nc_def_var (nc_file_id, "src_address", NC_INT, 1, &nc_numlinks_id, &nc_srcadd_id);
    ERR
    ncstat = nc_def_var (nc_file_id, "dst_address", NC_INT, 1, &nc_numlinks_id, &nc_dstadd_id);
    ERR

    nc_dims2_id = new int[2];
    nc_dims2_id[1] = nc_numwgts_id;
    nc_dims2_id[0] = nc_numlinks_id;

    ncstat = nc_def_var (nc_file_id, "remap_maxtrix", NC_DOUBLE, 2, nc_dims2_id, &nc_rmpmatrix_id);
    ERR


    /* end definition stage */
    ncstat = nc_enddef(nc_file_id);
    ERR


    /* compute integer masks */
    if (direction == 1)
    {
        src_mask_int = new int [grid1_size];
        dst_mask_int = new int [grid2_size];
        for (int i = 0; i < grid1_size; i++)
        {
            if (grid1_mask[i])
                src_mask_int[i] = 1;
            else
                src_mask_int[i] = 0;
        }
        for (int i = 0; i < grid2_size; i++)
        {
            if (grid2_mask[i])
                dst_mask_int[i] = 1;
            else
                dst_mask_int[i] = 0;
        }
    }
    else
    {
        src_mask_int = new int [grid2_size];
        dst_mask_int = new int [grid1_size];
        for (int i = 0; i < grid2_size; i++)
        {
            if (grid2_mask[i])
                src_mask_int[i] = 1;
            else
                src_mask_int[i] = 0;
        }
        for (int i = 0; i < grid1_size; i++)
        {
            if (grid1_mask[i])
                dst_mask_int[i] = 1;
            else
                dst_mask_int[i] = 0;
        }
    }


    /* change units of lat/lon coordinates if input units different from radians */
    if (streqls(grid1_units, "degrees", 0, 7) && direction == 1)
    {
        for (int i = 0; i < grid1_size; i++)
        {
            grid1_center_lat[i] /= deg2rad;
            grid1_center_lon[i] /= deg2rad;
        }
        int num_corners1 = grid1_size * grid1_corners_max;
        for (int i = 0; i < num_corners1; i++)
        {
            grid1_corner_lat[i] /= deg2rad;
            grid1_corner_lon[i] /= deg2rad;
        }
    }
    if (streqls(grid2_units, "degrees", 0, 7) && direction == 1)
    {
        for (int i = 0; i < grid2_size; i++)
        {
            grid2_center_lat[i] /= deg2rad;
            grid2_center_lon[i] /= deg2rad;
        }
        int num_corners2 = grid2_size * grid2_corners_max;
        for (int i = 0; i < num_corners2; i++)
        {
            grid2_corner_lat[i] /= deg2rad;
            grid2_corner_lon[i] /= deg2rad;
        }
    }


    /* write mapping data */
    // write mask
    if (direction == 1)
    {
        itmp1 = nc_srcgrddims_id;
        itmp2 = nc_dstgrddims_id;
    }
    else
    {
        itmp2 = nc_srcgrddims_id;
        itmp1 = nc_dstgrddims_id;
    }

    ncstat = nc_put_var_int (nc_file_id, itmp1, grid1_dims);
    ERR
    ncstat = nc_put_var_int (nc_file_id, itmp2, grid2_dims);
    ERR

    ncstat = nc_put_var_int (nc_file_id, nc_srcgrdimask_id, src_mask_int);
    ERR
    ncstat = nc_put_var_int (nc_file_id, nc_dstgrdimask_id, dst_mask_int);
    ERR
    // delete mask_int
    delete [] src_mask_int;
    delete [] dst_mask_int;

    // write center|corner latitude|longitude
    // write grid1 center|corner lat|lon
    if (direction == 1)
    {
        itmp1 = nc_srcgrdcntrlat_id;
        itmp2 = nc_srcgrdcntrlon_id;
        itmp3 = nc_srcgrdcrnrlat_id;
        itmp4 = nc_srcgrdcrnrlon_id;
    }
    else
    {
        itmp1 = nc_dstgrdcntrlat_id;
        itmp2 = nc_dstgrdcntrlon_id;
        itmp3 = nc_dstgrdcrnrlat_id;
        itmp4 = nc_dstgrdcrnrlon_id;
    }

    ncstat = nc_put_var_double (nc_file_id, itmp1, grid1_center_lat);
    ERR
    ncstat = nc_put_var_double (nc_file_id, itmp2, grid1_center_lon);
    ERR
    ncstat = nc_put_var_double (nc_file_id, itmp3, grid1_corner_lat);
    ERR
    ncstat = nc_put_var_double (nc_file_id, itmp4, grid1_corner_lon);
    ERR
    // write grid2 center|corner lat|lon
    if (direction == 1)
    {
        itmp1 = nc_dstgrdcntrlat_id;
        itmp2 = nc_dstgrdcntrlon_id;
        itmp3 = nc_dstgrdcrnrlat_id;
        itmp4 = nc_dstgrdcrnrlon_id;
    }
    else
    {
        itmp1 = nc_srcgrdcntrlat_id;
        itmp2 = nc_srcgrdcntrlon_id;
        itmp3 = nc_srcgrdcrnrlat_id;
        itmp4 = nc_srcgrdcrnrlon_id;
    }

    ncstat = nc_put_var_double (nc_file_id, itmp1, grid2_center_lat);
    ERR
    ncstat = nc_put_var_double (nc_file_id, itmp2, grid2_center_lon);
    ERR
    ncstat = nc_put_var_double (nc_file_id, itmp3, grid2_corner_lat);
    ERR
    ncstat = nc_put_var_double (nc_file_id, itmp4, grid2_corner_lon);
    ERR

    // write area and frac
    if (direction == 1)
    {
        itmp1 = nc_srcgrdarea_id;
        itmp2 = nc_srcgrdfrac_id;
        itmp3 = nc_dstgrdarea_id;
        itmp4 = nc_dstgrdfrac_id;
    }
    else
    {
        itmp1 = nc_dstgrdarea_id;
        itmp2 = nc_dstgrdfrac_id;
        itmp3 = nc_srcgrdarea_id;
        itmp4 = nc_srcgrdfrac_id;
    }

    if (luse_grid1_area)
        ncstat = nc_put_var_double (nc_file_id, itmp1, grid1_area_in);
    else
        ncstat = nc_put_var_double (nc_file_id, itmp1, grid1_area);
    ERR

    if (luse_grid2_area)
        ncstat = nc_put_var_double (nc_file_id, itmp2, grid2_area_in);
    else
        ncstat = nc_put_var_double (nc_file_id, itmp2, grid2_area);
    ERR
    
    ncstat = nc_put_var_double (nc_file_id, itmp3, grid1_frac);
    ERR
    ncstat = nc_put_var_double (nc_file_id, itmp4, grid2_frac);
    ERR

    // write src|dst grid add, remap matrix of mappings
    ncstat = nc_put_var_int (nc_file_id, nc_srcadd_id, grid1_add_map);
    ERR
    ncstat = nc_put_var_int (nc_file_id, nc_dstadd_id, grid2_add_map);
    ERR
    ncstat = nc_put_var_double (nc_file_id, nc_rmpmatrix_id, wts_map);
    ERR

    // delete nc_dims_id array
    delete [] nc_dims2_id;

    /* close netCDF file */
    nc_close (nc_file_id);
    ERR
}

// writes remap data to a netCDF file using NCAR-CSM conventions
void write_remap_csm(char *map_name, char *interp_file_str, int direction)
{
    log("--Using NCAR-CSM convention");
    char grid1_ctmp[CHARLEN];       // character temp for grid1 names
    char grid2_ctmp[CHARLEN];       // character temp for grid2 names

    int itmp1, itmp2, itmp3, itmp4; // integer temp
    int nc_numwgts1_id;              // extra netCDF id for additional weights
    int nc_src_isize_id;            // extra netCDF id for ni_a
    int nc_src_jsize_id;            // extra netCDF id for nj_a
    int nc_dst_isize_id;            // extra netCDF id for ni_b
    int nc_dst_jsize_id;            // extra netCDF id for nj_b
    int nc_rmpmatrix2_id;           // extra netCDF id for high-order remap matrix

    double *wts1;                   // CSM wants single array for 1st-order wts
    double *wts2;                   // write remaining weights in different array

    /* create netCDF file for mapping and define some global attributes */
    // create netCDF file named interp_file
    ncstat = nc_create (interp_file_str, NC_CLOBBER, &nc_file_id);
    ERR

    // map name
    ncstat = nc_put_att_text (nc_file_id, NC_GLOBAL, "title", len_trim(map_name), map_name);
    ERR

    // normalization option
    ncstat = nc_put_att_text (nc_file_id, NC_GLOBAL, "normalization", len_trim(norm_opt_write), norm_opt_write);
    ERR

    // map method
    ncstat = nc_put_att_text (nc_file_id, NC_GLOBAL, "map_method", len_trim(map_method_write), map_method_write);
    ERR

    // history
    ncstat = nc_put_att_text (nc_file_id, NC_GLOBAL, "history", len_trim(history), history);
    ERR

    // file convention
    strcpy(convention, "NCAR-CSM CONVENTION");
    ncstat = nc_put_att_text (nc_file_id, NC_GLOBAL, "conventions", len_trim(convention), convention);
    ERR

    // source and destination grid names
    if (direction == 1)
    {
        strcpy(grid1_ctmp, "domain_a");
        strcpy(grid2_ctmp, "domain_b");
    }
    else
    {
        strcpy(grid2_ctmp, "domain_a");
        strcpy(grid1_ctmp, "domain_b");
    }
    ncstat = nc_put_att_text (nc_file_id, NC_GLOBAL, trim(grid1_ctmp), len_trim(grid1_name), grid1_name);
    ERR
    ncstat = nc_put_att_text (nc_file_id, NC_GLOBAL, trim(grid2_ctmp), len_trim(grid2_name), grid2_name);
    ERR

    /* prepare netCDF dimention info */

    // define grid size dimentions
    if (direction == 1)
    {
        itmp1 = grid1_size;
        itmp2 = grid2_size;
    }
    else
    {
        itmp1 = grid2_size;
        itmp2 = grid1_size;
    }
    ncstat = nc_def_dim (nc_file_id, "n_a", itmp1, &nc_srcgrdsize_id);
    ERR
    ncstat = nc_def_dim (nc_file_id, "n_b", itmp2, &nc_dstgrdsize_id);
    ERR
    
    // define grid corner dimension
    if (direction == 1)
    {
        itmp1 = grid1_corners_max;
        itmp2 = grid2_corners_max;
    }
    else
    {
        itmp1 = grid2_corners_max;
        itmp2 = grid1_corners_max;
    }
    ncstat = nc_def_dim (nc_file_id, "nv_a", itmp1, &nc_srcgrdcorner_id);
    ERR
    ncstat = nc_def_dim (nc_file_id, "nv_b", itmp2, &nc_dstgrdcorner_id);
    ERR

    // define grid rank dimensions
    if (direction == 1)
    {
        itmp1 = grid1_rank;
        itmp2 = grid2_rank;
    }
    else
    {
        itmp1 = grid2_rank;
        itmp2 = grid1_rank;
    }
    ncstat = nc_def_dim (nc_file_id, "src_grid_rank", itmp1, &nc_srcgrdrank_id);
    ERR
    ncstat = nc_def_dim (nc_file_id, "dst_grid_rank", itmp2, &nc_dstgrdrank_id);
    ERR

    // define first two dims as if 2-d cartesian domain
    if (direction == 1)
    {
        itmp1 = grid1_dims[0];
        itmp3 = grid2_dims[0];
        if (grid1_rank > 1)
            itmp2 = grid1_dims[1];
        else
            itmp2 = 0;
        if (grid2_rank > 1)
            itmp4 = grid2_dims[1];
        else
            itmp4 = 0;
    }
    else
    {
        itmp3 = grid1_dims[0];
        itmp1 = grid2_dims[0];
        if (grid1_rank > 1)
            itmp4 = grid1_dims[1];
        else
            itmp4 = 0;
        if (grid2_rank > 1)
            itmp2 = grid2_dims[1];
        else
            itmp2 = 0;
    }

    ncstat = nc_def_dim (nc_file_id, "ni_a", itmp1, &nc_src_isize_id);
    ERR
    ncstat = nc_def_dim (nc_file_id, "nj_a", itmp2, &nc_src_jsize_id);
    ERR
    ncstat = nc_def_dim (nc_file_id, "ni_b", itmp3, &nc_src_isize_id);
    ERR
    ncstat = nc_def_dim (nc_file_id, "nj_b", itmp4, &nc_src_jsize_id);
    ERR

    // define map size dimensions
    ncstat = nc_def_dim (nc_file_id, "n_s", num_links_map, &nc_numlinks_id);
    ERR
    ncstat = nc_def_dim (nc_file_id, "num_wgts", num_wts, &nc_numwgts_id);
    ERR
    if (num_wts > 1)
    {
        ncstat = nc_def_dim (nc_file_id, "num_wgts1", num_wts-1, &nc_numwgts1_id);
        ERR
    }

    /* define all arrays for netCDF descriptors */
    
    // define grid center latitude array
    ncstat = nc_def_var (nc_file_id, "yc_a", NC_DOUBLE, 1, &nc_srcgrdsize_id, &nc_srcgrdcntrlat_id);
    ERR
    ncstat = nc_def_var (nc_file_id, "yc_b", NC_DOUBLE, 1, &nc_dstgrdsize_id, &nc_dstgrdcntrlat_id);
    ERR

    // define grid center longitude array
    ncstat = nc_def_var (nc_file_id, "xc_a", NC_DOUBLE, 1, &nc_srcgrdsize_id, &nc_srcgrdcntrlon_id);
    ERR
    ncstat = nc_def_var (nc_file_id, "xc_b", NC_DOUBLE, 1, &nc_dstgrdsize_id, &nc_dstgrdcntrlon_id);
    ERR

    // define grid corner latitude array
    ncstat = nc_def_var (nc_file_id, "yv_a", NC_DOUBLE, 1, &nc_srcgrdsize_id, &nc_srcgrdcrnrlat_id);
    ERR
    ncstat = nc_def_var (nc_file_id, "yv_b", NC_DOUBLE, 1, &nc_dstgrdsize_id, &nc_dstgrdcrnrlat_id);
    ERR

    // define grid corner longitude array
    ncstat = nc_def_var (nc_file_id, "xv_a", NC_DOUBLE, 1, &nc_srcgrdsize_id, &nc_srcgrdcrnrlon_id);
    ERR
    ncstat = nc_def_var (nc_file_id, "xv_b", NC_DOUBLE, 1, &nc_dstgrdsize_id, &nc_dstgrdcrnrlon_id);
    ERR

    /* CSM wants all in degrees */
    strcpy(grid1_units, "degrees");
    strcpy(grid2_units, "degrees");

    // define units for all coordinate arrays
    if (direction == 1)
    {
        strcpy(grid1_ctmp, grid1_units);
        strcpy(grid2_ctmp, grid2_units);
    }
    else
    {
        strcpy(grid1_ctmp, grid2_units);
        strcpy(grid2_ctmp, grid1_units);
    }

    ncstat = nc_put_att_text (nc_file_id, nc_srcgrdcntrlat_id, "units", 7, grid1_ctmp);
    ERR
    ncstat = nc_put_att_text (nc_file_id, nc_dstgrdcntrlat_id, "units", 7, grid2_ctmp);
    ERR
    ncstat = nc_put_att_text (nc_file_id, nc_srcgrdcntrlon_id, "units", 7, grid1_ctmp);
    ERR
    ncstat = nc_put_att_text (nc_file_id, nc_dstgrdcntrlon_id, "units", 7, grid2_ctmp);
    ERR
    ncstat = nc_put_att_text (nc_file_id, nc_srcgrdcrnrlat_id, "units", 7, grid1_ctmp);
    ERR
    ncstat = nc_put_att_text (nc_file_id, nc_dstgrdcrnrlat_id, "units", 7, grid2_ctmp);
    ERR
    ncstat = nc_put_att_text (nc_file_id, nc_srcgrdcrnrlon_id, "units", 7, grid1_ctmp);
    ERR
    ncstat = nc_put_att_text (nc_file_id, nc_dstgrdcrnrlon_id, "units", 7, grid2_ctmp);
    ERR
    
    // define grid mask
    ncstat = nc_def_var (nc_file_id, "mask_a", NC_INT, 1, &nc_srcgrdsize_id, &nc_srcgrdimask_id);
    ERR
    ncstat = nc_put_att_text (nc_file_id, nc_srcgrdimask_id, "units", 8, "unitless");
    ERR
    ncstat = nc_def_var (nc_file_id, "mask_b", NC_INT, 1, &nc_dstgrdsize_id, &nc_dstgrdimask_id);
    ERR
    ncstat = nc_put_att_text (nc_file_id, nc_dstgrdimask_id, "units", 8, "unitless");

    // define grid area arrays
    ncstat = nc_def_var (nc_file_id, "area_a", NC_DOUBLE, 1, &nc_srcgrdsize_id, &nc_srcgrdarea_id);
    ERR
    ncstat = nc_put_att_text (nc_file_id, nc_srcgrdarea_id, "units", 14, "square radians");
    ERR
    ncstat = nc_def_var (nc_file_id, "area_b", NC_DOUBLE, 1, &nc_dstgrdsize_id, &nc_dstgrdarea_id);
    ERR
    ncstat = nc_put_att_text (nc_file_id, nc_dstgrdarea_id, "units", 14, "square radians");
    ERR

    // define grid fraction arrays
    ncstat = nc_def_var (nc_file_id, "frac_a", NC_DOUBLE, 1, &nc_srcgrdsize_id, &nc_srcgrdfrac_id);
    ERR
    ncstat = nc_put_att_text (nc_file_id, nc_srcgrdfrac_id, "units", 8, "unitless");
    ERR
    ncstat = nc_def_var (nc_file_id, "frac_b", NC_DOUBLE, 1, &nc_dstgrdsize_id, &nc_dstgrdfrac_id);
    ERR
    ncstat = nc_put_att_text (nc_file_id, nc_dstgrdfrac_id, "units", 8, "unitless");
    ERR

    // define mapping arrays
    ncstat = nc_def_var (nc_file_id, "col", NC_INT, 1, &nc_numlinks_id, &nc_srcadd_id);
    ERR
    ncstat = nc_def_var (nc_file_id, "row", NC_INT, 1, &nc_numlinks_id, &nc_dstadd_id);
    ERR
    ncstat = nc_def_var (nc_file_id, "S", NC_DOUBLE, 1, &nc_numlinks_id, &nc_rmpmatrix_id);
    ERR
    nc_dims2_id = new int[2];
    if (num_wts > 1)
    {
        nc_dims2_id[1] = nc_numwgts1_id;
        nc_dims2_id[0] = nc_numlinks_id;
    }
    ncstat = nc_def_var (nc_file_id, "S2", NC_DOUBLE, 2, nc_dims2_id, &nc_rmpmatrix2_id);
    ERR

    /* end definition stage */
    ncstat = nc_close (nc_file_id);
    ERR

    /* compute integer masks */
    if (direction == 1)
    {
        src_mask_int = new int [grid1_size];
        dst_mask_int = new int [grid2_size];
        for (int i = 0; i < grid1_size; i++)
        {
            if (grid1_mask[i])
                src_mask_int[i] = 1;
            else
                src_mask_int[i] = 0;
        }
        for (int i = 0; i < grid2_size; i++)
        {
            if (grid2_mask[i])
                dst_mask_int[i] = 1;
            else
                dst_mask_int[i] = 0;
        }
    }
    else
    {
        src_mask_int = new int [grid2_size];
        dst_mask_int = new int [grid1_size];
        for (int i = 0; i < grid2_size; i++)
        {
            if (grid2_mask[i])
                src_mask_int[i] = 1;
            else
                src_mask_int[i] = 0;
        }
        for (int i = 0; i < grid1_size; i++)
        {
            if (grid1_mask[i])
                dst_mask_int[i] = 1;
            else
                dst_mask_int[i] = 0;
        }
    }


    /* change units of lat/lon coordinates if input units different from radians */
    if (streqls(grid1_units, "degrees", 0, 7) && direction == 1)
    {
        for (int i = 0; i < grid1_size; i++)
        {
            grid1_center_lat[i] /= deg2rad;
            grid1_center_lon[i] /= deg2rad;
        }
        int num_corners1 = grid1_size * grid1_corners_max;
        for (int i = 0; i < num_corners1; i++)
        {
            grid1_corner_lat[i] /= deg2rad;
            grid1_corner_lon[i] /= deg2rad;
        }
    }
    if (streqls(grid2_units, "degrees", 0, 7) && direction == 1)
    {
        for (int i = 0; i < grid2_size; i++)
        {
            grid2_center_lat[i] /= deg2rad;
            grid2_center_lon[i] /= deg2rad;
        }
        int num_corners2 = grid2_size * grid2_corners_max;
        for (int i = 0; i < num_corners2; i++)
        {
            grid2_corner_lat[i] /= deg2rad;
            grid2_corner_lon[i] /= deg2rad;
        }
    }


    /* write mapping data */
    // write mask
    if (direction == 1)
    {
        itmp1 = nc_srcgrddims_id;
        itmp2 = nc_dstgrddims_id;
    }
    else
    {
        itmp2 = nc_srcgrddims_id;
        itmp1 = nc_dstgrddims_id;
    }

    ncstat = nc_put_var_int (nc_file_id, itmp1, grid1_dims);
    ERR
    ncstat = nc_put_var_int (nc_file_id, itmp2, grid2_dims);
    ERR

    ncstat = nc_put_var_int (nc_file_id, nc_srcgrdimask_id, src_mask_int);
    ERR
    ncstat = nc_put_var_int (nc_file_id, nc_dstgrdimask_id, dst_mask_int);
    ERR
    // delete mask_int
    delete [] src_mask_int;
    delete [] dst_mask_int;

    // write center|corner latitude|longitude
    // write grid1 center|corner lat|lon
    if (direction == 1)
    {
        itmp1 = nc_srcgrdcntrlat_id;
        itmp2 = nc_srcgrdcntrlon_id;
        itmp3 = nc_srcgrdcrnrlat_id;
        itmp4 = nc_srcgrdcrnrlon_id;
    }
    else
    {
        itmp1 = nc_dstgrdcntrlat_id;
        itmp2 = nc_dstgrdcntrlon_id;
        itmp3 = nc_dstgrdcrnrlat_id;
        itmp4 = nc_dstgrdcrnrlon_id;
    }

    ncstat = nc_put_var_double (nc_file_id, itmp1, grid1_center_lat);
    ERR
    ncstat = nc_put_var_double (nc_file_id, itmp2, grid1_center_lon);
    ERR
    ncstat = nc_put_var_double (nc_file_id, itmp3, grid1_corner_lat);
    ERR
    ncstat = nc_put_var_double (nc_file_id, itmp4, grid1_corner_lon);
    ERR
    // write grid2 center|corner lat|lon
    if (direction == 1)
    {
        itmp1 = nc_dstgrdcntrlat_id;
        itmp2 = nc_dstgrdcntrlon_id;
        itmp3 = nc_dstgrdcrnrlat_id;
        itmp4 = nc_dstgrdcrnrlon_id;
    }
    else
    {
        itmp1 = nc_srcgrdcntrlat_id;
        itmp2 = nc_srcgrdcntrlon_id;
        itmp3 = nc_srcgrdcrnrlat_id;
        itmp4 = nc_srcgrdcrnrlon_id;
    }

    ncstat = nc_put_var_double (nc_file_id, itmp1, grid2_center_lat);
    ERR
    ncstat = nc_put_var_double (nc_file_id, itmp2, grid2_center_lon);
    ERR
    ncstat = nc_put_var_double (nc_file_id, itmp3, grid2_corner_lat);
    ERR
    ncstat = nc_put_var_double (nc_file_id, itmp4, grid2_corner_lon);
    ERR

    // write area and frac
    if (direction == 1)
    {
        itmp1 = nc_srcgrdarea_id;
        itmp2 = nc_srcgrdfrac_id;
        itmp3 = nc_dstgrdarea_id;
        itmp4 = nc_dstgrdfrac_id;
    }
    else
    {
        itmp1 = nc_dstgrdarea_id;
        itmp2 = nc_dstgrdfrac_id;
        itmp3 = nc_srcgrdarea_id;
        itmp4 = nc_srcgrdfrac_id;
    }

    if (luse_grid1_area)
        ncstat = nc_put_var_double (nc_file_id, itmp1, grid1_area_in);
    else
        ncstat = nc_put_var_double (nc_file_id, itmp1, grid1_area);
    ERR

    if (luse_grid2_area)
        ncstat = nc_put_var_double (nc_file_id, itmp2, grid2_area_in);
    else
        ncstat = nc_put_var_double (nc_file_id, itmp2, grid2_area);
    ERR
    
    ncstat = nc_put_var_double (nc_file_id, itmp3, grid1_frac);
    ERR
    ncstat = nc_put_var_double (nc_file_id, itmp4, grid2_frac);
    ERR

    // write src|dst grid add, remap matrix of mappings
    ncstat = nc_put_var_int (nc_file_id, nc_srcadd_id, grid1_add_map);
    ERR
    ncstat = nc_put_var_int (nc_file_id, nc_dstadd_id, grid2_add_map);
    ERR

    if (num_wts == 1)
    {
        ncstat = nc_put_var_double (nc_file_id, nc_rmpmatrix_id, wts_map);
        ERR
    }
    else
    {
        // allocate weights arrays
        wts1 = new double [num_links_map];
        wts2 = new double [(num_wts-1) * num_links_map];
        int wtsdex = 0;
        int wts2dex = 0;
        for (int i = 0; i < num_links_map; i++)
        {
            wts1[i] = wts_map[wtsdex];
            for (int j = 1; j < num_wts; j++)
            {
                wts2[wts2dex++] = wts_map[wtsdex + j];
            }
            wtsdex += num_wts;
        }
        ncstat = nc_put_var_double (nc_file_id, nc_rmpmatrix_id, wts1);
        ERR
        ncstat = nc_put_var_double (nc_file_id, nc_rmpmatrix2_id, wts2);
        ERR

        // delete array wts1, wts2
        delete [] wts1;
        delete [] wts2;
    }

    /* close netCDF file */
    nc_close (nc_file_id);
    ERR
}

/** sorts address and weight arrays based on the 
  * destination address with the source address as a secondary 
  * sorting criterion. the method is a standard heap sort.
  **/
void sort_add(int *add1, int *add2, double *weights, int num_links, int num_weights)
{
    log("sort add");
    // local variables
    int add1_tmp, add2_tmp;     // temp for address during swap
    int nwgt;
    int lvl, final_lvl;         // level indexes for heap sort levels
    int chk_lvl1, chk_lvl2, max_lvl;

    double wgttmp[num_weights]; // temp for holding wts during swap

    /** start at the lowest level (N/2) of the tree and sift lower
      * values to the bottom of the tree, prompting the largest numbers
      **/
    int wgtdex = num_links/2;       // initial header used for weights[]
    for (lvl = num_links/2; lvl > -1; lvl --)
    {
        final_lvl = lvl;
        add1_tmp = add1[lvl];
        add2_tmp = add2[lvl];
        for (int i = 0; i < num_weights; i++)
        {
            wgttmp[i] = weights[wgtdex + i];
        }
        wgtdex -= num_weights;      // calculate wgtdex for next

        // find the largest of the two daughters
        while (true)
        {
            chk_lvl1 = 2 * final_lvl;           // left child
            chk_lvl2 = 2 * final_lvl + 1;       // right child

            if ((add1[chk_lvl1] > add1[chk_lvl2]) ||
                    ((add1[chk_lvl1] == add1[chk_lvl2]) &&
                     (add2[chk_lvl1] > add2[chk_lvl2])))
            {
                max_lvl = chk_lvl1;
            }
            else
            {
                max_lvl = chk_lvl2;
            }

            // if the parent is greater than both daughters, the correct level has been found
            if ((add1_tmp > add1[max_lvl]) ||
                    ((add1_tmp == add1[max_lvl]) &&
                     (add2_tmp > add2[max_lvl])))
            {
                add1[final_lvl] = add1_tmp;
                add2[final_lvl] = add2_tmp;
                for (int i = 0; i < num_weights; i++)
                {
                    weights[final_lvl * num_weights + i] = wgttmp[i];
                }
                break;
            }

            // otherwise, promte the largest daughter and push down one level in the tree. if haven't reached the end of the tree, repeat the process. otherwise store last values and exit the loop
            else
            {
                add1[final_lvl] = add1[max_lvl];
                add2[final_lvl] = add2[max_lvl];
                for (int i = 0; i < num_weights; i++)
                {
                    weights[final_lvl * num_weights + i] = weights[max_lvl * num_weights + i];
                }

                final_lvl = max_lvl;
                if (final_lvl * 2 > num_links)
                {
                    add1[final_lvl] = add1_tmp;
                    add2[final_lvl] = add2_tmp;
                    for (int i = 0; i < num_weights; i++)
                    {
                        weights[final_lvl * num_weights + i] = wgttmp[i];
                    }
                    break;
                }
            }
        }

        /** now that the heap has been sorted, strip off the top (largest)
         *  value and promote the values below
         **/
        int strip_wgtdex = 0;
        for (lvl = num_links; lvl > 2; lvl --)
        {
            // move the top value and insert it into the correct place
            add1_tmp = add1[lvl];
            add1[lvl] = add1[0];

            add2_tmp = add2[lvl];
            add2[lvl] = add2[0];
            strip_wgtdex = lvl * num_weights;   // header index
            for (int i = 0; i < num_weights; i++)
            {
                wgttmp[i] = weights[strip_wgtdex + i];
                weights[strip_wgtdex + i] = weights[i];
            }

            // as above this loop sifts the tmp values down until proper level is reached
            final_lvl = 0;

            while (true)
            {
                // find the largest of the two daughters
                chk_lvl1 = 2 * final_lvl;
                chk_lvl2 = 2 * final_lvl + 1;
                if (chk_lvl2 >= lvl)
                    chk_lvl2 = chk_lvl1;

                if ((add1[chk_lvl1] > add1[chk_lvl2]) ||
                        ((add1[chk_lvl1] == add1[chk_lvl2]) &&
                         (add2[chk_lvl1] > add2[chk_lvl2])))
                {
                    max_lvl = chk_lvl1;
                }
                else
                {
                    max_lvl = chk_lvl2;
                }

                // if the parent is greater than both daughters, the correct level had been found
                if ((add1_tmp > add1[max_lvl]) ||
                        ((add1_tmp == add1[max_lvl]) &&
                         (add2_tmp > add2[max_lvl])))
                {
                    add1[final_lvl] = add1_tmp;
                    add2[final_lvl] = add2_tmp;
                    for (int i = 0; i < num_weights; i++)
                    {
                        weights[final_lvl * num_weights + i] = wgttmp[i];
                    }
                    break;
                }

                // otherwise, promote the largest daughter and push down one level in the tree. if haven't reached the end of the tree, repeat the process. otherwise store last values and exit the loop
                else
                {
                    add1[final_lvl] = add1[max_lvl];
                    add2[final_lvl] = add2[max_lvl];
                    for (int i = 0; i < num_weights; i++)
                    {
                        weights[final_lvl * num_weights + i] = weights[max_lvl * num_weights + i];
                    }

                    final_lvl = max_lvl;
                    if (final_lvl * 2 > lvl)
                    {
                        add1[final_lvl] = add1_tmp;
                        add2[final_lvl] = add2_tmp;
                        for (int i = 0; i < num_weights; i++)
                        {
                            weights[final_lvl * num_weights + i] = wgttmp[i];
                        }
                        break;
                    }   
                }
            }
        }
    }

    // swap the last two entries
    add1_tmp = add1[1];
    add1[1] = add1[0];
    add1[0] = add1_tmp;

    add2_tmp = add2[1];
    add2[1] = add2[0];
    add2[0] = add2_tmp;

    for (int i = 0; i < num_weights; i++)
    {
        wgttmp[i] = weights[num_weights + i];
        weights[num_weights + i] = weights[i];
        weights[i] = wgttmp[i];
    }
}
