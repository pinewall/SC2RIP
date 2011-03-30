#include "remap_write.h"

#include "kinds.h"
#include "constants.h"
#include "grids.h"
#include "remap_vars.h"
#include "nc_error.h"
#include "utils.h"
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

int nc_dims2_id[2];             // netCDF ids for 2d array dims


// call correct output routine based on output format choice (scrip or csm)
void write_remap(char *map_name, char *interp_file, int output_opt)
{
    printf("Write_Remap");

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
    sprintf(history, "Created: %s\n", cdate);


    // sort address and weight arrays
    sort_add(grid2_add_map, grid1_add_map, wts_map);

    // call appropriate output routine
    switch (output_opt)
    {
        case SCRIP_CONVENTION:
            write_remap_scrip(map_name, interp_file, 1);
            break;
        case CSM_CONVENTION:
            write_remap_csm  (map_name, interp_file, 1);
            break;
        default:
            log("Unknown output file convention");
            return;
    }

}

// writes remap data to a netCDF file using SCRIP conventions
void write_remap_scrip(char *map_name, char *interp_file, int direction)
{
    log("--Using SCRIP convention");
}

// writes remap data to a netCDF file using NCAR-CSM conventions
void write_remap_csm(char *map_name, char *interp_file, int direction)
{
    log("--Using NCAR-CSM convention");
}

/** sorts address and weight arrays based on the 
  * destination address with the source address as a secondary 
  * sorting criterion. the method is a standard heap sort.
  **/
void sort_add(int *add1, int *add2, double *weights)
{
    log("sort_add");
}
