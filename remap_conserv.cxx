#include "remap_conserv.h"
#include "debug.h"

/** this routine traces the perimeters of every grid cell on each
 *  grid checking for intersections with the other grid and computing
 *  line integrals for each subsegment
 **/
void remap_conserv()
{
    // local variables
    
    int grid1_add, grid2_add;   // linear index

    double norm_factor;         // factor for normalizing wts

    double *grid2_centroid_lat; // centroid coords on each grid
    double *grid2_centroid_lon;
    double *grid1_centroid_lat;
    double *grid1_centroid_lon;
    double *wts;

    /* initial centroid arrays */
    grid1_centroid_lat = new double [grid1_size];
    grid1_centroid_lon = new double [grid1_size];
    grid2_centroid_lat = new double [grid2_size];
    grid2_centroid_lon = new double [grid2_size];

    memset(grid1_centroid_lat, ZERO, sizeof(double) * grid1_size);
    memset(grid1_centroid_lon, ZERO, sizeof(double) * grid1_size);
    memset(grid2_centroid_lat, ZERO, sizeof(double) * grid2_size);
    memset(grid2_centroid_lon, ZERO, sizeof(double) * grid2_size);

    /* initial new of search arrays */
    srch_add = new int [SRCH_SIZE];
    srch_corner_lat = new double [SRCH_SIZE * SRCH_CORNER_MAX];
    srch_corner_lon = new double [SRCH_SIZE * SRCH_CORNER_MAX];
    memset(srch_add, -1, sizeof(int) * SRCH_SIZE);
    memset(srch_corner_lat, 0, sizeof(double) * SRCH_SIZE * SRCH_CORNER_MAX);
    memset(srch_corner_lon, 0, sizeof(double) * SRCH_SIZE * SRCH_CORNER_MAX);

    wts = new double[num_wts];

    /* integrate around each cell on grid1 */
    int error = 0;
    error = conserv_sweep(INTEGRATE_AROUND_SRC_GRID, grid1_centroid_lat, grid1_centroid_lon, grid2_centroid_lat, grid2_centroid_lon);
    if (error == 1)
    {
        printf("\n");
        return;
    }

    /* integrate around each cell on grid2 */
    error = conserv_sweep(INTEGRATE_AROUND_DST_GRID, grid1_centroid_lat, grid1_centroid_lon, grid2_centroid_lat, grid2_centroid_lon);
    if (error == 1)
    {
        printf("\n");
        return;
    }
    
    /* correct for situations where N/S pole not explicitly included in grid (i.e. as a grid corner point). if pole is missing from only one grid, need to correct only the area and centroid of that grid. if missing from both, do complete weight calculation */
    
    /* finish centroid computation */
    for (int n = 0; n < grid1_size; n++)
    {
        if (nonzero(grid1_area[n]))
        {
            grid1_centroid_lat[n] /= grid1_area[n];
            grid1_centroid_lon[n] /= grid1_area[n];
        }
    }

    for (int n = 0; n < grid2_size; n++)
    {
        if (nonzero(grid2_area[n]))
        {
            grid2_centroid_lat[n] /= grid2_area[n];
            grid2_centroid_lon[n] /= grid2_area[n];
        }
    }

    /* include centroids in weights and normalize using destination area if required */
    for (int n = 0; n < num_links_map; n++)
    {
        grid1_add = grid1_add_map[n];
        grid2_add = grid2_add_map[n];
        for (int nwgt = 0; nwgt < num_wts; nwgt ++)
        {
            wts[nwgt] = wts_map[n*num_wts + nwgt];
        }

        switch (norm_opt)
        {
            case NORM_OPT_DESTAREA:
                if (nonzero(grid2_area[grid2_add]))
                {
                    if (luse_grid2_area)
                        norm_factor = ONE / grid2_area_in[grid2_add];
                    else
                        norm_factor = ONE / grid2_area[grid2_add];
                }
                else
                    norm_factor = ZERO;
                break;
            case NORM_OPT_FRACAREA:
                if (nonzero(grid2_area[grid2_add]))
                {
                    if (luse_grid2_area)
                        norm_factor = grid2_area[grid2_add] /
                            (grid2_frac[grid2_add] * grid2_area_in[grid2_add]);
                    else
                        norm_factor = ONE / grid2_frac[grid2_add];
                }
                else
                    norm_factor = ZERO;
                break;
            case NORM_OPT_NONE:
                norm_factor = ONE;
                break;
            default:
                norm_factor = ONE;
        }

        wts_map[n*num_wts] = wts[0] * norm_factor;
        wts_map[n*num_wts + 1] = (wts[1] - wts[0]*grid1_centroid_lat[grid1_add]) * norm_factor;
        wts_map[n*num_wts + 2] = (wts[2] - wts[0]*grid1_centroid_lon[grid1_add]) * norm_factor;
    }    

    printf("Total number of links = %d\n", num_links_map);
    for (int n = 0; n < grid1_size; n++)
        if (nonzero(grid1_area[n]))
            grid1_frac[n] /= grid1_area[n];

    for (int n = 0; n < grid2_size; n++)
        if (nonzero(grid2_area[n]))
            grid2_frac[n] /= grid2_area[n];
    

    /* perform some error checking on final weights */
#if _CHECK_AREA_CENTROID_
    // check grid1 area and centroid lat
    for (int n = 0; n < grid1_size; n++)
    {
        if (grid1_area[n] < -0.01)
            printf("Grid 1 area error @ %d\t valued %3.6f\n", n, grid1_area[n]);
        if (grid1_centroid_lat[n] < -PIH - 0.01 ||
                grid1_centroid_lat[n] > PIH + 0.01)
            printf("Grid 1 centroid lat error @ %d\t valued %3.6f\n", n, grid1_centroid_lat[n]);
        grid1_centroid_lat[n] = ZERO;
        grid1_centroid_lon[n] = ZERO;
    }
    // check grid2 area and centroid lat
    for (int n = 0; n < grid2_size; n++)
    {
        if (grid2_area[n] < -0.01)
            printf("Grid 2 area error @ %d\t valued %3.6f\n", n, grid2_area[n]);
        if (grid2_centroid_lat[n] < -PIH - 0.01 ||
                grid2_centroid_lat[n] > PIH + 0.01)
            printf("Grid 2 centroid lat error @ %d\t valued %3.6f\n", n, grid2_centroid_lat[n]);
        grid2_centroid_lat[n] = ZERO;
        grid2_centroid_lon[n] = ZERO;
    }
    // check weights
    grid2_add = 0;
    for (int n = 0; n < num_links_map; n++)
    {
        grid1_add = grid1_add_map[n];
        grid2_add = grid2_add_map[n];

        if (wts_map[num_wts * n] < -0.01)
            printf("Map Weight < 0 @grid1_add=%d\tgrid2_add=%d\t valued %3.6f\n", grid1_add, grid2_add, wts_map[num_wts*2 * n]);
        if (norm_opt != NORM_OPT_NONE && wts_map[num_wts * n] > 1.01)
            printf("Map weight > 1 @grid1_add=%d\tgrid2_add=%d\t valued %3.6f\n", grid1_add, grid2_add, wts_map[num_wts*2 * n]);
        grid2_centroid_lat [grid2_add] += wts_map[num_wts * n];
    }

    for (int n = 0; n < grid2_size; n++)
    {
        switch (norm_opt)
        {
            case NORM_OPT_DESTAREA:
                norm_factor = grid2_frac[grid2_add];
                break;
            case NORM_OPT_FRACAREA:
                norm_factor = ONE;
            case NORM_OPT_NONE:
                if (luse_grid2_area)
                    norm_factor = grid2_area_in [grid2_add];
                else
                    norm_factor = grid2_area [grid2_add];
                break;
            default:
                break;
        }
        if (abs(grid2_centroid_lat[grid2_add] - norm_factor ) > 0.01)
            printf("Error: sum of wts for map @%d\t %3.6f\t vs %3.6f\n", grid2_add, grid2_centroid_lat[grid2_add], norm_factor);
    }
#endif
    delete [] grid1_centroid_lat;
    delete [] grid1_centroid_lon;
    delete [] grid2_centroid_lat;
    delete [] grid2_centroid_lon;
    delete [] wts;
    /* delete search arrays */
    //delete [] srch_add;
    //delete [] srch_corner_lat;
    //delete [] srch_corner_lon;
}

/* integrate around each cell on one grid */
int conserv_sweep(int choice, double *grid1_centroid_lat, double *grid1_centroid_lon, 
        double *grid2_centroid_lat, double *grid2_centroid_lon)
{
    // local variables of orignal remap_conserv
    int max_subseg = 1000;      // max number of subsegments per segment to prevent infinite loop
    int num_subseg;             // number of subsegments
    int n, nwgt;                // generic counter
    int corner, next_corn;      // corner of cell that segment starts from and ends on
    int min_add, max_add;               // store min/max index

    bool lcoinc;                // flag for coincident segments
    bool lrevers;               // flag for reversing direction of segment
    bool lbegin;                // flag for first integration of a segment
    bool *srch_mask;            // mask for restricting searches

    double intrsct_lat;         // lat of next intersect
    double intrsct_lon;         // lon of next intersect
    double beglat, endlat;      // endpoints of current seg
    double beglon, endlon;
    double *begseg;             // begin lat/lon for full segment
    double *conserv_weights;            // temp for weights

    // init variables
    intrsct_lat = ZERO;
    intrsct_lon = ZERO;
    begseg = new double[2];
    conserv_weights = new double[num_wts * 2];

    // print sweep info
    printf("grid%d sweep\n", choice);

    // pointers to corresponding array
    int src_grid_add, dst_grid_add;     // src/dst grid index
    unsigned int src_grid_size, dst_grid_size;   // src/dst grid size
    unsigned int src_grid_corners_max, dst_grid_corners_max;     // src/dst grid corner max
    unsigned int *src_grid_corners, *dst_grid_corners;
    int * src_grid_bin_addr, *dst_grid_bin_addr;        // src/dst srch bin lat range
    double *src_grid_corner_lat, *dst_grid_corner_lat;  // pointer to extern array
    double *src_grid_corner_lon, *dst_grid_corner_lon;
    double *src_grid_center_lat, *dst_grid_center_lat;
    double *src_grid_center_lon, *dst_grid_center_lon;
    double *src_grid_centroid_lat, *dst_grid_centroid_lat;
    double *src_grid_centroid_lon, *dst_grid_centroid_lon;
    bool *src_grid_mask, *dst_grid_mask;
    double *src_grid_area, *dst_grid_area;
    double *src_grid_frac, *dst_grid_frac;
    double *src_grid_bound_box, *dst_grid_bound_box;

    // define variables according to choice(=1 =2)
    if (choice == INTEGRATE_AROUND_SRC_GRID)
    {
        src_grid_size = grid1_size;
        dst_grid_size = grid2_size;
        src_grid_corners_max = grid1_corners_max;
        dst_grid_corners_max = grid2_corners_max;
        src_grid_corner_lat = grid1_corner_lat;
        dst_grid_corner_lat = grid2_corner_lat;
        src_grid_bin_addr = bin_addr1;
        dst_grid_bin_addr = bin_addr2;
        src_grid_corner_lon = grid1_corner_lon;
        dst_grid_corner_lon = grid2_corner_lon;
        src_grid_center_lat = grid1_center_lat;
        dst_grid_center_lat = grid2_center_lat;
        src_grid_center_lon = grid1_center_lon;
        dst_grid_center_lon = grid2_center_lon;
        src_grid_centroid_lat = grid1_centroid_lat;
        dst_grid_centroid_lat = grid2_centroid_lat;
        src_grid_centroid_lon = grid1_centroid_lon;
        dst_grid_centroid_lon = grid2_centroid_lon;
        src_grid_mask = grid1_mask;
        dst_grid_mask = grid2_mask;
        src_grid_area = grid1_area;
        dst_grid_area = grid2_area;
        src_grid_frac = grid1_frac;
        dst_grid_frac = grid2_frac;
        src_grid_bound_box = grid1_bound_box;
        dst_grid_bound_box = grid2_bound_box;
        src_grid_corners = grid1_corners;
        dst_grid_corners = grid2_corners;
    }
    else if (choice == INTEGRATE_AROUND_DST_GRID)
    {
        src_grid_size = grid2_size;
        dst_grid_size = grid1_size;
        src_grid_corners_max = grid2_corners_max;
        dst_grid_corners_max = grid1_corners_max;
        src_grid_bin_addr = bin_addr2;
        dst_grid_bin_addr = bin_addr1;
        src_grid_corner_lat = grid2_corner_lat;
        dst_grid_corner_lat = grid1_corner_lat;
        src_grid_corner_lon = grid2_corner_lon;
        dst_grid_corner_lon = grid1_corner_lon;
        src_grid_center_lat = grid2_center_lat;
        dst_grid_center_lat = grid1_center_lat;
        src_grid_center_lon = grid2_center_lon;
        dst_grid_center_lon = grid1_center_lon;
        src_grid_centroid_lat = grid2_centroid_lat;
        dst_grid_centroid_lat = grid1_centroid_lat;
        src_grid_centroid_lon = grid2_centroid_lon;
        dst_grid_centroid_lon = grid1_centroid_lon;
        src_grid_mask = grid2_mask;
        dst_grid_mask = grid1_mask;
        src_grid_area = grid2_area;
        dst_grid_area = grid1_area;
        src_grid_frac = grid2_frac;
        dst_grid_frac = grid1_frac;
        src_grid_bound_box = grid2_bound_box;
        dst_grid_bound_box = grid1_bound_box;
        src_grid_corners = grid2_corners;
        dst_grid_corners = grid1_corners;
    }
    else
    {
        printf("No such sweep choice\n");
        return 1;
    }
    
    // allocate srch_mask
    srch_mask = new bool [dst_grid_size];
    memset(srch_mask, false, sizeof(bool) * dst_grid_size);

    // integrate around on on each src_grid cell
    for (src_grid_add = 0; src_grid_add < src_grid_size; src_grid_add ++)
    {
        // restrict searches first using search bins
        min_add = dst_grid_size;
        max_add = -1;
        for (n = 0; n < num_srch_bins; n++)
        {
            if (src_grid_add >= src_grid_bin_addr[n*2] && 
                    src_grid_add <= src_grid_bin_addr[n*2 + 1])
            {
                min_add = MIN(min_add, dst_grid_bin_addr[n*2]);
                max_add = MAX(max_add, dst_grid_bin_addr[n*2 + 1]);
            }
        }
        
        // further restrict searches using bounding boxes
        num_srch_cells = 0;     // extern global var
        //printf("from %d to %d\n", min_add, max_add);
        for (dst_grid_add = min_add; dst_grid_add <= max_add; dst_grid_add ++)
        {
            srch_mask[dst_grid_add] = true;
            srch_mask[dst_grid_add] &= 
                (dst_grid_bound_box[dst_grid_add * BOUNDBOX_SIZE] <= src_grid_bound_box[src_grid_add * BOUNDBOX_SIZE + 1]);
            srch_mask[dst_grid_add] &=
                (dst_grid_bound_box[dst_grid_add * BOUNDBOX_SIZE + 1] >= src_grid_bound_box[src_grid_add * BOUNDBOX_SIZE]);
            srch_mask[dst_grid_add] &=
                ( dst_grid_bound_box[dst_grid_add * BOUNDBOX_SIZE + 2] <= src_grid_bound_box[src_grid_add *BOUNDBOX_SIZE + 3]);
                //le( dst_grid_bound_box[dst_grid_add * BOUNDBOX_SIZE + 2], src_grid_bound_box[src_grid_add *BOUNDBOX_SIZE + 3]);
            srch_mask[dst_grid_add] &=
                ( dst_grid_bound_box[dst_grid_add *BOUNDBOX_SIZE + 3] >= src_grid_bound_box[src_grid_add * BOUNDBOX_SIZE + 2]);
                //ge( dst_grid_bound_box[dst_grid_add *BOUNDBOX_SIZE + 3], src_grid_bound_box[src_grid_add * BOUNDBOX_SIZE + 2]);
            if (_3D_BOUND_BOX_)         // using 3D bound box option
            {
                srch_mask[dst_grid_add] &=
                    (dst_grid_bound_box[dst_grid_add * BOUNDBOX_SIZE] + 4 <= src_grid_bound_box[src_grid_add * BOUNDBOX_SIZE + 5]);
                srch_mask[dst_grid_add] &=
                    (dst_grid_bound_box[dst_grid_add * BOUNDBOX_SIZE + 5] >= src_grid_bound_box[src_grid_add * BOUNDBOX_SIZE + 4]);
            }

            if (srch_mask[dst_grid_add])
                num_srch_cells ++;
        }
#if _DEBUG_SRCH_
        printf("%5d%5d%5d\n", min_add+1, max_add+1, num_srch_cells);
#endif
        // create search arrays
        //printf("num_srch_cells of grid %d: %d\n", src_grid_add, num_srch_cells);
        /*
           * To use static array instead of new/delete each time
        srch_add = new int [num_srch_cells];
        srch_corner_lat = new double [dst_grid_corners_max * num_srch_cells];
        srch_corner_lon = new double [dst_grid_corners_max * num_srch_cells];
        */
        n = 0;
        for (dst_grid_add = min_add; dst_grid_add <= max_add; dst_grid_add ++)
        {
            if (srch_mask[dst_grid_add])
            {
                srch_add[n] = dst_grid_add;
                for (corner = 0; corner < dst_grid_corners[dst_grid_add]; corner ++)
                {
                    srch_corner_lat[n * dst_grid_corners_max + corner] = 
                        dst_grid_corner_lat[dst_grid_add * dst_grid_corners_max + corner];
                    srch_corner_lon[n * dst_grid_corners_max + corner] = 
                        dst_grid_corner_lon[dst_grid_add * dst_grid_corners_max + corner];
                }
                n++;
            }
        }

        // integrate around this cell
        int num_corner = src_grid_corners[src_grid_add];
        for (corner = 0; corner < num_corner; corner ++)
        {
            next_corn = (corner + 1) %  num_corner;

            // define endpoints of the current segment
            int corndex = src_grid_add * src_grid_corners_max;
            beglat = src_grid_corner_lat [corndex + corner];
            beglon = src_grid_corner_lon [corndex + corner];
            endlat = src_grid_corner_lat [corndex + next_corn];
            endlon = src_grid_corner_lon [corndex + next_corn];
            lrevers = false;

            // to ensure exact path taken during both sweeps, always integrate segments in the same direction (SW TO NE)
            if ((endlat < beglat) || 
                    (endlat == beglat && endlon < beglon ))
            {
                beglat = src_grid_corner_lat [corndex + next_corn];
                beglon = src_grid_corner_lon [corndex + next_corn];
                endlat = src_grid_corner_lat [corndex + corner];
                endlon = src_grid_corner_lon [corndex + corner];
                lrevers = true;
            }

            begseg[0] = beglat;
            begseg[1] = beglon;
            lbegin = true;
            num_subseg = 0;

            // if this is a constant-longitude segment, skip the rest since the line integral contribution will be zero
            if (nonzero(endlon - beglon))
            {
                // integrate along this segment, detecting intersections and computing the line integral for each sub-segment
                while (nonzero(beglat - endlat) || nonzero(beglon - endlon))
                {
                    // prevent infinite loops if integration gets stuck near cell or threshold boundary
                    num_subseg ++;
                    if (num_subseg > max_subseg)
                    {
                        printf("Integrate stalled: num_subseg exceeded limit @ src_grid %d\n", src_grid_add);
                        return 1;
                    }

                    // find next intersection of this segment with a grid line on dst grid
                    intersection(dst_grid_add, intrsct_lat, intrsct_lon, lcoinc, beglat, beglon, endlat, endlon, begseg, lbegin, lrevers);
                    lbegin = false;

                    // compute line integral for this subsegment
                    if (dst_grid_add != -1)
                    {
                        line_integral(conserv_weights, num_wts, beglon, intrsct_lon, beglat, intrsct_lat, 
                                src_grid_center_lat[src_grid_add],
                                src_grid_center_lon[src_grid_add],
                                dst_grid_center_lat[dst_grid_add],
                                dst_grid_center_lon[dst_grid_add]);
                    }
                    else
                    {
                        line_integral(conserv_weights, num_wts, beglon, intrsct_lon, beglat, intrsct_lat,
                                src_grid_center_lat[src_grid_add],
                                src_grid_center_lon[src_grid_add],
                                src_grid_center_lat[src_grid_add],
                                src_grid_center_lon[src_grid_add]);
                    }

                    // if integrating in reverse order, change sign of weights
                    if (lrevers)
                    {
                        for (int i = 0; i < num_wts * 2; i++)
                        {
                            conserv_weights[i] = - conserv_weights[i];
                        }
                    }

                    // store the appropriate addresses and weights. also add contributions to cell areas and centroids
                    if (choice == INTEGRATE_AROUND_SRC_GRID)
                    {
                        if (dst_grid_add != -1)
                        {
                            if (grid1_mask[src_grid_add])
                            {
                                //store_link_cnsrv(src_grid_add, dst_grid_add, conserv_weights, num_wts);
                                grid1_frac[src_grid_add] += conserv_weights[0];
                                grid2_frac[dst_grid_add] += conserv_weights[num_wts];
                            }
                        }
                    }
                    else if (choice == INTEGRATE_AROUND_DST_GRID)
                    {
                        if (dst_grid_add != -1 && !lcoinc)
                        {
                            if (grid1_mask[dst_grid_add])
                            {
                                //store_link_cnsrv(dst_grid_add, src_grid_add, conserv_weights, num_wts);
                                grid1_frac[dst_grid_add] += conserv_weights[0];
                                grid2_frac[src_grid_add] += conserv_weights[num_wts];
                            }
                        }
                    }

                    src_grid_area[src_grid_add] += conserv_weights[0];
                    src_grid_centroid_lat[src_grid_add] += conserv_weights[1];
                    src_grid_centroid_lon[src_grid_add] += conserv_weights[2];

                    // reset beglat and beglon for next subsegment
                    beglat = intrsct_lat;
                    beglon = intrsct_lon;
                }
            }
        }
    }
    //printf("before delete weights\n");
    delete [] srch_mask;
    delete [] begseg;
    delete [] conserv_weights;

    return 0;
}   
