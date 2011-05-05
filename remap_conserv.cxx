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

    /* integrate around each cell on grid1 */
    // allocate srch_mask
    srch_mask = new bool [grid2_size];
    memset(srch_mask, false, sizeof(bool) * grid2_size);

    // integrate around on on each src_grid cell
    for (grid1_add = 0; grid1_add < grid1_size; grid1_add ++)
    {
#if _DEBUG_WEIGHTS_3_
        printf("index: %d\n", grid1_add + 1);
#endif
        // restrict searches first using search bins
        min_add = grid2_size;
        max_add = -1;
        for (n = 0; n < num_srch_bins; n++)
        {
            if (grid1_add >= bin_addr1[n*2] && 
                    grid1_add <= bin_addr1[n*2 + 1])
            {
                min_add = MIN(min_add, bin_addr2[n*2]);
                max_add = MAX(max_add, bin_addr2[n*2 + 1]);
            }
        }
        
        // further restrict searches using bounding boxes
        num_srch_cells = 0;     // extern global var
        //printf("from %d to %d\n", min_add, max_add);
        for (grid2_add = min_add; grid2_add <= max_add; grid2_add ++)
        {
            srch_mask[grid2_add] = true;
            srch_mask[grid2_add] &= 
                (grid2_bound_box[grid2_add * BOUNDBOX_SIZE] <= grid1_bound_box[grid1_add * BOUNDBOX_SIZE + 1]);
            srch_mask[grid2_add] &=
                (grid2_bound_box[grid2_add * BOUNDBOX_SIZE + 1] >= grid1_bound_box[grid1_add * BOUNDBOX_SIZE]);
            srch_mask[grid2_add] &=
                ( grid2_bound_box[grid2_add * BOUNDBOX_SIZE + 2] <= grid1_bound_box[grid1_add *BOUNDBOX_SIZE + 3]);
                //le( grid2_bound_box[grid2_add * BOUNDBOX_SIZE + 2], grid1_bound_box[grid1_add *BOUNDBOX_SIZE + 3]);
            srch_mask[grid2_add] &=
                ( grid2_bound_box[grid2_add *BOUNDBOX_SIZE + 3] >= grid1_bound_box[grid1_add * BOUNDBOX_SIZE + 2]);
                //ge( grid2_bound_box[grid2_add *BOUNDBOX_SIZE + 3], grid1_bound_box[grid1_add * BOUNDBOX_SIZE + 2]);
            if (_3D_BOUND_BOX_)         // using 3D bound box option
            {
                srch_mask[grid2_add] &=
                    (grid2_bound_box[grid2_add * BOUNDBOX_SIZE] + 4 <= grid1_bound_box[grid1_add * BOUNDBOX_SIZE + 5]);
                srch_mask[grid2_add] &=
                    (grid2_bound_box[grid2_add * BOUNDBOX_SIZE + 5] >= grid1_bound_box[grid1_add * BOUNDBOX_SIZE + 4]);
            }

            if (srch_mask[grid2_add])
                num_srch_cells ++;
        }
        // create search arrays
        n = 0;
        for (grid2_add = min_add; grid2_add <= max_add; grid2_add ++)
        {
            if (srch_mask[grid2_add])
            {
                srch_add[n] = grid2_add;
                for (corner = 0; corner < grid2_corners[grid2_add]; corner ++)
                {
                    srch_corner_lat[n * grid2_corners_max + corner] = 
                        grid2_corner_lat[grid2_add * grid2_corners_max + corner];
                    srch_corner_lon[n * grid2_corners_max + corner] = 
                        grid2_corner_lon[grid2_add * grid2_corners_max + corner];
                }
                n++;
            }
        }

        // integrate around this cell
        int num_corner = grid1_corners[grid1_add];
        for (corner = 0; corner < num_corner; corner ++)
        {
            next_corn = (corner + 1) %  num_corner;

            // define endpoints of the current segment
            int corndex = grid1_add * grid1_corners_max;
            beglat = grid1_corner_lat [corndex + corner];
            beglon = grid1_corner_lon [corndex + corner];
            endlat = grid1_corner_lat [corndex + next_corn];
            endlon = grid1_corner_lon [corndex + next_corn];
            lrevers = false;

            // to ensure exact path taken during both sweeps, always integrate segments in the same direction (SW TO NE)
            if ((endlat < beglat) || 
                    (endlat == beglat && endlon < beglon ))
            {
                beglat = grid1_corner_lat [corndex + next_corn];
                beglon = grid1_corner_lon [corndex + next_corn];
                endlat = grid1_corner_lat [corndex + corner];
                endlon = grid1_corner_lon [corndex + corner];
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
                        printf("Integrate stalled: num_subseg exceeded limit @ src_grid %d\n", grid1_add);
                    }

                    // find next intersection of this segment with a grid line on dst grid

#if _DEBUG_INPUT_OUTPUT_
                    // input compare for debug
                    printf("%6d:(%3.3f %3.3f) (%3.3f %3.3f %3.3f %3.3f) (%3.3f %3.3f)",
                            grid2_add + 1, intrsct_lat, intrsct_lon,
                            beglat, beglon, endlat, endlon, begseg[0], begseg[1]);
                    if (lrevers)
                        printf(" true\n");
                    else
                        printf(" false\n");
#endif
                    intersection(grid2_add, intrsct_lat, intrsct_lon, lcoinc, beglat, beglon, endlat, endlon, begseg, lbegin, lrevers);
                    lbegin = false;
                    
                    // compute line integral for this subsegment
                    if (grid2_add != -1)
                    {
                        line_integral(conserv_weights, num_wts, beglon, intrsct_lon, beglat, intrsct_lat, 
                                grid1_center_lat[grid1_add],
                                grid1_center_lon[grid1_add],
                                grid2_center_lat[grid2_add],
                                grid2_center_lon[grid2_add]);
                    }
                    else
                    {
                        //printf("unhappy\n");
                        line_integral(conserv_weights, num_wts, beglon, intrsct_lon, beglat, intrsct_lat,
                                grid1_center_lat[grid1_add],
                                grid1_center_lon[grid1_add],
                                grid1_center_lat[grid1_add],
                                grid1_center_lon[grid1_add]);
                    }
                    //printf("%3.6f %3.6f vs %3.6f %3.6f\n", intrsct_lat, intrsct_lon, beglat, beglon);
                    // if integrating in reverse order, change sign of weights
                    if (lrevers)
                    {
                        for (int i = 0; i < num_wts * 2; i++)
                        {
                            conserv_weights[i] = - conserv_weights[i];
                        }
                    }
                    
                    // store the appropriate addresses and weights. also add contributions to cell areas and centroids
                    if (grid2_add != -1)
                    {
                        if (grid1_mask[grid1_add])
                        {
#if _DEBUG_STORE_LINK_
                            printf("%6d--%6d: %3.8f %3.8f %3.8f %3.8f %3.8f %3.8f\n", grid1_add, grid2_add, 
                                    conserv_weights[0], conserv_weights[1], conserv_weights[2],
                                    conserv_weights[num_wts], conserv_weights[num_wts+1], conserv_weights[num_wts+2]);
#endif
                            store_link_cnsrv(grid1_add, grid2_add, conserv_weights, 3);
                            grid1_frac[grid1_add] += conserv_weights[0];
                            grid2_frac[grid2_add] += conserv_weights[num_wts];
                        }
                    }

#if _DEBUG_LINKS_
                    printf("%6d--%6d: %3.8f\t%3.8f\t%3.8f\n", grid1_add+1, grid2_add+1, 
                            conserv_weights[0], conserv_weights[1], conserv_weights[2]);
#endif
                    grid1_area[grid1_add] += conserv_weights[0];
                    grid1_centroid_lat[grid1_add] += conserv_weights[1];
                    grid1_centroid_lon[grid1_add] += conserv_weights[2];

                    // reset beglat and beglon for next subsegment
                    beglat = intrsct_lat;
                    beglon = intrsct_lon;
                }
            }
        }
    }
    delete [] srch_mask;

    /* integrate around each cell on grid2 */
    // allocate srch_mask
    srch_mask = new bool [grid1_size];
    memset(srch_mask, false, sizeof(bool) * grid1_size);

    // integrate around on on each grid2 cell
    for (grid2_add = 0; grid2_add < grid2_size; grid2_add ++)
    {
#if _DEBUG_STORE_LINK_
        //printf("index: %d\n", grid2_add + 1);
#endif
        // restrict searches first using search bins
        min_add = grid1_size;
        max_add = -1;
        for (n = 0; n < num_srch_bins; n++)
        {
            if (grid2_add >= bin_addr2[n*2] && 
                    grid2_add <= bin_addr2[n*2 + 1])
            {
                min_add = MIN(min_add, bin_addr1[n*2]);
                max_add = MAX(max_add, bin_addr1[n*2 + 1]);
            }
        }
        
        // further restrict searches using bounding boxes
        num_srch_cells = 0;     // extern global var
        //printf("from %d to %d\n", min_add, max_add);
        for (grid1_add = min_add; grid1_add <= max_add; grid1_add ++)
        {
            srch_mask[grid1_add] = true;
            srch_mask[grid1_add] &= 
                (grid1_bound_box[grid1_add * BOUNDBOX_SIZE] <= grid2_bound_box[grid2_add * BOUNDBOX_SIZE + 1]);
            srch_mask[grid1_add] &=
                (grid1_bound_box[grid1_add * BOUNDBOX_SIZE + 1] >= grid2_bound_box[grid2_add * BOUNDBOX_SIZE]);
            srch_mask[grid1_add] &=
                ( grid1_bound_box[grid1_add * BOUNDBOX_SIZE + 2] <= grid2_bound_box[grid2_add *BOUNDBOX_SIZE + 3]);
                //le( grid1_bound_box[grid1_add * BOUNDBOX_SIZE + 2], grid2_bound_box[grid2_add *BOUNDBOX_SIZE + 3]);
            srch_mask[grid1_add] &=
                ( grid1_bound_box[grid1_add *BOUNDBOX_SIZE + 3] >= grid2_bound_box[grid2_add * BOUNDBOX_SIZE + 2]);
                //ge( grid1_bound_box[grid1_add *BOUNDBOX_SIZE + 3], grid2_bound_box[grid2_add * BOUNDBOX_SIZE + 2]);
            if (_3D_BOUND_BOX_)         // using 3D bound box option
            {
                srch_mask[grid1_add] &=
                    (grid1_bound_box[grid1_add * BOUNDBOX_SIZE] + 4 <= grid2_bound_box[grid2_add * BOUNDBOX_SIZE + 5]);
                srch_mask[grid1_add] &=
                    (grid1_bound_box[grid1_add * BOUNDBOX_SIZE + 5] >= grid2_bound_box[grid2_add * BOUNDBOX_SIZE + 4]);
            }

            if (srch_mask[grid1_add])
                num_srch_cells ++;
        }
        // create search arrays
        n = 0;
        for (grid1_add = min_add; grid1_add <= max_add; grid1_add ++)
        {
            if (srch_mask[grid1_add])
            {
                srch_add[n] = grid1_add;
                for (corner = 0; corner < grid1_corners[grid1_add]; corner ++)
                {
                    srch_corner_lat[n * grid1_corners_max + corner] = 
                        grid1_corner_lat[grid1_add * grid1_corners_max + corner];
                    srch_corner_lon[n * grid1_corners_max + corner] = 
                        grid1_corner_lon[grid1_add * grid1_corners_max + corner];
                }
                n++;
            }
        }

        // integrate around this cell
        int num_corner = grid2_corners[grid2_add];
        for (corner = 0; corner < num_corner; corner ++)
        {
            next_corn = (corner + 1) %  num_corner;

            // define endpoints of the current segment
            int corndex = grid2_add * grid2_corners_max;
            beglat = grid2_corner_lat [corndex + corner];
            beglon = grid2_corner_lon [corndex + corner];
            endlat = grid2_corner_lat [corndex + next_corn];
            endlon = grid2_corner_lon [corndex + next_corn];
            lrevers = false;

            // to ensure exact path taken during both sweeps, always integrate segments in the same direction (SW TO NE)
            if ((endlat < beglat) || 
                    (endlat == beglat && endlon < beglon ))
            {
                beglat = grid2_corner_lat [corndex + next_corn];
                beglon = grid2_corner_lon [corndex + next_corn];
                endlat = grid2_corner_lat [corndex + corner];
                endlon = grid2_corner_lon [corndex + corner];
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
                        printf("Integrate stalled: num_subseg exceeded limit @ grid2 %d\n", grid2_add);
                    }

                    // find next intersection of this segment with a grid line on dst grid

#if _DEBUG_INPUT_OUTPUT_
                    // input compare for debug
                    if (false)
                    {
                    printf("%6d:(%3.3f %3.3f) (%3.3f %3.3f %3.3f %3.3f) (%3.3f %3.3f)",
                            grid1_add + 1, intrsct_lat, intrsct_lon,
                            beglat, beglon, endlat, endlon, begseg[0], begseg[1]);
                    if (lrevers)
                        printf(" true\n");
                    else
                        printf(" false\n");
                    }
#endif
                    intersection(grid1_add, intrsct_lat, intrsct_lon, lcoinc, beglat, beglon, endlat, endlon, begseg, lbegin, lrevers);
                    lbegin = false;
                    
                    // compute line integral for this subsegment
                    if (grid1_add != -1)
                    {
                        line_integral(conserv_weights, num_wts, beglon, intrsct_lon, beglat, intrsct_lat, 
                                grid1_center_lat[grid1_add],
                                grid1_center_lon[grid1_add],
                                grid2_center_lat[grid2_add],
                                grid2_center_lon[grid2_add]);
                    }
                    else
                    {
                        //printf("unhappy\n");
                        line_integral(conserv_weights, num_wts, beglon, intrsct_lon, beglat, intrsct_lat,
                                grid2_center_lat[grid2_add],
                                grid2_center_lon[grid2_add],
                                grid2_center_lat[grid2_add],
                                grid2_center_lon[grid2_add]);
                    }
                    //printf("%3.6f %3.6f vs %3.6f %3.6f\n", intrsct_lat, intrsct_lon, beglat, beglon);
                    // if integrating in reverse order, change sign of weights
                    if (lrevers)
                    {
                        for (int i = 0; i < num_wts * 2; i++)
                        {
                            conserv_weights[i] = - conserv_weights[i];
                        }
                    }

                    // store the appropriate addresses and weights. also add contributions to cell areas and centroids
                    if (grid1_add != -1 && !lcoinc)
                    {
                        if (grid1_mask[grid1_add])
                        {
#if _DEBUG_STORE_LINK_
                            printf("%6d--%6d: %3.8f %3.8f %3.8f %3.8f %3.8f %3.8f\n", grid1_add, grid2_add, 
                                    conserv_weights[0], conserv_weights[1], conserv_weights[2],
                                    conserv_weights[num_wts], conserv_weights[num_wts+1], conserv_weights[num_wts+2]);
#endif
                            store_link_cnsrv(grid1_add, grid2_add, conserv_weights, num_wts);
                            grid1_frac[grid1_add] += conserv_weights[0];
                            grid2_frac[grid2_add] += conserv_weights[num_wts];
                        }
                    }

#if _DEBUG_LINKS_
                    printf("%6d--%6d: %3.8f--%3.8f--%3.8f\n", grid2_add, grid1_add, 
                            conserv_weights[0], conserv_weights[1], conserv_weights[2]);
#endif
                    grid2_area[grid2_add] += conserv_weights[num_wts];
                    grid2_centroid_lat[grid2_add] += conserv_weights[num_wts + 1];
                    grid2_centroid_lon[grid2_add] += conserv_weights[num_wts + 2];

                    // reset beglat and beglon for next subsegment
                    beglat = intrsct_lat;
                    beglon = intrsct_lon;
                }
            }
        }
    }

    delete [] srch_mask;
    
#if _DEBUG_INITIAL_WEIGHTS_ 
    // check weights
    for (int n = 0; n < num_links_map; n++)
    {
        grid1_add = grid1_add_map[n];
        grid2_add = grid2_add_map[n];
        printf("%6d %6d  %3.8f %3.8f %3.8f\n", grid1_add+1, grid2_add+1, wts_map[num_wts * n], wts_map[num_wts * n + 1], wts_map[num_wts * n + 2]);
    }
#endif
    /* correct for situations where N/S pole not explicitly included in grid (i.e. as a grid corner point). if pole is missing from only one grid, need to correct only the area and centroid of that grid. if missing from both, do complete weight calculation */
    int grid1_pole = -1, grid2_pole = -1;
    conserv_weights[0] = PI2;
    conserv_weights[1] = PI * PI;
    conserv_weights[2] = ZERO;
    conserv_weights[3] = PI2;
    conserv_weights[4] = PI * PI;
    conserv_weights[5] = ZERO;

    // for north pole
    for (int i = 0; i < grid1_size; i++)
    {
        if (grid1_area[i] < - THREE * PIH &&
            grid1_center_lat[i] > ZERO)
        {
            grid1_pole = i;
            break;
        }
    }
    for (int i = 0; i < grid2_size; i++)
    {
        if (grid2_area[i] < -THREE * PIH &&
            grid2_center_lat[i] > ZERO)
        {
            grid2_pole = i;
            break;
        }
    }
    if (grid1_pole != -1)
    {
        grid1_area[grid1_pole] += conserv_weights[0];
        grid1_centroid_lat[grid1_pole] += conserv_weights[1];
        grid1_centroid_lon[grid1_pole] += conserv_weights[2];
    }
    if (grid2_pole != -1)
    {
        grid2_area[grid2_pole] += conserv_weights[num_wts];
        grid2_centroid_lat[grid2_pole] += conserv_weights[num_wts + 1];
        grid2_centroid_lon[grid2_pole] += conserv_weights[num_wts + 2];
    }
    if (grid1_pole != -1 && grid2_pole != -1)
    {
        store_link_cnsrv(grid1_pole, grid2_pole, conserv_weights, num_wts);
        grid1_frac[grid1_pole] += conserv_weights[0];
        grid2_frac[grid2_pole] += conserv_weights[num_wts];
    }

    // for south pole
    grid1_pole = -1;
    grid2_pole = -1;
    conserv_weights[0] = PI2;
    conserv_weights[1] = - PI * PI;
    conserv_weights[2] = ZERO;
    conserv_weights[3] = PI2;
    conserv_weights[4] = - PI * PI;
    conserv_weights[5] = ZERO;

    for (int i = 0; i < grid1_size; i++)
    {
        if (grid1_area[i] < - THREE * PIH &&
            grid1_center_lat[i] < ZERO)
        {
            grid1_pole = i;
            break;
        }
    }
    for (int i = 0; i < grid2_size; i++)
    {
        if (grid2_area[i] < -THREE * PIH &&
            grid2_center_lat[i] < ZERO)
        {
            grid2_pole = i;
            break;
        }
    }
    if (grid1_pole != -1)
    {
        grid1_area[grid1_pole] += conserv_weights[0];
        grid1_centroid_lat[grid1_pole] += conserv_weights[1];
        grid1_centroid_lon[grid1_pole] += conserv_weights[2];
    }
    if (grid2_pole != -1)
    {
        grid2_area[grid2_pole] += conserv_weights[num_wts];
        grid2_centroid_lat[grid2_pole] += conserv_weights[num_wts + 1];
        grid2_centroid_lon[grid2_pole] += conserv_weights[num_wts + 2];
    }
    if (grid1_pole != -1 && grid2_pole != -1)
    {
        store_link_cnsrv(grid1_pole, grid2_pole, conserv_weights, num_wts);
        grid1_frac[grid1_pole] += conserv_weights[0];
        grid2_frac[grid2_pole] += conserv_weights[num_wts];
    }
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
            conserv_weights[nwgt] = wts_map[n*num_wts + nwgt];
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
                if (nonzero(grid2_frac[grid2_add]))
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

        wts_map[n*num_wts] = conserv_weights[0] * norm_factor;
        wts_map[n*num_wts + 1] = (conserv_weights[1] - conserv_weights[0]*grid1_centroid_lat[grid1_add]) * norm_factor;
        wts_map[n*num_wts + 2] = (conserv_weights[2] - conserv_weights[0]*grid1_centroid_lon[grid1_add]) * norm_factor;
    }    

    printf("Total number of links = %d\n", num_links_map);
    for (int n = 0; n < grid1_size; n++)
        if (nonzero(grid1_area[n]))
            grid1_frac[n] /= grid1_area[n];

    for (int n = 0; n < grid2_size; n++)
        if (nonzero(grid2_area[n]))
            grid2_frac[n] /= grid2_area[n];

#if _DEBUG_FINAL_WEIGHTS_ 
// check weights
for (int n = 0; n < num_links_map; n++)
{
    grid1_add = grid1_add_map[n];
    grid2_add = grid2_add_map[n];
    printf("%6d %6d  %3.8f %3.8f %3.8f\n", grid1_add+1, grid2_add+1, wts_map[num_wts * n], wts_map[num_wts * n + 1], wts_map[num_wts * n + 2]);
}
#endif


#if _CHECK_AREA_CENTROID_
    /* perform some error checking on final weights */
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
    for (int n = 0; n < num_links_map; n++)
    {
        grid1_add = grid1_add_map[n];
        grid2_add = grid2_add_map[n];

        if (wts_map[num_wts * n] < -0.01)
            printf("Map Weight < 0 @grid1_add=%d\tgrid2_add=%d\t valued %3.6f\n", grid1_add, grid2_add, wts_map[num_wts * n]);
        if (norm_opt != NORM_OPT_NONE && wts_map[num_wts * n] > 1.01)
            printf("Map weight > 1 @grid1_add=%d\tgrid2_add=%d\t valued %3.6f\n", grid1_add, grid2_add, wts_map[num_wts * n]);

        grid2_centroid_lat [grid2_add] += wts_map[num_wts * n];
    }

    for (int n = 0; n < grid2_size; n++)
    {
        switch (norm_opt)
        {
            case NORM_OPT_DESTAREA:
                norm_factor = grid2_frac[n];
                break;
            case NORM_OPT_FRACAREA:
                norm_factor = ONE;
                break;
            case NORM_OPT_NONE:
                if (luse_grid2_area)
                    norm_factor = grid2_area_in [n];
                else
                    norm_factor = grid2_area [n];
                break;
            default:
                break;
        }
        if (abs(grid2_centroid_lat[n] - norm_factor ) > 0.01)
            printf("Error: sum of wts for map @%d\t %3.6f\t vs %3.6f\n", grid2_add, grid2_centroid_lat[grid2_add], norm_factor);
    }
#endif
    delete [] grid1_centroid_lat;
    delete [] grid1_centroid_lon;
    delete [] grid2_centroid_lat;
    delete [] grid2_centroid_lon;
    delete [] begseg;
    delete [] conserv_weights;
}   
