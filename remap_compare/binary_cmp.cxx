#include <stdio.h>
#include <iostream>
#include <math.h>
#include <string.h>

using namespace std;

int main(int argc, char **argv)
{
    FILE *fp;
    int buf_size;
    if (argc < 4)
    {
        printf("Usage: ./binary_cmp file1 file2 choice\n");
        return 1;
    }
    double * scrip, * me;
/*
    if (strcmp(argv[3],"1") == 0)
        buf_size = 63112 * 3;
    else if (strcmp(argv[3], "2") == 0)
        buf_size = (8192+24576) * 10;
    printf("buf_size:%d\n", buf_size);
    scrip = new double[buf_size];
    me = new double[buf_size];
*/

    printf("opening %s\n", argv[1]);
    fp = fopen(argv[1], "rb");
    if (fp == (FILE *)0)
    {
        printf("file %s does not exist!\n", argv[1]);
        return 1;
    }
    fseek(fp, 0L, SEEK_END);
    buf_size = ftell(fp) / sizeof(double);
    fseek(fp, 0L, SEEK_SET);    // set back
    scrip = new double[buf_size];
    fread(scrip, sizeof(double), buf_size, fp);

    printf("opening %s\n", argv[2]);
    fp = fopen(argv[2], "rb");
    if (fp == (FILE *)0)
    {
        printf("file %s does not exist!\n", argv[2]);
        return 1;
    }
    fseek(fp, 0L, SEEK_END);
    if (ftell(fp) / sizeof(double) != buf_size)
    {
        printf("file.scrip and file.me do not match! please CHECK!\n");
        return 1;
    }
    fseek(fp, 0L, SEEK_SET);    // set back
    me = new double[buf_size];
    fread(me, sizeof(double), buf_size, fp);
    fclose(fp);
    
    if (strcmp(argv[3], "1") == 0)
    {
        // check weights info
        double max_error1, max_error2, max_error3;
        int indx_error1, indx_error2, indx_error3;
        double error;
        int num_links = 63112;
        max_error1 = max_error2 = max_error3 = -1.0;
        indx_error1 = indx_error2 = indx_error3 = -1;
        for (int i = 0; i < num_links; i++)
        {
            error = - scrip[i*3] + me[i*3];
            if (fabs(error) > max_error1)
            {
                indx_error1 = i;
                max_error1 = fabs(error);
            }
            error = me[i*3 + 1] - scrip[i*3 + 1];
            if (fabs(error) > max_error2)
            {
                indx_error2 = i;
                max_error2 = fabs(error);
            }
            error = me[i*3 + 2] - scrip[i*3 + 2];
            if (fabs(error) > max_error3)
            {
                indx_error3 = i;
                max_error3 = fabs(error);
            }
            //printf("%2.30f\t%2.30f\n", scrip[i], me[i]);
        }
        printf("error max 1 @%d\n\tscrip:\t%2.30f\n\tme:\t%2.30f\n", indx_error1, scrip[indx_error1 * 3], me[indx_error1 * 3]);
        printf("error max 2 @%d\n\tscrip:\t%2.30f\n\tme:\t%2.30f\n", indx_error2, scrip[indx_error1 * 3 + 1], me[indx_error1 * 3 + 1]);
        printf("error max 3 @%d\n\tscrip:\t%2.30f\n\tme:\t%2.30f\n", indx_error3, scrip[indx_error1 * 3 + 2], me[indx_error1 * 3 + 2]);
    }

    else if (strcmp(argv[3], "2") == 0)
    {
        // check grid info including center/corner lat/lon
        int grid1_size = 8192, grid2_size = 24576;
        double *head_scrip, *head_me;
        head_scrip = scrip;
        head_me = me;
        double cntr_lat_error, cntr_lon_error, crnr_lat_error, crnr_lon_error;
        int indx_cntr_lat, indx_cntr_lon, indx_crnr_lat, indx_crnr_lon;
        double error;
        cntr_lat_error = cntr_lon_error = crnr_lat_error = crnr_lon_error = -1.0;
        indx_cntr_lat = indx_cntr_lon = indx_crnr_lat = indx_crnr_lon = -1;

        //printf("-PIH\n\t%2.30f\n\t%2.30f\n", head_scrip[0], head_me[0]);
        //printf("deg2rad\n\t%2.30f\n\t%2.30f\n", head_scrip[1], head_me[1]);
        //head_scrip += 2;
        //head_me += 2;

        // grid1
        for (int i = 0; i < grid1_size; i++)
        {
            error = head_me[i] - head_scrip[i];
            if (fabs(error) > indx_cntr_lat)
            {
                cntr_lat_error = fabs(error);
                indx_cntr_lat = i;
            }
        }
        printf("grid1 cntr_lat max_error @ %d\n\t%2.30f\n\t%2.30f\n", indx_cntr_lat, head_scrip[indx_cntr_lat], head_me[indx_cntr_lat]);
        head_scrip += grid1_size;
        head_me += grid1_size;

        for (int i = 0; i < grid1_size; i++)
        {
            error = head_me[i] - head_scrip[i];
            if (fabs(error) > indx_cntr_lon)
            {
                cntr_lon_error = fabs(error);
                indx_cntr_lon = i;
            }
        }
        printf("grid1 cntr_lon max_error @ %d\n\t%2.30f\n\t%2.30f\n", indx_cntr_lon, head_scrip[indx_cntr_lon], head_me[indx_cntr_lon]);
        head_scrip += grid1_size;
        head_me += grid1_size;

        int grid1_corners = grid1_size * 4;
        for (int i = 0; i < grid1_corners; i++)
        {
            error = head_me[i] - head_scrip[i];
            if (fabs(error) > indx_crnr_lat)
            {
                crnr_lat_error = fabs(error);
                indx_crnr_lat = i;
            }
        }
        printf("grid1 crnr_lat max_error @ %d\n\t%2.30f\n\t%2.30f\n", indx_crnr_lat, head_scrip[indx_crnr_lat], head_me[indx_crnr_lat]);
        head_scrip += grid1_corners;
        head_me += grid1_corners;

        for (int i = 0; i < grid1_corners; i++)
        {
            error = head_me[i] - head_scrip[i];
            if (fabs(error) > indx_crnr_lon)
            {
                crnr_lon_error = fabs(error);
                indx_crnr_lon = i;
            }
        }
        printf("grid1 crnr_lon max_error @ %d\n\t%2.30f\n\t%2.30f\n", indx_crnr_lon, head_scrip[indx_crnr_lon], head_me[indx_crnr_lon]);
        head_scrip += grid1_corners;
        head_me += grid1_corners;
        
        // grid2
        cntr_lat_error = cntr_lon_error = crnr_lat_error = crnr_lon_error = -1.0;
        for (int i = 0; i < grid2_size; i++)
        {
            error = head_me[i] - head_scrip[i];
            if (fabs(error) > indx_cntr_lat)
            {
                cntr_lat_error = fabs(error);
                indx_cntr_lat = i;
            }
        }
        printf("grid2 cntr_lat max_error @ %d\n\t%2.30f\n\t%2.30f\n", indx_cntr_lat, head_scrip[indx_cntr_lat], head_me[indx_cntr_lat]);
        head_scrip += grid2_size;
        head_me += grid2_size;

        for (int i = 0; i < grid2_size; i++)
        {
            error = head_me[i] - head_scrip[i];
            if (fabs(error) > indx_cntr_lon)
            {
                cntr_lon_error = fabs(error);
                indx_cntr_lon = i;
            }
        }
        printf("grid2 cntr_lon max_error @ %d\n\t%2.30f\n\t%2.30f\n", indx_cntr_lon, head_scrip[indx_cntr_lon], head_me[indx_cntr_lon]);
        head_scrip += grid2_size;
        head_me += grid2_size;

        int grid2_corners = grid2_size * 4;
        for (int i = 0; i < grid2_corners; i++)
        {
            error = head_me[i] - head_scrip[i];
            if (fabs(error) > indx_crnr_lat)
            {
                crnr_lat_error = fabs(error);
                indx_crnr_lat = i;
            }
        }
        printf("grid2 crnr_lat max_error @ %d\n\t%2.30f\n\t%2.30f\n", indx_crnr_lat, head_scrip[indx_crnr_lat], head_me[indx_crnr_lat]);
        head_scrip += grid2_corners;
        head_me += grid2_corners;

        for (int i = 0; i < grid2_corners; i++)
        {
            error = head_me[i] - head_scrip[i];
            if (fabs(error) > indx_crnr_lon)
            {
                crnr_lon_error = fabs(error);
                indx_crnr_lon = i;
            }
        }
        printf("grid2 crnr_lon max_error @%d\n\t%2.30f\n\t%2.30f\n", indx_crnr_lon, head_scrip[indx_crnr_lon], head_me[indx_crnr_lon]);
    }
    delete [] me;
    delete [] scrip;
    return 0;
}
