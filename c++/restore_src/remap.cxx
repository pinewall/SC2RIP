#include "remap.h"
#include "remap_vars.h"

// 1-order remap: src grid to dst grid, using weights, cell values and grad values
void remap(double *dst_array, double *src_array, int *dst_add, int *src_add, double *map_wts, int size)
{
    remap(dst_array, src_array, dst_add, src_add, map_wts, size, 0, 0, 0);
}

// 2-order conserv remap
void remap(double* dst_array, double *src_array, int *dst_add, int *src_add, double *map_wts, int size, 
        double *src_grad1, double *src_grad2)
{
    remap(dst_array, src_array, dst_add, src_add, map_wts, size, src_grad1, src_grad2, 0);
}

// 2-order remap src grid to dst grid, using weights, cell values and grad values
void remap(double *dst_array, double *src_array, int *dst_add, int *src_add, double *map_wts, int size, 
        double *src_grad1, double *src_grad2, double *src_grad3)
{
    int iorder = 1;         // conservative remapping order; 1-order default
    if (src_grad1 != 0 && src_grad2 != 0 && src_grad3 != 0)
        iorder = 2;         // provided grad values, uses 2-order conservative remapping

    int src, dst;
    // first order remapping
    if (iorder == 1)
    {
        for (int i = 0; i < size; i++)
        {
            src = src_add[i];
            dst = dst_add[i];
            dst_array[dst] += src_array[src] * map_wts[i*num_wts]; 
        }
    }

    // second order remapping
    else if (iorder == 2)       // 2-order conservative
    {
        if (num_wts == 3)
        {
            for (int i = 0; i < size; i++)
            {
                src = src_add[i];
                dst = dst_add[i];
                dst_array[dst] += (src_array[src] * map_wts[i*num_wts] +
                                    src_grad1[src] * map_wts[i*num_wts+1] + 
                                    src_grad2[src] * map_wts[i*num_wts+2]);

            }
        }
        else if (num_wts == 4)  // 2-order bicubic
        {
            for (int i = 0; i < size; i++)
            {
                src = src_add[i];
                dst = dst_add[i];
                dst_array[dst] += (src_array[src] * map_wts[i*num_wts] +
                                    src_grad1[src] * map_wts[i*num_wts+1] + 
                                    src_grad2[src] * map_wts[i*num_wts+2] +
                                    src_grad3[src] * map_wts[i*num_wts+3]);

            }
        }
    }
}
