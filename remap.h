#ifndef _REMAP_
#define _REMAP_ 1

#include "kinds.h"
#include "constants.h"

// 1-order remap: src grid to dst grid, using weights, cell values and grad values
void remap(double *dst_array, double *src_array, int *dst_add, int *src_add, double *map_wts, int size);

// 2-order conserv remap:
void remap(double* dst_array, double *src_array, int *dst_add, int *src_add, double *map_wts, int size, 
        double *src_grad1, double *src_grad2);

// 2-order remap: src grid to dst grid, using weights, cell values and grad values
void remap(double *dst_array, double *src_array, int *dst_add, int *src_add, double *map_wts, int size, 
        double *src_grad1, double *src_grad2, double *src_grad3);
#endif
