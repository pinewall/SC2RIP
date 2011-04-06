#ifndef _STORE_LINK_CNSRV_
#define _STORE_LINK_CNSRV_ 1

/** this routine stores the address and weight for this link in
  * the appropriate address and weight arrays and resizes those
  * arrays if necessary
  **/
void store_link_cnsrv(int add1, int add2, double *weights, int num_weights);

#endif
