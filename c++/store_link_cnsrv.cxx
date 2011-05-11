#include "store_link_cnsrv.h"
#include "remap_vars.h"
#include "utils.h"

bool first_call_store = true;
int *src_link_add, *dst_link_add;   // min|max link add to restrict search

/** this routine stores the address and weight for this link in
  * the appropriate address and weight arrays and resizes those
  * arrays if necessary
  **/
void store_link_cnsrv(int add1, int add2, double *weights, int num_weights)
{
    // local variables
    int nlink, min_link, max_link;      // index for link
    int map_id;

    bool allzero = true;
    for (int i = 0; i < num_weights; i++)
    {
        allzero &= zero(weights[i]);
    }

    // if all weights are zero, do not bother storing the link
    if (allzero)
        return;

    // restrict the range of links to search for existing links
    if (first_call_store)
    {
        src_link_add = new int [grid1_size * 2];
        dst_link_add = new int [grid2_size * 2];
        memset(src_link_add, -1, sizeof(int) * grid1_size * 2);
        memset(dst_link_add, -1, sizeof(int) * grid2_size * 2);
        first_call_store = false;
        min_link = 1;       // initial values
        max_link = 0;
    }
    else
    {
        min_link = MIN(src_link_add[add1 * 2], dst_link_add[add2 * 2]);
        max_link = MAX(src_link_add[add1*2+1], dst_link_add[add2*2+1]);
        if (min_link == -1)
        {
            min_link = 1;
            max_link = 0;
        }
    }

    // if the link already exists, add the weight to the current weight arrays
    for (nlink = min_link; nlink <= max_link; nlink ++)
    //for (nlink = 0; nlink < num_links_map; nlink ++)
    {
        if ((add1 == grid1_add_map[nlink]) &&
            (add2 == grid2_add_map[nlink]))
        {
            for (int i = 0; i < num_wts; i++)
            {
                wts_map[nlink * num_wts + i] += weights[i];
            }
            return;
        }
    }

    // if the link does not yet exist, increment number of links and check to see if remap arrays need to be increased to accomodate the new link. then store the link
    map_id = num_links_map;
    num_links_map ++;
    if (num_links_map > max_links_map)
        resize_remap_vars(resize_increment);
    grid1_add_map[map_id] = add1;
    grid2_add_map[map_id] = add2;
    
    for (int i = 0; i < num_wts; i++)
    {
        wts_map[map_id*num_wts + i] = weights[i];
    }

    if (src_link_add[add1 * 2] == -1)
        src_link_add[add1 * 2] = map_id;
    if (dst_link_add[add2 * 2] == -1)
        dst_link_add[add2 * 2] = map_id;
    src_link_add[add1 * 2 + 1] = map_id;
    dst_link_add[add2 * 2 + 1] = map_id;
}
