#include "iounits.h"

bool first_call = true;
bool * unit_free = new bool[99];
int get_unit()
{
    int n;
    if (first_call)
    {
        for (n = 0; n < 99; n++)
            unit_free[n] = true;
        unit_free[_stdin] = false;
        unit_free[_stdout] = false;
        unit_free[_stderr] = false;
        first_call = false;
    }

    n = 0;
    while(! unit_free[n])
    {
        n ++;
    }
    unit_free[n] = false;
    return n;
}

void release_unit(int iunit)
{
    unit_free[iunit] = true;
    //close file in source code
}
