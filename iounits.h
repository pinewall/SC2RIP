#ifndef _IOUNITS_
#define _IOUNITS_ 1

#include "kinds.h"

extern bool *unit_free;       // flags to determine whether unit is free for use
const int _stdin = 5;    // reserves unit for standard input
const int _stdout = 6;   // reserves unit for standard output
const int _stderr = 6;   // reserves unit for standard error

int get_unit();         // return the next available I/O unit number
void release_unit(int iunit);   //releaase the specified unit and close the file

#endif
