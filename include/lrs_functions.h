// code borrowed from LRS's chdemo.c, copyright David Avis 2001

// before including lrslib.h, GMP must be defined
#include <lrslib.h>

#define MAXCOL 1000     /* maximum number of colums */

void makecyclic(lrs_dic*,lrs_dat*);
int ch();
// vim: ts=2
