#ifndef _LABNH_H
#define _LABNH_H

int nhinit();
void nhfinish();

void nh_verbosity(int verbosity);
void nh_sampling(int sampling);

double nh_equ(double ra, double dec);
double nh_gal(double lii, double bii);

#endif

