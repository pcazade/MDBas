#ifndef RANDH_INCLUDED
#define RANDH_INCLUDED

void init_rand(unsigned int seed);
double get_rand();
void get_BoxMuller(double *u, double *v);

#endif