/**
 * \file rand.h
 * \brief Prototypes for file rand.c
 * \author Pierre-Andre Cazade and Florent Hedin
 * \version alpha-branch
 * \date 2012
 */

#ifndef RANDH_INCLUDED
#define RANDH_INCLUDED

void init_rand(unsigned int seed);
double get_rand();
void get_BoxMuller(double *u, double *v);

#endif