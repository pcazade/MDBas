/**
 * \file numderiv.h
 * \brief Prototypes for file numderiv.c
 * \author Pierre-Andre Cazade and Florent Hedin
 * \version alpha-branch
 * \date 2012
 */

#ifndef NUMDERIVH_INCLUDED
#define NUMDERIVH_INCLUDED

void numforce(ATOM atom[],DELTA *nForce,FORCEFIELD *ff,ENERGY *ener,SIMULPARAMS *simulCond,PBC *box,int npoints,double h);

#endif