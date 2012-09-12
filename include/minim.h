/**
 * \file minim.h
 * \brief Prototypes for file minim.c
 * \author Pierre-Andre Cazade and Florent Hedin
 * \version alpha-branch
 * \date 2012
 */

#ifndef MINIMH_INCLUDED
#define MINIMH_INCLUDED

void minimise(ATOM atom[],FORCEFIELD *ff,ENERGY *ener,SIMULPARAMS *simulCond);

void steepestDescent(ATOM atom[],FORCEFIELD *ff,ENERGY *ener,SIMULPARAMS *simulCond,PBC *box);

void conjugateGradients(ATOM atom[],FORCEFIELD *ff,ENERGY *ener,SIMULPARAMS *simulCond);

#endif