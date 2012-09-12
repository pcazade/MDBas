/**
 * \file energy.h
 * \brief Contains highest level functions for evaluating energy and forces of a system.
 * \author Pierre-Andre Cazade and Florent Hedin
 * \version alpha-branch
 * \date 2012
 */

#ifndef ENERGYH_INCLUDED
#define ENERGYH_INCLUDED

void init_energy_ptrs(SIMULPARAMS *simulCond);
void energy(ATOM atom[],FORCEFIELD *ff,ENERGY *ener,SIMULPARAMS *simulCond,PBC *box);
void nonbond_energy(ATOM atom[],FORCEFIELD *ff,ENERGY *ener,SIMULPARAMS *simulCond,PBC *box);
void nonbond14_energy(ATOM atom[],FORCEFIELD *ff,ENERGY *ener,SIMULPARAMS *simulCond,PBC *box);

/* pointers to function for electrostatic energy functions */
double (*coulomb)(ATOM atom[],FORCEFIELD *ff,SIMULPARAMS *simulCond,PBC *box,int i,int j,double r,double *delec);
double (*coulomb14)(ATOM atom[],FORCEFIELD *ff,SIMULPARAMS *simulCond,PBC *box,int i,int j,double r,double *delec);

/* pointers to function for vdw energy functions */
double (*vdw)(ATOM atom[],FORCEFIELD *ff,SIMULPARAMS *simulCond,PBC *box,int i,int j,double r,double *dvdw);
double (*vdw14)(ATOM atom[],FORCEFIELD *ff,SIMULPARAMS *simulCond,PBC *box,int i,int j,double r,double *dvdw);

#endif