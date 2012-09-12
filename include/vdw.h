/**
 * \file vdw.h
 * \brief Prototypes for file vdw.c
 * \author Pierre-Andre Cazade and Florent Hedin
 * \version alpha-branch
 * \date 2012
 */

#ifndef VDWH_INCLUDED
#define VDWH_INCLUDED

double vdw_none(ATOM atom[],FORCEFIELD *ff,SIMULPARAMS *simulCond,PBC *box,int i,int j,double r,double *dvdw);
void vdw_full(ATOM atom[],FORCEFIELD *ff,ENERGY *ener,SIMULPARAMS *simulCond,PBC *box);
double vdw_switch(ATOM atom[],FORCEFIELD *ff,SIMULPARAMS *simulCond,PBC *box,int i,int j,double r,double *dvdw);

double vdw14_none(ATOM atom[],FORCEFIELD *ff,SIMULPARAMS *simulCond,PBC *box,int i,int j,double r,double *dvdw);
double vdw14_full(ATOM atom[],FORCEFIELD *ff,SIMULPARAMS *simulCond,PBC *box,int i,int j,double r,double *dvdw);
double vdw14_switch(ATOM atom[],FORCEFIELD *ff,SIMULPARAMS *simulCon,PBC *boxd,int i,int j,double r,double *dvdw);

#endif