#ifndef ELECH_INCLUDED
#define ELECH_INCLUDED

double coulomb_none(ATOM *atom,FORCEFIELD *ff,SIMULPARAMS *simulCond,PBC *box,int i,int j,double r,double *delec);
void coulomb_full(ATOM *atom,FORCEFIELD *ff,ENERGY *ener,SIMULPARAMS *simulCond,PBC *box);
double coulomb_shift1(ATOM *atom,FORCEFIELD *ff,SIMULPARAMS *simulCond,PBC *box,int i,int j,double r,double *delec);
double coulomb_shift2(ATOM *atom,FORCEFIELD *ff,SIMULPARAMS *simulCond,PBC *box,int i,int j,double r,double *delec);
double coulomb_switch(ATOM *atom,FORCEFIELD *ff,SIMULPARAMS *simulCond,PBC *box,int i,int j,double r,double *delec);

double coulomb14_none(ATOM *atom,FORCEFIELD *ff,SIMULPARAMS *simulCond,PBC *box,int i,int j,double r,double *delec);
double coulomb14_full(ATOM *atom,FORCEFIELD *ff,SIMULPARAMS *simulCond,PBC *box,int i,int j,double r,double *delec);
double coulomb14_shift1(ATOM *atom,FORCEFIELD *ff,SIMULPARAMS *simulCond,PBC *box,int i,int j,double r,double *delec);
double coulomb14_shift2(ATOM *atom,FORCEFIELD *ff,SIMULPARAMS *simulCond,PBC *box,int i,int j,double r,double *delec);
double coulomb14_switch(ATOM *atom,FORCEFIELD *ff,SIMULPARAMS *simulCond,PBC *box,int i,int j,double r,double *delec);

#endif