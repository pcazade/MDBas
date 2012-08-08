#ifndef NUMDERIVH_INCLUDED
#define NUMDERIVH_INCLUDED

void numforce(ATOM *atom,DELTA *nForce,FORCEFIELD *ff,ENERGY *ener,SIMULPARAMS *simulCond,PBC *box,int npoints,double h);

#endif