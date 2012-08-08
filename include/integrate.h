#ifndef INTEGRATEH_INCLUDED
#define INTEGRATEH_INCLUDED

void lf_integrate(ATOM *atom, ENERGY *ener, SIMULPARAMS *simulCond,CONSTRAINT *constList,PBC *box);
void lf_nve(ATOM *atom, ENERGY *ener, SIMULPARAMS *simulCond,CONSTRAINT *constList,PBC *box);
void lf_nvt_b(ATOM *atom, ENERGY *ener, SIMULPARAMS *simulCond,CONSTRAINT *constList,PBC *box);

void vv_integrate(ATOM *atom, ENERGY *ener, SIMULPARAMS *simulCond,CONSTRAINT *constList,PBC *box,int stage);
void vv_nve(ATOM *atom, ENERGY *ener, SIMULPARAMS *simulCond,CONSTRAINT *constList,PBC *box,int stage);
void vv_nvt_b(ATOM *atom, ENERGY *ener, SIMULPARAMS *simulCond,CONSTRAINT *constList,PBC *box,int stage);

#endif