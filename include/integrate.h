#ifndef INTEGRATEH_INCLUDED
#define INTEGRATEH_INCLUDED

void lf_integrate(ATOM *atom, ENERGYFORCE *enerFor, SIMULPARAMS *simulCond,CONSTRAINT *constList);
void lf_nve(ATOM *atom, ENERGYFORCE *enerFor, SIMULPARAMS *simulCond,CONSTRAINT *constList);
void lf_nvt_b(ATOM *atom, ENERGYFORCE *enerFor, SIMULPARAMS *simulCond,CONSTRAINT *constList);

void vv_integrate(ATOM *atom, ENERGYFORCE *enerFor, SIMULPARAMS *simulCond,CONSTRAINT *constList,int stage);
void vv_nve(ATOM *atom, ENERGYFORCE *enerFor, SIMULPARAMS *simulCond,CONSTRAINT *constList,int stage);

#endif