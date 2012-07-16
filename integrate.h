#ifndef INTEGRATEH_INCLUDED
#define INTEGRATEH_INCLUDED

void lf_nve(ATOM *atom, ENERGYFORCE *enerFor, SIMULPARAMS *simulCond);
void vv_nve(ATOM *atom, ENERGYFORCE *enerFor, SIMULPARAMS *simulCond,int stage);

#endif