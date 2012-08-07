#ifndef VDWH_INCLUDED
#define VDWH_INCLUDED

void vdw_full(ATOM *atom,FORCEFIELD *ff,ENERGY *ener,SIMULPARAMS *simulCond);
void vdw_switch(ATOM *atom,FORCEFIELD *ff,ENERGY *ener,SIMULPARAMS *simulCond);
void vdw14_full(ATOM *atom,FORCEFIELD *ff,ENERGY *ener,SIMULPARAMS *simulCond);
void vdw14_switch(ATOM *atom,FORCEFIELD *ff,ENERGY *ener,SIMULPARAMS *simulCond);

#endif