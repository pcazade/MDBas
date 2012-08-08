#ifndef VDWH_INCLUDED
#define VDWH_INCLUDED

void vdw_full(ATOM *atom,FORCEFIELD *ff,ENERGY *ener,SIMULPARAMS *simulCond,PBC *box);
void vdw_switch(ATOM *atom,FORCEFIELD *ff,ENERGY *ener,SIMULPARAMS *simulCond,PBC *box);
void vdw14_full(ATOM *atom,FORCEFIELD *ff,ENERGY *ener,SIMULPARAMS *simulCond,PBC *box);
void vdw14_switch(ATOM *atom,FORCEFIELD *ff,ENERGY *ener,SIMULPARAMS *simulCon,PBC *boxd);

#endif