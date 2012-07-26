#ifndef LISTH_INCLUDED
#define LISTH_INCLUDED

void makelist(SIMULPARAMS *simulCond,ATOM *atom,FORCEFIELD *ff,CONSTRAINT *constList);
void exclude_list(SIMULPARAMS *simulCond,ATOM *atom,FORCEFIELD *ff,CONSTRAINT *constList);
void verlet_list(SIMULPARAMS *simulCond,ATOM *atom,FORCEFIELD *ff);
void verlet_list_update(SIMULPARAMS *simulCond,ATOM *atom,FORCEFIELD *ff);

#endif