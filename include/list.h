#ifndef LISTH_INCLUDED
#define LISTH_INCLUDED

void exclude_list(SIMULPARAMS *simulCond,ATOM *atom,FORCEFIELD *ff);
void verlet_list(SIMULPARAMS *simulCond,ATOM *atom,FORCEFIELD *ff);
void verlet_list_update(SIMULPARAMS *simulCond,ATOM *atom,FORCEFIELD *ff);

#endif