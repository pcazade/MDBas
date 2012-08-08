#ifndef INTERNALH_INCLUDED
#define INTERNALH_INCLUDED

void bond_energy(ATOM *atom,FORCEFIELD *ff,ENERGY *ener,SIMULPARAMS *simulCond,PBC *box);
void ub_energy(ATOM *atom,FORCEFIELD *ff,ENERGY *ener,SIMULPARAMS *simulCond,PBC *box);
void angle_energy(ATOM *atom,FORCEFIELD *ff,ENERGY *ener,SIMULPARAMS *simulCond,PBC *box);
void dihedral_energy(ATOM *atom,FORCEFIELD *ff,ENERGY *ener,SIMULPARAMS *simulCond,PBC *box);
void improper_energy(ATOM *atom,FORCEFIELD *ff,ENERGY *ener,SIMULPARAMS *simulCond,PBC *box);

#endif