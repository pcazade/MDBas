#ifndef INTERNALH_INCLUDED
#define INTERNALH_INCLUDED

void bond_energy(ATOM *atom,FORCEFIELD *ff,ENERGY *ener,SIMULPARAMS *simulCond);
void ub_energy(ATOM *atom,FORCEFIELD *ff,ENERGY *ener,SIMULPARAMS *simulCond);
void angle_energy(ATOM *atom,FORCEFIELD *ff,ENERGY *ener,SIMULPARAMS *simulCond);
void dihedral_energy(ATOM *atom,FORCEFIELD *ff,ENERGY *ener,SIMULPARAMS *simulCond);
void improper_energy(ATOM *atom,FORCEFIELD *ff,ENERGY *ener,SIMULPARAMS *simulCond);

#endif