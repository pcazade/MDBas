#include <stdio.h>
#include <stdlib.h>
#include "global.h"
#include "elec.h"
#include "vdw.h"
#include "internal.h"
#include "io.h"

void energy(ATOM *atom,FORCEFIELD *ff,ENERGY *ener,SIMULPARAMS *simulCond,PBC *box)
{
  
  int i;
  
  ener->elec=0.;
  ener->vdw=0.;
  ener->bond=0.;
  ener->ub=0.;
  ener->ang=0.;
  ener->dihe=0.;
  ener->impr=0.;
  
  for(i=0;i<simulCond->natom;i++)
  {
    atom[i].fx=0.;
    atom[i].fy=0.;
    atom[i].fz=0.;
  }
  
  /* Performing electrostatic interactions */
  
  if(simulCond->elecType==NOELEC)
  {
    printf("No electrostatics requested\n");
  }
  else if (simulCond->elecType==FULL)
  {
    coulomb_full(atom,ff,ener,simulCond,box);
    if(simulCond->nb14)
      coulomb14_full(atom,ff,ener,simulCond,box);
  }
  else if(simulCond->elecType==SHIFT1)
  {
    coulomb_shift1(atom,ff,ener,simulCond,box);
    
    if(simulCond->nb14)
      coulomb14_shift1(atom,ff,ener,simulCond,box);
  }
  else if(simulCond->elecType==SHIFT2)
  {
    coulomb_shift2(atom,ff,ener,simulCond,box);
    
    if(simulCond->nb14)
      coulomb14_shift2(atom,ff,ener,simulCond,box);
  }
  else if (simulCond->elecType==SWITCH)
  {
    coulomb_switch(atom,ff,ener,simulCond,box);
    
    if(simulCond->nb14)
      coulomb14_switch(atom,ff,ener,simulCond,box);
  }
  else
  {
    error(201);
  }
  
/* Performing van der Waals interactions */

  if(simulCond->vdwType==NOVDW)
  {
    printf("No non-bonded requested.\n");
  }
  else if (simulCond->vdwType==VFULL)
  {
    vdw_full(atom,ff,ener,simulCond,box);
    if(simulCond->nb14)
      vdw14_full(atom,ff,ener,simulCond,box);
  }
  else if(simulCond->vdwType==VSWITCH)
  {
    vdw_switch(atom,ff,ener,simulCond,box);
    
    if(simulCond->nb14)
      vdw14_switch(atom,ff,ener,simulCond,box);
  }
  else
  {
    error(202);
  }
  
  /* Performing bond terms */
  
  if(ff->nBond>0)
    bond_energy(atom,ff,ener,simulCond,box);
  
  /* Performing angle terms */
  
  if(ff->nAngle>0)
    angle_energy(atom,ff,ener,simulCond,box);
  
  /* Performing Urey-Bradley terms */
  
   if(ff->nUb>0)
    ub_energy(atom,ff,ener,simulCond,box);
   
  /* Performing diherdral terms */
  
  if(ff->nDihedral>0)
    dihedral_energy(atom,ff,ener,simulCond,box);
  
  /* Performing improper terms */
  
  if(ff->nImproper>0)
    improper_energy(atom,ff,ener,simulCond,box);
  
  /* Calculate potential energy */
    
  ener->pot=ener->elec+ener->vdw+ener->bond+
    ener->ang+ener->ub+ener->dihe+ener->impr;

}
