#include <stdio.h>
#include <stdlib.h>
#include "global.h"
#include "elec.h"
#include "vdw.h"
#include "internal.h"
#include "io.h"
#include "list.h"

void energy(ATOM *atom,FORCEFIELD *ff,ENERGYFORCE *enerFor,SIMULPARAMS *simulCond)
{
  
  int i;
  
/* Performing elctrostatic interactions */  
//   if(lqpoly)
//       update_charges(ff);
  if(simulCond->firstener==1)
  {
    printf("Building exclusion and Verlet lists.\n");
    exclude_list(simulCond,atom,ff);
    verlet_list(simulCond,atom,ff);
    
    simulCond->firstener=0;
    
  }
  else if(simulCond->step>0)
  {
    printf("Updating Verlet list.\n");
    verlet_list_update(simulCond,atom,ff);
  }
  
  enerFor->energyElec=0.;
  enerFor->energyVdw=0.;
  enerFor->energyBond=0.;
  enerFor->energyUb=0.;
  enerFor->energyAng=0.;
  enerFor->energyDih=0.;
  enerFor->energyImpr=0.;
  
  for(i=0;i<atom->natom;i++)
  {
    atom->fx[i]=0.;
    atom->fy[i]=0.;
    atom->fz[i]=0.;
  }
  
  /* Performing electrostatic interactions */
  
  if(simulCond->elecType==NOELEC)
  {
    printf("No electrostatics requested\n");
  }
  else if (simulCond->elecType==FULL)
  {
    coulomb_full(atom,ff,enerFor,simulCond);
    if(simulCond->nb14)
      coulomb14_full(atom,ff,enerFor,simulCond);
  }
  else if(simulCond->elecType==SHIFT1)
  {
    coulomb_shift1(atom,ff,enerFor,simulCond);
    
    if(simulCond->nb14)
      coulomb14_shift1(atom,ff,enerFor,simulCond);
  }
  else if(simulCond->elecType==SHIFT2)
  {
    coulomb_shift2(atom,ff,enerFor,simulCond);
    
    if(simulCond->nb14)
      coulomb14_shift2(atom,ff,enerFor,simulCond);
  }
  else if (simulCond->elecType==SWITCH)
  {
    coulomb_switch(atom,ff,enerFor,simulCond);
    
    if(simulCond->nb14)
      coulomb14_switch(atom,ff,enerFor,simulCond);
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
    vdw_full(atom,ff,enerFor,simulCond);
    if(simulCond->nb14)
      vdw14_full(atom,ff,enerFor,simulCond);
  }
  else if(simulCond->vdwType==VSWITCH)
  {
    vdw_switch(atom,ff,enerFor,simulCond);
    
    if(simulCond->nb14)
      vdw14_switch(atom,ff,enerFor,simulCond);
  }
  else
  {
    error(202);
  }
  
  /* Performing bond terms */
  
  if(ff->nBond>0)
    bond_energy(atom,ff,enerFor,simulCond);
  
  /* Performing angle terms */
  
  if(ff->nAngle>0)
    angle_energy(atom,ff,enerFor,simulCond);
  
  /* Performing Urey-Bradley terms */
  
   if(ff->nUb>0)
    ub_energy(atom,ff,enerFor,simulCond);
   
  /* Performing diherdral terms */
  
  if(ff->nDihedral>0)
    dihedral_energy(atom,ff,enerFor,simulCond);
  
  /* Performing improper terms */
  
  if(ff->nImproper>0)
    improper_energy(atom,ff,enerFor,simulCond);
  
  /* Calculate potential energy */
    
  enerFor->energyPot=enerFor->energyElec+enerFor->energyVdw+enerFor->energyBond+
    enerFor->energyAng+enerFor->energyUb+enerFor->energyDih+enerFor->energyImpr;

}
