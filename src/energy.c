#include <stdio.h>
#include <stdlib.h>
#include "global.h"
#include "energy.h"
#include "elec.h"
#include "vdw.h"
#include "internal.h"
#include "io.h"
#include "utils.h"

void init_energy_ptrs(SIMULPARAMS *simulCond)
{
  
  switch(simulCond->elecType)
  {
    case NOELEC:
      coulomb = &(coulomb_none);
      coulomb14 = &(coulomb14_none);
    break;
    
    case FULL:
      coulomb = &(coulomb_none);
      coulomb14 = &(coulomb14_none);
    break;
    
    case SHIFT1:
      coulomb = &(coulomb_shift1);
      coulomb14 = &(coulomb14_shift1);
    break;
    
    case SHIFT2:
      coulomb = &(coulomb_shift2);
      coulomb14 = &(coulomb14_shift2);
    break;
    
    case SWITCH:
      coulomb = &(coulomb_switch);
      coulomb14 = &(coulomb14_switch);
    break;
    
    default:
      error(201);
    break;
  }
  
  switch(simulCond->vdwType)
  {
    case NOVDW:
      vdw = &(vdw_none);
      vdw14 = &(vdw14_none);
    break;
    
    case VFULL:
      vdw = &(vdw_none);
      vdw14 = &(vdw14_none);
    break;
       
    case VSWITCH:
      vdw = &(vdw_switch);
      vdw14 = &(vdw14_switch);
    break;
    
    default:
      error(202);
    break;
  }

}

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
  
  /*if(simulCond->elecType==NOELEC)
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
  }*/
  
/* Performing van der Waals interactions */

  /*if(simulCond->vdwType==NOVDW)
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
  }*/
  
  /* Performing non-bonding terms */
  
  nonbond_energy(atom,ff,ener,simulCond,box);
  if(simulCond->nb14)
    nonbond14_energy(atom,ff,ener,simulCond,box);
  
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

void nonbond_energy(ATOM *atom,FORCEFIELD *ff,ENERGY *ener,SIMULPARAMS *simulCond,PBC *box)
{
  
  int i,j,k;
  double elec=0.,evdw=0.,delec=0.,dvdw=0.;
  double r,fx,fy,fz,fxi,fyi,fzi;
  double delta[3];
  
  #ifdef _OPENMP

  #pragma omp parallel default(none) shared(atom,ff,simulCond,box,coulomb,vdw) private(i,j,k,delec,dvdw,r,fx,fy,fz,fxi,fyi,fzi,delta) reduction(+:elec,evdw)
  {
  #pragma omp for schedule(dynamic) nowait
  #endif
    for(i=0;i<simulCond->natom;i++)
    {
      fxi=0.;
      fyi=0.;
      fzi=0.;
      
      for(k=0;k<ff->verPair[i];k++)
      {
	j=ff->verList[i][k];
	
	r=distance(i,j,atom,delta,simulCond,box);
	
	if(r<=simulCond->cutoff)
	{
	  elec+=(*coulomb)(atom,ff,/*ener,*/simulCond,box,i,j,r,&delec);
	  evdw+=(*vdw)(atom,ff,/*ener,*/simulCond,box,i,j,r,&dvdw);
	  
          fx=(delec+dvdw)*delta[0]/r;
	  fy=(delec+dvdw)*delta[1]/r;
	  fz=(delec+dvdw)*delta[2]/r;
	  
	  fxi+=fx;
	  fyi+=fy;
	  fzi+=fz;
	  
	  #ifdef _OPENMP
	  #pragma omp atomic
	  #endif
	  atom[j].fx+=-fx;
	  
	  #ifdef _OPENMP
	  #pragma omp atomic
	  #endif
	  atom[j].fy+=-fy;
	  
	  #ifdef _OPENMP
	  #pragma omp atomic
	  #endif
	  atom[j].fz+=-fz;
	  
	}
	
      }
      
      #ifdef _OPENMP
      #pragma omp atomic
      #endif
      atom[i].fx+=fxi;
      
      #ifdef _OPENMP
      #pragma omp atomic
      #endif
      atom[i].fy+=fyi;
      
      #ifdef _OPENMP
      #pragma omp atomic
      #endif
      atom[i].fz+=fzi;
      
    } //end of parallel for
    
  #ifdef _OPENMP
  } //end of parallel region
  #endif
  
  ener->elec+=elec;
  ener->vdw+=evdw;
  
}

void nonbond14_energy(ATOM *atom,FORCEFIELD *ff,ENERGY *ener,SIMULPARAMS *simulCond,PBC *box)
{
  
  int i,j,k;
  double elec=0.,evdw=0.,delec=0.,dvdw=0.;
  double r,fx,fy,fz;
  double delta[3];
  
  for(k=0;k<ff->npr14;k++)
  {
    i=ff->ver14[k][0];
    j=ff->ver14[k][1];
    
    delec=0.;
    dvdw=0.;
    
    r=distance(i,j,atom,delta,simulCond,box);
      
    elec+=(*coulomb14)(atom,ff,/*ener,*/simulCond,box,i,j,r,&delec);
    evdw+=(*vdw14)(atom,ff,/*ener,*/simulCond,box,i,j,r,&dvdw);
      
    fx=(delec+dvdw)*delta[0]/r;
    fy=(delec+dvdw)*delta[1]/r;
    fz=(delec+dvdw)*delta[2]/r;
  
    atom[i].fx+=fx;
    atom[i].fy+=fy;
    atom[i].fz+=fz;
  
    atom[j].fx+=-fx;
    atom[j].fy+=-fy;
    atom[j].fz+=-fz;
      
  }
  
  ener->elec+=elec;
  ener->vdw+=evdw;

}
