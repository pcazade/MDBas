/**
 * \file energy.c
 * \brief Contains highest level functions for evaluating energy and forces of a system.
 * \author Pierre-Andre Cazade and Florent Hedin
 * \version alpha-branch
 * \date 2012
 */

#include <stdio.h>
#include <stdlib.h>

#include "global.h"
#include "energy.h"
#include "elec.h"
#include "vdw.h"
#include "internal.h"
#include "io.h"
#include "utils.h"

/**
 * \param simulCond Pointer to structure SIMULPARAMS containing parameters of the current simulation.
 * 
 * \brief Set the pointers to function used for non-bonded energy evaluation.
 */
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

/**
 * \param atom Array of structure ATOM (coordinates, forces, etc...).
 * \param ff Pointer to structure FORCEFIELD containing forcefield parameters.
 * \param ener Pointer to structure ENERGY containing values of the different energies.
 * \param simulCond Pointer to structure SIMULPARAMS containing parameters of the current simulation.
 * \param box Pointer to structure PBC containing Periodic Boundaries Conditions parameters.
 * 
 * \brief Main energy function, collecting total energy by calling the required subfunctions.
 */
void energy(ATOM atom[], FORCEFIELD *ff, ENERGY *ener, SIMULPARAMS *simulCond, PBC *box)
{
  
  int i;
  
  ener->elec=0.;
  ener->vdw=0.;
  ener->bond=0.;
  ener->ub=0.;
  ener->ang=0.;
  ener->dihe=0.;
  ener->impr=0.;
  
  ener->virelec=0.;
  ener->virvdw=0.;
  ener->virbond=0.;
  ener->virub=0.;
  ener->virpot=0.;
  
  box->stress1=0.;
  box->stress2=0.;
  box->stress3=0.;
  box->stress4=0.;
  box->stress5=0.;
  box->stress6=0.;
  box->stress7=0.;
  box->stress8=0.;
  box->stress9=0.;
  
  for(i=0;i<simulCond->natom;i++)
  {
    atom[i].fx=0.;
    atom[i].fy=0.;
    atom[i].fz=0.;
  }
  
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
    
  ener.virpot=ener.virbond+ener.virub+ener.virelec+ener.virvdw;

}

/**
 * \param atom Array of structure ATOM (coordinates, forces, etc...).
 * \param ff Pointer to structure FORCEFIELD containing forcefield parameters.
 * \param ener Pointer to structure ENERGY containing values of the different energies.
 * \param simulCond Pointer to structure SIMULPARAMS containing parameters of the current simulation.
 * \param box Pointer to structure PBC containing Periodic Boundaries Conditions parameters.
 * 
 * \brief Energy function collecting all terms of non-bonded energies.
 */
void nonbond_energy(ATOM atom[] ,FORCEFIELD *ff, ENERGY *ener, SIMULPARAMS *simulCond, PBC *box)
{
  
  int i,j,k;
  double elec=0.,evdw=0.,delec=0.,dvdw=0.,virelec=0.,virvdw=0.;
  double r,fx,fy,fz,fxi,fyi,fzi;
  double delta[3],stress[6]={0.};
  
  #ifdef _OPENMP

  #pragma omp parallel default(none) shared(atom,ff,simulCond,box,coulomb,vdw) private(i,j,k,delec,dvdw,r,fx,fy,fz,fxi,fyi,fzi,delta) reduction(+:elec,evdw,virelec,virvdw)
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
	  
	  virelec+=delec*r;
	  virvdw+=dvdw*r;
	  
          fx=(delec+dvdw)*delta[0]/r;
	  fy=(delec+dvdw)*delta[1]/r;
	  fz=(delec+dvdw)*delta[2]/r;
	  
	  fxi+=fx;
	  fyi+=fy;
	  fzi+=fz;
	  
	  /*stress[0]-=fx*delta[0];
	  stress[1]-=fy*delta[0];
	  stress[2]-=fz*delta[0];
	  stress[3]-=fy*delta[1];
	  stress[4]-=fz*delta[1];
	  stress[5]-=fz*delta[2];*/
	  
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
  
  ener->virelec+=virelec;
  ener->virvdw+=virvdw;
  
  /*box->stress1+=stress[0];
  box->stress2+=stress[1];
  box->stress3+=stress[2];
  box->stress4+=stress[1];
  box->stress5+=stress[3];
  box->stress6+=stress[4];
  box->stress7+=stress[2];
  box->stress8+=stress[4];
  box->stress9+=stress[5];*/
  
}

/**
 * \param atom Array of structure ATOM (coordinates, forces, etc...).
 * \param ff Pointer to structure FORCEFIELD containing forcefield parameters.
 * \param ener Pointer to structure ENERGY containing values of the different energies.
 * \param simulCond Pointer to structure SIMULPARAMS containing parameters of the current simulation.
 * \param box Pointer to structure PBC containing Periodic Boundaries Conditions parameters.
 * 
 * \brief Energy function collecting all terms of 1-4 non-bonded energies.
 */
void nonbond14_energy(ATOM atom[] ,FORCEFIELD *ff, ENERGY *ener, SIMULPARAMS *simulCond, PBC *box)
{
  
  int i,j,k;
  double elec=0.,evdw=0.,delec=0.,dvdw=0.,virelec=0.,virvdw=0.;
  double r,fx,fy,fz;
  double delta[3],stress[6]={0.};
  
  for(k=0;k<ff->npr14;k++)
  {
    i=ff->ver14[k][0];
    j=ff->ver14[k][1];
    
    delec=0.;
    dvdw=0.;
    
    r=distance(i,j,atom,delta,simulCond,box);
      
    elec+=(*coulomb14)(atom,ff,/*ener,*/simulCond,box,i,j,r,&delec);
    evdw+=(*vdw14)(atom,ff,/*ener,*/simulCond,box,i,j,r,&dvdw);
    
    virelec+=delec*r;
    virvdw+=dvdw*r;
      
    fx=(delec+dvdw)*delta[0]/r;
    fy=(delec+dvdw)*delta[1]/r;
    fz=(delec+dvdw)*delta[2]/r;
    
    /*stress[0]-=fx*delta[0];
    stress[1]-=fy*delta[0];
    stress[2]-=fz*delta[0];
    stress[3]-=fy*delta[1];
    stress[4]-=fz*delta[1];
    stress[5]-=fz*delta[2];*/
  
    atom[i].fx+=fx;
    atom[i].fy+=fy;
    atom[i].fz+=fz;
  
    atom[j].fx+=-fx;
    atom[j].fy+=-fy;
    atom[j].fz+=-fz;
      
  }
  
  ener->elec+=elec;
  ener->vdw+=evdw;
  
  ener->virelec+=virelec;
  ener->virvdw+=virvdw;
  
  /*box->stress1+=stress[0];
  box->stress2+=stress[1];
  box->stress3+=stress[2];
  box->stress4+=stress[1];
  box->stress5+=stress[3];
  box->stress6+=stress[4];
  box->stress7+=stress[2];
  box->stress8+=stress[4];
  box->stress9+=stress[5];*/

}
