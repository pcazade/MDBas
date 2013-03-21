/**
 * \file energy.c
 * \brief Contains highest level functions for evaluating energy and forces of a system.
 * \author Pierre-Andre Cazade and Florent Hedin
 * \version alpha-branch
 * \date 2012
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "global.h"
#include "energy.h"
#include "elec.h"
#include "vdw.h"
#include "ewald.h"
#include "spme.h"
#include "energy.h"
#include "internal.h"
#include "io.h"
#include "utils.h"
#include "errors.h"

#if (defined TIMING && defined __unix__ && !defined __STRICT_ANSI__)
#define TIMER
#include "timing.h"
#endif

/**
 * \param simulCond Pointer to structure SIMULPARAMS containing parameters of the current simulation.
 * 
 * \brief Set the pointers to function used for non-bonded energy evaluation.
 */
void init_energy_ptrs(CTRL *ctrl)
{
  
  switch(ctrl->elecType)
  {
    case NOELEC:
      ptr_coulomb = &(coulomb_none);
      ptr_coulomb14 = &(coulomb14_none);
    break;
    
    case FULL:
      ptr_coulomb = &(coulomb_none);
      ptr_coulomb14 = &(coulomb14_none);
    break;
    
    case SHIFT1:
      ptr_coulomb = &(coulomb_shift1);
      ptr_coulomb14 = &(coulomb14_shift1);
    break;
    
    case SHIFT2:
      ptr_coulomb = &(coulomb_shift2);
      ptr_coulomb14 = &(coulomb14_shift2);
    break;
    
    case SWITCH:
      ptr_coulomb = &(coulomb_switch);
      ptr_coulomb14 = &(coulomb14_switch);
    break;
    
    default:
      my_error(UNKNOWN_ELEC_ERROR,__FILE__,__LINE__,0);
    break;
  }
  
  switch(ctrl->vdwType)
  {
    case NOVDW:
      ptr_vdw = &(vdw_none);
      ptr_vdw14 = &(vdw14_none);
    break;
    
    case VFULL:
      ptr_vdw = &(vdw_none);
      ptr_vdw14 = &(vdw14_none);
    break;
       
    case VSWITCH:
      ptr_vdw = &(vdw_switch);
      ptr_vdw14 = &(vdw14_switch);
    break;
    
    default:
      my_error(UNKNOWN_VDW_ERROR,__FILE__,__LINE__,0);
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
void energy(CTRL *ctrl,PARAM *param,ENERGY *ener,EWALD *ewald,PBC *box,NEIGH *neigh,
	    BOND bond[],BOND ub[],ANGLE angle[],DIHE dihe[],DIHE impr[],
	    const double x[],const double y[], const double z[],
	    double vx[],double vy[], double vz[],double fx[],double fy[],
	    double fz[],const double q[],const double eps[],const double sig[],
	    const double eps14[],const double sig14[],const int frozen[],
	    const int neighList[],const int neighPair[],const int neighOrder[],
	    const int neighList14[],int **exclList,const int exclPair[])
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
  
  for(i=0;i<param->nAtom;i++)
  {
    fx[i]=0.;
    fy[i]=0.;
    fz[i]=0.;
  }
  
#ifdef TIMER
  update_timer_begin(TIMER_ENERGY_TOT,__func__);
#endif
  
  /* Performing non-bonding terms */
  
  if(ctrl->keyEwald!=0)
  {
    ewald_energy(ctrl,param,ener,ewald,box,x,y,z,fx,fy,fz,q,eps,sig,
		 neighList,neighPair,neighOrder,exclList,exclPair);
    if(ctrl->keyNb14)
      ewald14_energy(param,ener,ewald,box,neigh,x,y,z,fx,fy,fz,q,eps14,sig14,neighList14);
  }
  else
  {
    nonbond_energy(param,ener,box,x,y,z,fx,fy,fz,q,eps,sig,neighList,neighPair,neighOrder);
    if(ctrl->keyNb14)
      nonbond14_energy(param,ener,box,neigh,x,y,z,fx,fy,fz,q,eps14,sig14,neighList14);
  }
  
  /* Performing bond terms */
  
  if(param->nBond>0)
    bond_energy(param,ener,box,bond,x,y,z,fx,fy,fz);
  
  /* Performing angle terms */
  
  if(param->nAngle>0)
    angle_energy(param,ener,box,angle,x,y,z,fx,fy,fz);
  
  /* Performing Urey-Bradley terms */
  
   if(param->nUb>0)
    ub_energy(param,ener,box,ub,x,y,z,fx,fy,fz);
   
  /* Performing diherdral terms */
  
  if(param->nDihedral>0)
    dihedral_energy(param,ener,box,dihe,x,y,z,fx,fy,fz);
  
  /* Performing improper terms */
  
  if(param->nImproper>0)
    improper_energy(param,ener,box,impr,x,y,z,fx,fy,fz);
  
  /* Calculate potential energy */
    
  ener->pot=ener->elec+ener->vdw+ener->bond+
    ener->ang+ener->ub+ener->dihe+ener->impr;
    
  ener->virpot=ener->virbond+ener->virub+ener->virelec+ener->virvdw;
  
#ifdef TIMER
  update_timer_end(TIMER_ENERGY_TOT,__func__);
#endif
  
  if(param->nFrozen>0)
    freeze_atoms(param,vx,vy,vz,fx,fy,fz,frozen);

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
void nonbond_energy(PARAM *param,ENERGY *ener,PBC *box,const double x[],const double y[],
		    const double z[],double fx[],double fy[],double fz[],const double q[],
		    const double eps[],const double sig[],const int neighList[],
		    const int neighPair[],const int neighOrder[])
{
  
  int i,j,k,l,offset;
  double elec=0.,evdw=0.,delec=0.,dvdw=0.,virelec=0.,virvdw=0.;
  double r,r2,rt,fxj,fyj,fzj,fxi,fyi,fzi;
  double qel,veps,vsig;
  double delta[3]/*,stress[6]={0.}*/;
  
#ifdef TIMER
  update_timer_begin(TIMER_ENERGY_NB,__func__);
#endif
  
  offset=0;
  #ifdef _OPENMP
  #pragma omp parallel default(none) shared(atom,param,box,vdw,neigh,ptr_coulomb,ptr_vdw) private(i,j,k,delec,dvdw,r,fx,fy,fz,fxi,fyi,fzi,delta) reduction(+:elec,evdw,virelec,virvdw)
  {
  #pragma omp for schedule(dynamic,1000) nowait
  #endif
    for(l=0;l<param->nAtom;l++)
    {
      i=neighOrder[l];
      
      fxi=0.;
      fyi=0.;
      fzi=0.;
      
      for(k=0;k<neighPair[i];k++)
      {
	j=neighList[offset++];
	
	delta[0]=x[j]-x[i];
	delta[1]=y[j]-y[i];
	delta[2]=z[j]-z[i];
	
	r2=dist(box,delta);
	r=sqrt(r2);
	rt=1./r;
	
	if(r2<=param->cutOff2)
	{
	  qel=q[i]*q[j];
	  elec+=(*ptr_coulomb)(param,&delec,qel,r2,rt);
	  
	  veps=eps[i]*eps[j];
	  vsig=sig[i]+sig[j];
	  evdw+=(*ptr_vdw)(param,&dvdw,veps,vsig,r2,rt);
	  
	  virelec+=delec*r;
	  virvdw+=dvdw*r;
	  
          fxj=(delec+dvdw)*delta[0]*rt;
	  fyj=(delec+dvdw)*delta[1]*rt;
	  fzj=(delec+dvdw)*delta[2]*rt;
	  
	  fxi+=fxj;
	  fyi+=fyj;
	  fzi+=fzj;
	  
	  /*stress[0]-=fx*delta[0];
	  stress[1]-=fy*delta[0];
	  stress[2]-=fz*delta[0];
	  stress[3]-=fy*delta[1];
	  stress[4]-=fz*delta[1];
	  stress[5]-=fz*delta[2];*/
	  
	  #ifdef _OPENMP
	  #pragma omp atomic
	  #endif
	  fx[j]+=-fxj;
	  
	  #ifdef _OPENMP
	  #pragma omp atomic
	  #endif
	  fy[j]+=-fyj;
	  
	  #ifdef _OPENMP
	  #pragma omp atomic
	  #endif
	  fz[j]+=-fzj;
	  
	}
	
      }
      
      #ifdef _OPENMP
      #pragma omp atomic
      #endif
      fx[i]+=fxi;
      
      #ifdef _OPENMP
      #pragma omp atomic
      #endif
      fy[i]+=fyi;
      
      #ifdef _OPENMP
      #pragma omp atomic
      #endif
      fz[i]+=fzi;
      
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
  
#ifdef TIMER
  update_timer_end(TIMER_ENERGY_NB,__func__);
#endif
  
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
void nonbond14_energy(PARAM *param,ENERGY *ener,PBC *box,NEIGH *neigh,
		      const double x[],const double y[],const double z[],double fx[],double fy[],
		      double fz[],const double q[],const double eps[],const double sig[],
		      const int neighList14[])
{
  
  int i,j,k;
  double elec=0.,evdw=0.,delec=0.,dvdw=0.,virelec=0.,virvdw=0.;
  double r,r2,rt,fxj,fyj,fzj;
  double qel,veps,vsig;
  double delta[3]/*,stress[6]={0.}*/;
  
#ifdef TIMER
  update_timer_begin(TIMER_ENERGY_NB14,__func__);
#endif
  
  for(k=0;k<neigh->nPair14;k++)
  {
    i=neighList14[2*k];
    j=neighList14[2*k+1];
    
    delec=0.;
    dvdw=0.;
    
    delta[0]=x[j]-x[i];
    delta[1]=y[j]-y[i];
    delta[2]=z[j]-z[i];
    
    r2=dist(box,delta);
    r=sqrt(r2);
    rt=1./r;
    
    qel=q[i]*q[j];
    elec+=(*ptr_coulomb14)(param,&delec,qel,r2,rt);
    
    veps=eps[i]*eps[j];
    vsig=sig[i]+sig[j];
    evdw+=(*ptr_vdw14)(param,&dvdw,veps,vsig,r2,rt);
    
    virelec+=delec*r;
    virvdw+=dvdw*r;
      
    fxj=(delec+dvdw)*delta[0]*rt;
    fyj=(delec+dvdw)*delta[1]*rt;
    fzj=(delec+dvdw)*delta[2]*rt;
    
    /*stress[0]-=fx*delta[0];
    stress[1]-=fy*delta[0];
    stress[2]-=fz*delta[0];
    stress[3]-=fy*delta[1];
    stress[4]-=fz*delta[1];
    stress[5]-=fz*delta[2];*/
  
    fx[i]+=fxj;
    fy[i]+=fyj;
    fz[i]+=fzj;
  
    fx[j]+=-fxj;
    fy[j]+=-fyj;
    fz[j]+=-fzj;
      
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
  
#ifdef TIMER
  update_timer_end(TIMER_ENERGY_NB14,__func__);
#endif

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
void ewald_energy(CTRL *ctrl,PARAM *param,ENERGY *ener,EWALD *ewald,PBC *box,const double x[],
		  const double y[],const double z[],double fx[],double fy[],
		  double fz[],const double q[],const double eps[],const double sig[],
		  const int neighList[],const int neighPair[],const int neighOrder[],
		  int **exclList,const int exclPair[])
{
  
  int i,j,k,l,offset;
  double eEwaldRec=0.,virEwaldRec=0.,eEwaldDir=0.,dEwaldDir=0.,virEwaldDir=0.;
  double eEwaldCorr=0.,dEwaldCorr=0.,virEwaldCorr=0.;
  double evdw=0.,dvdw=0.,virvdw=0.;
  double r,r2,rt,fxj,fyj,fzj,fxi,fyi,fzi;
  double qel,veps,vsig;
  double delta[3],stress1[6]={0.};
  
#ifdef TIMER
  update_timer_begin(TIMER_ENERGY_NB,__func__);
#endif
  
  if(ctrl->keyEwald==1)
    eEwaldRec=ewald_rec(param,ewald,box,x,y,z,fx,fy,fz,q,stress1,&virEwaldRec);
  else if(ctrl->keyEwald==2)
    eEwaldRec=spme_energy(param,ewald,box,x,y,z,fx,fy,fz,q,stress1,&virEwaldRec);
  
  offset=0;
  #ifdef _OPENMP
  #pragma omp parallel default(none) shared(atom,param,box,vdw,neigh,ptr_coulomb,ptr_vdw) private(i,j,k,delec,dvdw,r,fx,fy,fz,fxi,fyi,fzi,delta) reduction(+:elec,evdw,virelec,virvdw)
  {
  #pragma omp for schedule(dynamic,1000) nowait
  #endif
    for(l=0;l<param->nAtom;l++)
    {
      i=neighOrder[l];
      
      fxi=0.;
      fyi=0.;
      fzi=0.;
      
      for(k=0;k<neighPair[i];k++)
      {
	j=neighList[offset++];
	
	delta[0]=x[j]-x[i];
	delta[1]=y[j]-y[i];
	delta[2]=z[j]-z[i];
	
	r2=dist(box,delta);
	r=sqrt(r2);
	rt=1./r;
	
	if(r2<=param->cutOff2)
	{ 
	  qel=param->chargeConst*q[i]*q[j];
	  eEwaldDir+=ewald_dir(ewald,&dEwaldDir,qel,r,rt);
	  
	  veps=eps[i]*eps[j];
	  vsig=sig[i]+sig[j];
	  evdw+=(*ptr_vdw)(param,&dvdw,veps,vsig,r2,rt);
	  
	  virEwaldDir+=dEwaldDir*r;
	  virvdw+=dvdw*r;
	  
          fxj=(dEwaldDir+dvdw)*delta[0]*rt;
	  fyj=(dEwaldDir+dvdw)*delta[1]*rt;
	  fzj=(dEwaldDir+dvdw)*delta[2]*rt;
	  
	  fxi+=fxj;
	  fyi+=fyj;
	  fzi+=fzj;
	  
	  /*stress[0]-=fx*delta[0];
	  stress[1]-=fy*delta[0];
	  stress[2]-=fz*delta[0];
	  stress[3]-=fy*delta[1];
	  stress[4]-=fz*delta[1];
	  stress[5]-=fz*delta[2];*/
	  
	  #ifdef _OPENMP
	  #pragma omp atomic
	  #endif
	  fx[j]+=-fxj;
	  
	  #ifdef _OPENMP
	  #pragma omp atomic
	  #endif
	  fy[j]+=-fyj;
	  
	  #ifdef _OPENMP
	  #pragma omp atomic
	  #endif
	  fz[j]+=-fzj;
	  
	}
	
      }
      
      #ifdef _OPENMP
      #pragma omp atomic
      #endif
      fx[i]+=fxi;
      
      #ifdef _OPENMP
      #pragma omp atomic
      #endif
      fy[i]+=fyi;
      
      #ifdef _OPENMP
      #pragma omp atomic
      #endif
      fz[i]+=fzi;
      
    } //end of parallel for
    
  #ifdef _OPENMP
  } //end of parallel region
  #endif
  
  for(i=0;i<param->nAtom;i++)
  {
    
    fxi=0.;
    fyi=0.;
    fzi=0.;
    
    for(k=0;k<exclPair[i];k++)
    {
      j=exclList[i][k];
      
      if(j>i)
      {
      
	delta[0]=x[j]-x[i];
	delta[1]=y[j]-y[i];
	delta[2]=z[j]-z[i];
	
	r2=dist(box,delta);
	r=sqrt(r2);
	rt=1./r;
	
	qel=param->chargeConst*q[i]*q[j];
	eEwaldCorr+=ewald_corr(ewald,&dEwaldCorr,qel,r,rt);
	
	virEwaldCorr+=dEwaldCorr*r;
	
	fxj=dEwaldCorr*delta[0]*rt;
	fyj=dEwaldCorr*delta[1]*rt;
	fzj=dEwaldCorr*delta[2]*rt;
	
	fxi+=fxj;
	fyi+=fyj;
	fzi+=fzj;
	
	/*stress[0]-=fx*delta[0];
	stress[1]-=fy*delta[0];
	stress[2]-=fz*delta[0];
	stress[3]-=fy*delta[1];
	stress[4]-=fz*delta[1];
	stress[5]-=fz*delta[2];*/
	
	fx[j]+=-fxj;
	fy[j]+=-fyj;
	fz[j]+=-fzj;
	
      }
      
    }
    
    fx[i]+=fxi;
    fy[i]+=fyi;
    fz[i]+=fzi;
    
  }
  
  ener->elec+=eEwaldDir+eEwaldRec+eEwaldCorr;
  ener->vdw+=evdw;
  
  ener->virelec+=virEwaldDir+virEwaldRec+virEwaldCorr;
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
  
#ifdef TIMER
  update_timer_end(TIMER_ENERGY_NB,__func__);
#endif
  
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
void ewald14_energy(PARAM *param,ENERGY *ener,EWALD *ewald,PBC *box,NEIGH *neigh,const double x[],
		    const double y[],const double z[],double fx[],double fy[],
		    double fz[],const double q[],const double eps[],const double sig[],
		    const int neighList14[])
{
  
  int i,j,k;
  double eEwaldDir=0.,dEwaldDir=0.,virEwaldDir=0.;
  double eEwaldCorr=0.,dEwaldCorr=0.,virEwaldCorr=0.;
  double evdw=0.,dvdw=0.,virvdw=0.;
  double r,r2,rt,fxj,fyj,fzj;
  double qel,veps,vsig;
  double delta[3]/*,stress[6]={0.}*/;
  
#ifdef TIMER
  update_timer_begin(TIMER_ENERGY_NB14,__func__);
#endif
  
  for(k=0;k<neigh->nPair14;k++)
  {
    i=neighList14[2*k];
    j=neighList14[2*k+1];
    
    dEwaldCorr=0.;
    dvdw=0.;
    
    delta[0]=x[j]-x[i];
    delta[1]=y[j]-y[i];
    delta[2]=z[j]-z[i];
    
    r2=dist(box,delta);
    r=sqrt(r2);
    rt=1./r;
    
    qel=param->chargeConst*q[i]*q[j];
    
    eEwaldDir+=ewald_dir14(param,ewald,&dEwaldDir,qel,r,rt);
    eEwaldCorr+=ewald_corr14(param,ewald,&dEwaldCorr,qel,r,rt);
    
    veps=eps[i]*eps[j];
    vsig=sig[i]+sig[j];
    evdw+=(*ptr_vdw14)(param,&dvdw,veps,vsig,r2,rt);
    
    virEwaldDir+=dEwaldDir*r;
    virEwaldCorr+=dEwaldCorr*r;
    virvdw+=dvdw*r;
      
    fxj=(dEwaldDir+dEwaldCorr+dvdw)*delta[0]*rt;
    fyj=(dEwaldDir+dEwaldCorr+dvdw)*delta[1]*rt;
    fzj=(dEwaldDir+dEwaldCorr+dvdw)*delta[2]*rt;
    
    /*stress[0]-=fx*delta[0];
    stress[1]-=fy*delta[0];
    stress[2]-=fz*delta[0];
    stress[3]-=fy*delta[1];
    stress[4]-=fz*delta[1];
    stress[5]-=fz*delta[2];*/
  
    fx[i]+=fxj;
    fy[i]+=fyj;
    fz[i]+=fzj;
  
    fx[j]+=-fxj;
    fy[j]+=-fyj;
    fz[j]+=-fzj;
      
  }
  
  ener->elec+=eEwaldDir+eEwaldCorr;
  ener->vdw+=evdw;
  
  ener->virelec+=virEwaldDir+virEwaldCorr;
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
  
#ifdef TIMER
  update_timer_end(TIMER_ENERGY_NB14,__func__);
#endif

}
