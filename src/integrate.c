/*
 * Copyright (c) 2013 Pierre-Andre Cazade
 * Copyright (c) 2013 Florent hedin
 * 
 * This file is part of MDBas.
 *
 * MDBas is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MDBas is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with MDBas.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * \file integrate.c
 * \brief Contains functions performing the integration
 * \author Pierre-Andre Cazade and Florent Hedin
 * \version alpha-branch
 * \date 2012
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "global.h"
#include "utils.h"
#include "shake.h"
#include "integrate.h"
#include "memory.h"

#if (defined TIMING && defined __unix__ && !defined __STRICT_ANSI__)
#define TIMER
#include "timing.h"
#endif

#ifdef _OPENMP
#undef _OPENMP
#endif

/** Pointer to the output file. **/
extern FILE *outFile;

static double *ddx,*ddy,*ddz;
static double *xo,*yo,*zo;
static double *xt,*yt,*zt;
static double *vxo,*vyo,*vzo;
static double *vxu,*vyu,*vzu;

void integrators_allocate_arrays(CTRL *ctrl, PARAM *param, PARALLEL *parallel)
{
  ddx=ddy=ddz=NULL;
  xo=yo=zo=xt=yt=zt=vxo=vyo=vzo=vxu=vyu=vzu=NULL;

  if(ctrl->integrator == LEAPFROG)
  {
    if (ctrl->ens == NVE)
    {
      xo=(double*)my_malloc(parallel->nAtProc*sizeof(*xo));
      yo=(double*)my_malloc(parallel->nAtProc*sizeof(*yo));
      zo=(double*)my_malloc(parallel->nAtProc*sizeof(*zo));
      vxu=(double*)my_malloc(param->nAtom*sizeof(*vxu));
      vyu=(double*)my_malloc(param->nAtom*sizeof(*vyu));
      vzu=(double*)my_malloc(param->nAtom*sizeof(*vzu));
      if(parallel->nCtProc>0)
      {
	ddx=(double*)my_malloc(parallel->nCtProc*sizeof(*ddx));
	ddy=(double*)my_malloc(parallel->nCtProc*sizeof(*ddy));
	ddz=(double*)my_malloc(parallel->nCtProc*sizeof(*ddz));
      }
    }
    else
    {
      xo=(double*)my_malloc(parallel->nAtProc*sizeof(*xo));
      yo=(double*)my_malloc(parallel->nAtProc*sizeof(*yo));
      zo=(double*)my_malloc(parallel->nAtProc*sizeof(*zo));
      vxo=(double*)my_malloc(parallel->nAtProc*sizeof(*vxo));
      vyo=(double*)my_malloc(parallel->nAtProc*sizeof(*vyo));
      vzo=(double*)my_malloc(parallel->nAtProc*sizeof(*vzo));
      vxu=(double*)my_malloc(param->nAtom*sizeof(*vxu));
      vyu=(double*)my_malloc(param->nAtom*sizeof(*vyu));
      vzu=(double*)my_malloc(param->nAtom*sizeof(*vzu));
      if(parallel->nCtProc>0)
      {
	ddx=(double*)my_malloc(parallel->nCtProc*sizeof(*ddx));
	ddy=(double*)my_malloc(parallel->nCtProc*sizeof(*ddy));
	ddz=(double*)my_malloc(parallel->nCtProc*sizeof(*ddz));
	
        xt=(double*)my_malloc(param->nAtom*sizeof(*xt));
        yt=(double*)my_malloc(param->nAtom*sizeof(*yt));
        zt=(double*)my_malloc(param->nAtom*sizeof(*zt));
      }
    }
  }
  else if (ctrl->integrator == VELOCITY)
  {
    if (ctrl->ens == NVE || ctrl->ens == NVT_B || ctrl->ens == NVT_H)
    {
        if(parallel->nCtProc>0)
	{
	  ddx=(double*)my_malloc(parallel->nCtProc*sizeof(*ddx));
	  ddy=(double*)my_malloc(parallel->nCtProc*sizeof(*ddy));
	  ddz=(double*)my_malloc(parallel->nCtProc*sizeof(*ddz));
	}
    }
    else
    {
      xo=(double*)my_malloc(parallel->nAtProc*sizeof(*xo));
      yo=(double*)my_malloc(parallel->nAtProc*sizeof(*yo));
      zo=(double*)my_malloc(parallel->nAtProc*sizeof(*zo));
      vxo=(double*)my_malloc(parallel->nAtProc*sizeof(*vxo));
      vyo=(double*)my_malloc(parallel->nAtProc*sizeof(*vyo));
      vzo=(double*)my_malloc(parallel->nAtProc*sizeof(*vzo));
      if(parallel->nCtProc>0)
      {
	ddx=(double*)my_malloc(parallel->nCtProc*sizeof(*ddx));
	ddy=(double*)my_malloc(parallel->nCtProc*sizeof(*ddy));
	ddz=(double*)my_malloc(parallel->nCtProc*sizeof(*ddz));
      }
    }
  }
}

void integrators_free_arrays(CTRL *ctrl, PARAM *param, PARALLEL *parallel)
{
  if(ctrl->integrator == LEAPFROG)
  {
    if (ctrl->ens == NVE)
    {
      free(xo);
      free(yo);
      free(zo);
      free(vxu);
      free(vyu);
      free(vzu);
      if(parallel->nCtProc>0)
      {
	free(ddx);
	free(ddy);
	free(ddz);
      }
    }
    else
    {
      free(xo);
      free(yo);
      free(zo);
      free(vxo);
      free(vyo);
      free(vzo);
      free(vxu);
      free(vyu);
      free(vzu);
      if(parallel->nCtProc>0)
      {
	free(ddx);
	free(ddy);
	free(ddz);
	
        free(xt);
        free(yt);
        free(zt);
      }
    }
  }
  else if (ctrl->integrator == VELOCITY)
  {
    if (ctrl->ens == NVE || ctrl->ens == NVT_B || ctrl->ens == NVT_H)
    {
      if(parallel->nCtProc>0)
      {
	free(ddx);
	free(ddy);
	free(ddz);
      }
    }
    else
    {
      free(xo);
      free(yo);
      free(zo);
      free(vxo);
      free(vyo);
      free(vzo);
      if(parallel->nCtProc>0)
      {
	free(ddx);
	free(ddy);
	free(ddz);
      }
    }
  }

}

void lf_integrate(CTRL *ctrl,PARAM *param,ENERGY *ener,PBC *box,
		  BATH *bath,CONSTRAINT constList[],
		  double *x,double *y,double *z,
		  double *vx,double *vy,double *vz,
		  double *fx,double *fy,double *fz,
		  double *mass,double *rmass,int *parallel->tConst)
{
  
#ifdef TIMER
  update_timer_begin(TIMER_INTEGRATE,__func__);
#endif
  
  switch (ctrl->ens)
  {
    case NVE:
      lf_nve(param,ener,box,constList,x,y,z,vx,vy,vz,fx,fy,fz,mass,rmass,parallel->tConst);
      break;
    case NVT_B:
      lf_nvt_b(param,ener,box,bath,constList,x,y,z,vx,vy,vz,fx,fy,fz,mass,rmass,parallel->tConst);
      break;
    case NPT_B:
      lf_npt_b(param,ener,box,bath,constList,x,y,z,vx,vy,vz,fx,fy,fz,mass,rmass,parallel->tConst);
      break;
    case NVT_H:
      lf_nvt_h(param,ener,box,bath,constList,x,y,z,vx,vy,vz,fx,fy,fz,mass,rmass,parallel->tConst);
      break;
    case NPT_H:
      lf_npt_h(param,ener,box,bath,constList,x,y,z,vx,vy,vz,fx,fy,fz,mass,rmass,parallel->tConst);
      break;
    default:
      lf_nve(param,ener,box,constList,x,y,z,vx,vy,vz,fx,fy,fz,mass,rmass,parallel->tConst);
      break;
  }
  
#ifdef TIMER
  update_timer_end(TIMER_INTEGRATE,__func__);
#endif
  
}

void lf_nve(PARAM *param,ENERGY *ener,PBC *box,CONSTRAINT constList[],
	    double *x,double *y,double *z,
	    double *vx,double *vy,double *vz,
	    double *fx,double *fy,double *fz,
	    double *mass,double *rmass,int *nAparallel->tConst)
{
  int i,ia,ib,l;
  double virshake=0.,stress[6]={0.},stresk[6]={0.};
  
  if(parallel->nCtProc>0)
  {
    l=0;
    #ifdef _OPENMP
    #pragma omp parallel for default(none) shared(param,atom,constList,dd) private(i,ia,ib)
    #endif
    for(i=parallel->fConst;i<parallel->lConst;i++)
    {
      ia=constList[i].a;
      ib=constList[i].b;
      
      ddx[l]=x[ib]-x[ia];
      ddy[l]=y[ib]-y[ia];
      ddz[l]=z[ib]-z[ia];
      
      l++;
    }
    
    image_array(box,ddx,ddy,ddz,parallel->tConst);
    
    l=0;
    #ifdef _OPENMP
    #pragma omp parallel for default(none) shared(param,xo,yo,zo,atom) private(i)
    #endif
    for(i=parallel->fAtom;i<parallel->lAtom;i++)
    {
      
// Store old coordinates.
      
      xo[l]=x[i];
      yo[l]=y[i];
      zo[l]=z[i];
      
      l++;
      
    }
    
  }

// move atoms by leapfrog algorithm
  #ifdef _OPENMP
  #pragma omp parallel for default(none) shared(vxu,vyu,vzu,param,atom) private(i)
  #endif
  for(i=parallel->fAtom;i<parallel->lAtom;i++)
  {
    
// update velocities
    
    vxu[i]=vx[i]+param->timeStep*fx[i]*rmass[i];
    vyu[i]=vy[i]+param->timeStep*fy[i]*rmass[i];
    vzu[i]=vz[i]+param->timeStep*fz[i]*rmass[i];
    
// update positions
    
    x[i]+=param->timeStep*vxu[i];
    y[i]+=param->timeStep*vyu[i];
    z[i]+=param->timeStep*vzu[i];
    
  }
  
  if(param->nConst>0)
  {
    
    if(parallel->nProc>1)
    {
      update_double_para(param,parallel,x,buffer);
      update_double_para(param,parallel,y,buffer);
      update_double_para(param,parallel,z,buffer);
    }
    
// Apply constraint with Shake algorithm.

    lf_shake(param,box,constList,x,y,z,ddx,ddy,ddz,rmass,nAparallel->tConst,stress,&virshake);
    
    l=0;
    #ifdef _OPENMP
    #pragma omp parallel for default(none) shared(vxu,vyu,vzu,xo,yo,zo,param,atom) private(i)
    #endif
    for(i=parallel->fAtom;i<parallel->lAtom;i++)
    {
        
// Corrected velocities
    
      vxu[i]=(x[i]-xo[l])*param->rTimeStep;
      vyu[i]=(y[i]-yo[l])*param->rTimeStep;
      vzu[i]=(z[i]-zo[l])*param->rTimeStep;
    
// Corrected Forces
    
      fx[i]=(vxu[i]-vx[i])*mass[i]*param->rTimeStep;
      fy[i]=(vyu[i]-vy[i])*mass[i]*param->rTimeStep;
      fz[i]=(vzu[i]-vz[i])*mass[i]*param->rTimeStep;
    
    }
  }
  
// calculate full timestep velocity
  #ifdef _OPENMP
  #pragma omp parallel for default(none) shared(vxu,vyu,vzu,param,atom) private(i)
  #endif
  for(i=parallel->fAtom;i<parallel->lAtom;i++)
  {
    
    vx[i]=0.5*(vx[i]+vxu[i]);
    vy[i]=0.5*(vy[i]+vyu[i]);
    vz[i]=0.5*(vz[i]+vzu[i]);
    
  }
  
// calculate kinetic energy
  
  ener->kin=kinetic(param,vx,vy,vz,mass);
  
  ener->virshake=virshake;
  
//   stress_kinetic(atom,simulCond,stresk);
  stress_kinetic(param,vx,vy,vz,mass,stresk);
  
  box->stress1+=stress[0]+stresk[0];
  box->stress2+=stress[1]+stresk[1];
  box->stress3+=stress[2]+stresk[2];
  box->stress4+=stress[1]+stresk[1];
  box->stress5+=stress[3]+stresk[3];
  box->stress6+=stress[4]+stresk[4];
  box->stress7+=stress[2]+stresk[2];
  box->stress8+=stress[4]+stresk[4];
  box->stress9+=stress[5]+stresk[5];
  
// periodic boundary condition
  
//   image_update(atom,simulCond,box);
  image_update(param,box,x,y,z);
  
// updated velocity
  #ifdef _OPENMP
  #pragma omp parallel for default(none) shared(vxu,vyu,vzu,param,atom) private(i)
  #endif
  for(i=parallel->fAtom;i<parallel->lAtom;i++)
  {
    
    vx[i]=vxu[i];
    vy[i]=vyu[i];
    vz[i]=vzu[i];
    
  }
  
  if(parallel->nProc>1)
  {
    update_double_para(param,parallel,x,buffer);
    update_double_para(param,parallel,y,buffer);
    update_double_para(param,parallel,z,buffer);
    
    update_double_para(param,parallel,vx,buffer);
    update_double_para(param,parallel,vy,buffer);
    update_double_para(param,parallel,vz,buffer);
    
    if(parallel->nCtProc>0)
    {
      update_double_para(param,parallel,fx,buffer);
      update_double_para(param,parallel,fy,buffer);
      update_double_para(param,parallel,fz,buffer);
    }
    
  }
}

void lf_nvt_b(PARAM *param,ENERGY *ener,PBC *box,BATH *bath,CONSTRAINT constList[],
	      double *x,double *y,double *z,
	      double *vx,double *vy,double *vz,
	      double *fx,double *fy,double *fz,
	      double *mass,double *rmass,int *nAparallel->tConst)
{
  int i,k,ia,ib,bercycle;
  double lambda,rts2;
  double virshake=0.,virshakt=0.,stress[6]={0.},strest[6]={0.},stresk[6]={0.};
  
  l=0;
  #ifdef _OPENMP
  #pragma omp parallel for default(none) shared(xo,yo,zo,vxo,vyo,vzo,param,atom) private(i)
  #endif
  for(i=parallel->fAtom;i<parallel->lAtom;i++)
  {
    
// Store old coordinates and old velocities.

    xo[l]=x[i];
    yo[l]=y[i];
    zo[l]=z[i];
    
    vxo[l]=vx[i];
    vyo[l]=vy[i];
    vzo[l]=vz[i];
    
    l++;
    
  }
  
  if(param->nConst>0)
  {
//     dd=(DELTA*)my_malloc(param->nConst*sizeof(*dd));
    
    #ifdef _OPENMP
    #pragma omp parallel for default(none) shared(param,atom,constList,dd) private(i,ia,ib)
    #endif
    for(i=0;i<param->nConst;i++)
    {
      ia=constList[i].a;
      ib=constList[i].b;
      
      ddx[i]=x[ib]-x[ia];
      ddy[i]=y[ib]-y[ia];
      ddz[i]=z[ib]-z[ia];
    }
    
//     image_array(param->nConst,dd,simulCond,box);
    image_array(box,ddx,ddy,ddz,param->nConst);
    
//     xt=(double*)my_malloc(param->nAtom*sizeof(*xt));
//     yt=(double*)my_malloc(param->nAtom*sizeof(*yt));
//     zt=(double*)my_malloc(param->nAtom*sizeof(*zt));
    
  }
  
  rts2=1./X2(param->timeStep);
  
  #ifdef _OPENMP
  #pragma omp parallel for default(none) shared(param,atom) private(i)
  #endif
  for(i=parallel->fAtom;i<parallel->lAtom;i++)
  { 
    vx[i]+=0.5*param->timeStep*fx[i]*rmass[i];
    vy[i]+=0.5*param->timeStep*fy[i]*rmass[i];
    vz[i]+=0.5*param->timeStep*fz[i]*rmass[i];
  }
  
  ener->kin=kinetic(param,vx,vy,vz,mass);
  
  if(param->nConst>0)
    bercycle=2;
  else
    bercycle=3;
  
  for(k=0;k<bercycle;k++)
  {
   
    lambda=sqrt(1.0+param->timeStep/bath->tauT*(param->kinTemp0/ener->kin-1.0));
    
// move atoms by leapfrog algorithm
    
    #ifdef _OPENMP
    #pragma omp parallel for default(none) shared(param,atom,xo,yo,zo,xt,yt,zt,vxo,vyo,vzo,vxu,vyu,vzu,lambda) private(i)
    #endif
    for(i=parallel->fAtom;i<parallel->lAtom;i++)
    {
      
// update velocities
      
      vxu[i]=(vxo[i]+param->timeStep*fx[i]*rmass[i])*lambda;
      vyu[i]=(vyo[i]+param->timeStep*fy[i]*rmass[i])*lambda;
      vzu[i]=(vzo[i]+param->timeStep*fz[i]*rmass[i])*lambda;
      
// update positions
      
      x[i]=xo[i]+param->timeStep*vxu[i];
      y[i]=yo[i]+param->timeStep*vyu[i];
      z[i]=zo[i]+param->timeStep*vzu[i];
      
// Temporary storage of the uncorrected positions
      
      if(param->nConst>0)
      {
	xt[i]=x[i];
	yt[i]=y[i];
	zt[i]=z[i];
      }
      
    }
    
    if( (param->nConst>0) && (k==0) )
    {
// Apply constraint with Shake algorithm.
      
      //lf_shake(atom,simulCond,constList,dd,box,&virshakt,strest);
      lf_shake(param,box,constList,x,y,z,ddx,ddy,ddz,rmass,nAparallel->tConst,strest,&virshakt);
      
      virshake+=virshakt;
      for(i=0;i<6;i++)
	stress[i]+=strest[i];
      
      #ifdef _OPENMP
      #pragma omp parallel for default(none) shared(param,atom,xt,yt,zt,vxu,vyu,vzu,rts2) private(i)
      #endif
      for(i=parallel->fAtom;i<parallel->lAtom;i++)
      {
        
// Corrected velocities
      
	vxu[i]+=(x[i]-xt[i])*param->rTimeStep;
	vyu[i]+=(y[i]-yt[i])*param->rTimeStep;
	vzu[i]+=(z[i]-zt[i])*param->rTimeStep;
      
// Corrected Forces
      
	fx[i]+=(x[i]-xt[i])*mass[i]*rts2;
	fy[i]+=(y[i]-yt[i])*mass[i]*rts2;
	fz[i]+=(z[i]-zt[i])*mass[i]*rts2;
      
      }
    }
    
// calculate full timestep velocity
    #ifdef _OPENMP
    #pragma omp parallel for default(none) shared(vxo,vyo,vzo,vxu,vyu,vzu,param,atom) private(i)
    #endif
    for(i=parallel->fAtom;i<parallel->lAtom;i++)
    {
      
      vx[i]=0.5*(vxo[i]+vxu[i]);
      vy[i]=0.5*(vyo[i]+vyu[i]);
      vz[i]=0.5*(vzo[i]+vzu[i]);
      
    }
    
// calculate kinetic energy
    
    ener->kin=kinetic(param,vx,vy,vz,mass);
    
  }
  
  ener->virshake=virshake;
  
//   stress_kinetic(atom,simulCond,stresk);
  stress_kinetic(param,vx,vy,vz,mass,stresk);
  
  box->stress1+=stress[0]+stresk[0];
  box->stress2+=stress[1]+stresk[1];
  box->stress3+=stress[2]+stresk[2];
  box->stress4+=stress[1]+stresk[1];
  box->stress5+=stress[3]+stresk[3];
  box->stress6+=stress[4]+stresk[4];
  box->stress7+=stress[2]+stresk[2];
  box->stress8+=stress[4]+stresk[4];
  box->stress9+=stress[5]+stresk[5];
  
// periodic boundary condition
  
//   image_update(atom,simulCond,box);
  image_update(param,box,x,y,z);
  
// updated velocity
  
  #ifdef _OPENMP
  #pragma omp parallel for default(none) shared(vxu,vyu,vzu,param,atom) private(i)
  #endif
  for(i=parallel->fAtom;i<parallel->lAtom;i++)
  {
    
    vx[i]=vxu[i];
    vy[i]=vyu[i];
    vz[i]=vzu[i];
    
  }
  
// Free temporary arrays
  
//   free(vxu);
//   free(vyu);
//   free(vzu);
//   
//   free(xo);
//   free(yo);
//   free(zo);
//   
//   free(vxo);
//   free(vyo);
//   free(vzo);
//   
//   if(param->nConst>0)
//   {
//     free(dd);
//     
//     free(xt);
//     free(yt);
//     free(zt);
//   }
      
}

void lf_npt_b(PARAM *param,ENERGY *ener,PBC *box,BATH *bath,CONSTRAINT constList[],
	      double *x,double *y,double *z,
	      double *vx,double *vy,double *vz,
	      double *fx,double *fy,double *fz,
	      double *mass,double *rmass,int *nAparallel->tConst)
{
  int i,k,ia,ib,bercycle;
  double lambda,gamma,cbrga,pp,rts2;
//   double *xo=NULL,*yo=NULL,*zo=NULL;
//   double *vxo=NULL,*vyo=NULL,*vzo=NULL;
//   double *xt=NULL,*yt=NULL,*zt=NULL;
//   double *vxu=NULL,*vyu=NULL,*vzu=NULL;
  double volume,cell0[9];
  double virshake=0.,virshakt=0.,stress[6]={0.},strest[6]={0.},stresk[6]={0.};
//   DELTA *dd=NULL;
//   
//   vxu=(double*)my_malloc(param->nAtom*sizeof(*vxu));
//   vyu=(double*)my_malloc(param->nAtom*sizeof(*vyu));
//   vzu=(double*)my_malloc(param->nAtom*sizeof(*vzu));
//   
//   xo=(double*)my_malloc(param->nAtom*sizeof(*xo));
//   yo=(double*)my_malloc(param->nAtom*sizeof(*yo));
//   zo=(double*)my_malloc(param->nAtom*sizeof(*zo));
//   
//   vxo=(double*)my_malloc(param->nAtom*sizeof(*vxo));
//   vyo=(double*)my_malloc(param->nAtom*sizeof(*vyo));
//   vzo=(double*)my_malloc(param->nAtom*sizeof(*vzo));
    
  #ifdef _OPENMP
  #pragma omp parallel for default(none) shared(xo,yo,zo,vxo,vyo,vzo,param,atom) private(i)
  #endif
  for(i=parallel->fAtom;i<parallel->lAtom;i++)
  {
    
// Store old coordinates and old velocities.

    xo[i]=x[i];
    yo[i]=y[i];
    zo[i]=z[i];
    
    vxo[i]=vx[i];
    vyo[i]=vy[i];
    vzo[i]=vz[i];
    
  }
  
  //Store initial box parameters
  
  cell0[0]=box->a1;
  cell0[1]=box->a2;
  cell0[2]=box->a3;
  cell0[3]=box->b1;
  cell0[4]=box->b2;
  cell0[5]=box->b3;
  cell0[6]=box->c1;
  cell0[7]=box->c2;
  cell0[8]=box->c3;
  
  volume=box->vol;
  
  if(param->nConst>0)
  {
//     dd=(DELTA*)my_malloc(param->nConst*sizeof(*dd));
    
    #ifdef _OPENMP
    #pragma omp parallel for default(none) shared(param,atom,constList,dd) private(i,ia,ib)
    #endif
    for(i=0;i<param->nConst;i++)
    {
      ia=constList[i].a;
      ib=constList[i].b;
      
      ddx[i]=x[ib]-x[ia];
      ddy[i]=y[ib]-y[ia];
      ddz[i]=z[ib]-z[ia];
    }
    
//     image_array(param->nConst,dd,simulCond,box);
    image_array(box,ddx,ddy,ddz,param->nConst);
    
//     xt=(double*)my_malloc(param->nAtom*sizeof(*xt));
//     yt=(double*)my_malloc(param->nAtom*sizeof(*yt));
//     zt=(double*)my_malloc(param->nAtom*sizeof(*zt));
    
  }
  
  rts2=1./X2(param->timeStep);
  
  #ifdef _OPENMP
  #pragma omp parallel for default(none) shared(param,atom) private(i)
  #endif
  for(i=parallel->fAtom;i<parallel->lAtom;i++)
  { 
    vx[i]+=0.5*param->timeStep*fx[i]*rmass[i];
    vy[i]+=0.5*param->timeStep*fy[i]*rmass[i];
    vz[i]+=0.5*param->timeStep*fz[i]*rmass[i];
  }
  
  ener->kin=kinetic(param,vx,vy,vz,mass);
  
  pp=(2.*ener->kin-ener->virpot-virshake)/(3.*volume);
  gamma=1.+watercomp*param->timeStep*(pp-param->press0)/bath->tauP;
  cbrga=cbrt(gamma);
  
  lambda=sqrt(1.0+param->timeStep/bath->tauT*(param->kinTemp0/ener->kin-1.0));
    
  if(param->nConst>0)
    bercycle=4;
  else
    bercycle=5;
  
  for(k=0;k<bercycle;k++)
  {
    
// move atoms by leapfrog algorithm
    
    #ifdef _OPENMP
    #pragma omp parallel for default(none) shared(param,atom,xo,yo,zo,xt,yt,zt,vxo,vyo,vzo,vxu,vyu,vzu,lambda,cbrga) private(i)
    #endif
    for(i=parallel->fAtom;i<parallel->lAtom;i++)
    {
      
// update velocities
      
      vxu[i]=(vxo[i]+param->timeStep*fx[i]*rmass[i])*lambda;
      vyu[i]=(vyo[i]+param->timeStep*fy[i]*rmass[i])*lambda;
      vzu[i]=(vzo[i]+param->timeStep*fz[i]*rmass[i])*lambda;
      
// update positions
      
      x[i]=cbrga*xo[i]+param->timeStep*vxu[i];
      y[i]=cbrga*yo[i]+param->timeStep*vyu[i];
      z[i]=cbrga*zo[i]+param->timeStep*vzu[i];
      
// Temporary storage of the uncorrected positions
      
      if(param->nConst>0)
      {
	xt[i]=x[i];
	yt[i]=y[i];
	zt[i]=z[i];
      }
      
    }
    
    if( (param->nConst>0) && (k==0) )
    {
      
//       scale_box(box,cbrga,cell0);
      scale_box(box,cell0,cbrga);
      
// Apply constraint with Shake algorithm.
      
//       lf_shake(atom,simulCond,constList,dd,box,&virshakt,strest);
      lf_shake(param,box,constList,x,y,z,ddx,ddy,ddz,rmass,nAparallel->tConst,strest,&virshakt);
      
      virshake+=virshakt;
      for(i=0;i<6;i++)
	stress[i]+=strest[i];
      
      #ifdef _OPENMP
      #pragma omp parallel for default(none) shared(param,atom,xt,yt,zt,vxu,vyu,vzu,rts2) private(i)
      #endif
      for(i=parallel->fAtom;i<parallel->lAtom;i++)
      {
        
// Corrected velocities
      
	vxu[i]+=(x[i]-xt[i])*param->rTimeStep;
	vyu[i]+=(y[i]-yt[i])*param->rTimeStep;
	vzu[i]+=(z[i]-zt[i])*param->rTimeStep;
      
// Corrected Forces
      
	fx[i]+=(x[i]-xt[i])*mass[i]*rts2;
	fy[i]+=(y[i]-yt[i])*mass[i]*rts2;
	fz[i]+=(z[i]-zt[i])*mass[i]*rts2;
      
      }
    }
    
// calculate full timestep velocity
    #ifdef _OPENMP
    #pragma omp parallel for default(none) shared(vxo,vyo,vzo,vxu,vyu,vzu,param,atom) private(i)
    #endif
    for(i=parallel->fAtom;i<parallel->lAtom;i++)
    {
      
      vx[i]=0.5*(vxo[i]+vxu[i]);
      vy[i]=0.5*(vyo[i]+vyu[i]);
      vz[i]=0.5*(vzo[i]+vzu[i]);
      
    }
    
// calculate kinetic energy
    
    ener->kin=kinetic(param,vx,vy,vz,mass);
    
    pp=(2.*ener->kin-ener->virpot-virshake)/(3.*volume);
    gamma=1.+watercomp*param->timeStep*(pp-param->press0)/bath->tauP;
    cbrga=cbrt(gamma);
    
    lambda=sqrt(1.0+param->timeStep/bath->tauT*(param->kinTemp0/ener->kin-1.0));
    
  }
  
  ener->virshake=virshake;
  
//   stress_kinetic(atom,simulCond,stresk);
  stress_kinetic(param,vx,vy,vz,mass,stresk);
  
  box->stress1+=stress[0]+stresk[0];
  box->stress2+=stress[1]+stresk[1];
  box->stress3+=stress[2]+stresk[2];
  box->stress4+=stress[1]+stresk[1];
  box->stress5+=stress[3]+stresk[3];
  box->stress6+=stress[4]+stresk[4];
  box->stress7+=stress[2]+stresk[2];
  box->stress8+=stress[4]+stresk[4];
  box->stress9+=stress[5]+stresk[5];
  
//   scale_box(box,cbrga,cell0);
  scale_box(box,cell0,cbrga);
  
// periodic boundary condition
  
//   image_update(atom,simulCond,box);
  image_update(param,box,x,y,z);
  
// updated velocity
  
  #ifdef _OPENMP
  #pragma omp parallel for default(none) shared(vxu,vyu,vzu,param,atom) private(i)
  #endif
  for(i=parallel->fAtom;i<parallel->lAtom;i++)
  {
    
    vx[i]=vxu[i];
    vy[i]=vyu[i];
    vz[i]=vzu[i];
    
  }
  
// Free temporary arrays
/*  
  free(vxu);
  free(vyu);
  free(vzu);
  
  free(xo);
  free(yo);
  free(zo);
  
  free(vxo);
  free(vyo);
  free(vzo);
  
  if(param->nConst>0)
  {
    free(dd);
    
    free(xt);
    free(yt);
    free(zt);
  }*/
      
}

void lf_nvt_h(PARAM *param,ENERGY *ener,PBC *box,BATH *bath,CONSTRAINT constList[],
	      double *x,double *y,double *z,
	      double *vx,double *vy,double *vz,
	      double *fx,double *fy,double *fz,
	      double *mass,double *rmass,int *nAparallel->tConst)
{
  int i,k,ia,ib,nosecycle;
  double lambda,lambdb,lambdc,rts2,qmass;
//   double *xo=NULL,*yo=NULL,*zo=NULL;
//   double *vxo=NULL,*vyo=NULL,*vzo=NULL;
//   double *xt=NULL,*yt=NULL,*zt=NULL;
//   double *vxu=NULL,*vyu=NULL,*vzu=NULL;
  double virshake=0.,virshakt=0.,stress[6]={0.},strest[6]={0.},stresk[6]={0.};
//   DELTA *dd=NULL;
//   
//   vxu=(double*)my_malloc(param->nAtom*sizeof(*vxu));
//   vyu=(double*)my_malloc(param->nAtom*sizeof(*vyu));
//   vzu=(double*)my_malloc(param->nAtom*sizeof(*vzu));
//   
//   xo=(double*)my_malloc(param->nAtom*sizeof(*xo));
//   yo=(double*)my_malloc(param->nAtom*sizeof(*yo));
//   zo=(double*)my_malloc(param->nAtom*sizeof(*zo));
//   
//   vxo=(double*)my_malloc(param->nAtom*sizeof(*vxo));
//   vyo=(double*)my_malloc(param->nAtom*sizeof(*vyo));
//   vzo=(double*)my_malloc(param->nAtom*sizeof(*vzo));
    
  
  #ifdef _OPENMP
  #pragma omp parallel for default(none) shared(xo,yo,zo,vxo,vyo,vzo,param,atom) private(i)
  #endif
  for(i=parallel->fAtom;i<parallel->lAtom;i++)
  {
    
// Store old coordinates and old velocities.

    xo[i]=x[i];
    yo[i]=y[i];
    zo[i]=z[i];
    
    vxo[i]=vx[i];
    vyo[i]=vy[i];
    vzo[i]=vz[i];
    
  }
  
  if(param->nConst>0)
  {
//     dd=(DELTA*)my_malloc(param->nConst*sizeof(*dd));
    
    #ifdef _OPENMP
    #pragma omp parallel for default(none) shared(param,atom,dd,constList) private(i,ia,ib)
    #endif
    for(i=0;i<param->nConst;i++)
    {
      ia=constList[i].a;
      ib=constList[i].b;
      
      ddx[i]=x[ib]-x[ia];
      ddy[i]=y[ib]-y[ia];
      ddz[i]=z[ib]-z[ia];
    }
    
//     image_array(param->nConst,dd,simulCond,box);
    image_array(box,ddx,ddy,ddz,param->nConst);
    
//     xt=(double*)my_malloc(param->nAtom*sizeof(*xt));
//     yt=(double*)my_malloc(param->nAtom*sizeof(*yt));
//     zt=(double*)my_malloc(param->nAtom*sizeof(*zt));
    
  }
  
  //   Mass parameter for Nose-Hoover thermostat
  
  qmass=2.0*param->kinTemp0*X2(bath->tauT);
  
  rts2=1./X2(param->timeStep);
  
  #ifdef _OPENMP
  #pragma omp parallel for default(none) shared(param,atom) private(i)
  #endif
  for(i=parallel->fAtom;i<parallel->lAtom;i++)
  { 
    vx[i]+=0.5*param->timeStep*fx[i]*rmass[i];
    vy[i]+=0.5*param->timeStep*fy[i]*rmass[i];
    vz[i]+=0.5*param->timeStep*fz[i]*rmass[i];
  }
  
  ener->kin=kinetic(param,vx,vy,vz,mass);
  
  lambdb=2.0*(ener->kin-param->kinTemp0)/qmass;
  lambdc=bath->chiT+param->timeStep*lambdb;
  lambda=0.5*(bath->chiT+lambdc);
  
  if(param->nConst>0)
    nosecycle=3;
  else
    nosecycle=4;
  
  for(k=0;k<nosecycle;k++)
  {
    
// move atoms by leapfrog algorithm
    
    #ifdef _OPENMP
    #pragma omp parallel for default(none) shared(xt,yt,zt,xo,yo,zo,vxo,vyo,vzo,vxu,vyu,vzu,param,atom,lambda) private(i)
    #endif
    for(i=parallel->fAtom;i<parallel->lAtom;i++)
    {
      
// update velocities
      
      vxu[i]=vxo[i]+param->timeStep*(fx[i]*rmass[i]-vx[i]*lambda);
      vyu[i]=vyo[i]+param->timeStep*(fy[i]*rmass[i]-vy[i]*lambda);
      vzu[i]=vzo[i]+param->timeStep*(fz[i]*rmass[i]-vz[i]*lambda);
      
// update positions
      
      x[i]=xo[i]+param->timeStep*vxu[i];
      y[i]=yo[i]+param->timeStep*vyu[i];
      z[i]=zo[i]+param->timeStep*vzu[i];
      
// Temporary storage of the uncorrected positions
      
      if(param->nConst>0)
      {
	xt[i]=x[i];
	yt[i]=y[i];
	zt[i]=z[i];
      }
      
    }
    
    if( (param->nConst>0) && (k==0) )
    {
// Apply constraint with Shake algorithm.
      
//       lf_shake(atom,simulCond,constList,dd,box,&virshakt,strest);
      lf_shake(param,box,constList,x,y,z,ddx,ddy,ddz,rmass,nAparallel->tConst,strest,&virshakt);
      
      virshake+=virshakt;
      for(i=0;i<6;i++)
	stress[i]+=strest[i];
      
      #ifdef _OPENMP
      #pragma omp parallel for default(none) shared(xt,yt,zt,vxu,vyu,vzu,param,atom,rts2) private(i)
      #endif
      for(i=parallel->fAtom;i<parallel->lAtom;i++)
      {
        
// Corrected velocities
      
	vxu[i]+=(x[i]-xt[i])*param->rTimeStep;
	vyu[i]+=(y[i]-yt[i])*param->rTimeStep;
	vzu[i]+=(z[i]-zt[i])*param->rTimeStep;
      
// Corrected Forces
      
	fx[i]+=(x[i]-xt[i])*mass[i]*rts2;
	fy[i]+=(y[i]-yt[i])*mass[i]*rts2;
	fz[i]+=(z[i]-zt[i])*mass[i]*rts2;
      
      }
    }
    
// calculate full timestep velocity

    #ifdef _OPENMP
    #pragma omp parallel for default(none) shared(vxo,vyo,vzo,vxu,vyu,vzu,param,atom) private(i)
    #endif
    for(i=parallel->fAtom;i<parallel->lAtom;i++)
    {
      
      vx[i]=0.5*(vxo[i]+vxu[i]);
      vy[i]=0.5*(vyo[i]+vyu[i]);
      vz[i]=0.5*(vzo[i]+vzu[i]);
      
    }
    
// calculate kinetic energy
    
    ener->kin=kinetic(param,vx,vy,vz,mass);
    
    lambdb=2.0*(ener->kin-param->kinTemp0)/qmass;
    lambdc=bath->chiT+param->timeStep*lambdb;
    lambda=0.5*(bath->chiT+lambdc);
    
  }
  
  ener->virshake=virshake;
  
//   stress_kinetic(atom,simulCond,stresk);
  stress_kinetic(param,vx,vy,vz,mass,stresk);
  
  box->stress1+=stress[0]+stresk[0];
  box->stress2+=stress[1]+stresk[1];
  box->stress3+=stress[2]+stresk[2];
  box->stress4+=stress[1]+stresk[1];
  box->stress5+=stress[3]+stresk[3];
  box->stress6+=stress[4]+stresk[4];
  box->stress7+=stress[2]+stresk[2];
  box->stress8+=stress[4]+stresk[4];
  box->stress9+=stress[5]+stresk[5];
  
  bath->chiT=lambdc;
  
  ener->conint+=param->timeStep*lambda*qmass/X2(bath->tauT);
  ener->consv=ener->conint+0.5*qmass*X2(lambda);
  
// periodic boundary condition
  
//   image_update(atom,simulCond,box);
  image_update(param,box,x,y,z);
  
// updated velocity
  
  #ifdef _OPENMP
  #pragma omp parallel for default(none) shared(vxu,vyu,vzu,param,atom) private(i)
  #endif
  for(i=parallel->fAtom;i<parallel->lAtom;i++)
  {
    
    vx[i]=vxu[i];
    vy[i]=vyu[i];
    vz[i]=vzu[i];
    
  }
  
// Free temporary arrays
  
//   free(vxu);
//   free(vyu);
//   free(vzu);
//   
//   free(xo);
//   free(yo);
//   free(zo);
//   
//   free(vxo);
//   free(vyo);
//   free(vzo);
//   
//   if(param->nConst>0)
//   {
//     free(dd);
//     
//     free(xt);
//     free(yt);
//     free(zt);
//   }
      
}

void lf_npt_h(PARAM *param,ENERGY *ener,PBC *box,BATH *bath,CONSTRAINT constList[],
	      double *x,double *y,double *z,
	      double *vx,double *vy,double *vz,
	      double *fx,double *fy,double *fz,
	      double *mass,double *rmass,int *nAparallel->tConst)
{
  int i,k,ia,ib,nosecycle;
  double lambda,lambdb,lambdc,rts2,qmass;
  double gamma,gammb,gammc,cbrga,pmass;
//   double *xo=NULL,*yo=NULL,*zo=NULL;
//   double *vxo=NULL,*vyo=NULL,*vzo=NULL;
//   double *xt=NULL,*yt=NULL,*zt=NULL;
//   double *vxu=NULL,*vyu=NULL,*vzu=NULL;
  double volume,masst=0.,cell0[9],com[3]={0.},vom[3]={0.};
  double virshake=0.,virshakt=0.,stress[6]={0.},strest[6]={0.},stresk[6]={0.};
//   DELTA *dd=NULL;
  
//   vxu=(double*)my_malloc(param->nAtom*sizeof(*vxu));
//   vyu=(double*)my_malloc(param->nAtom*sizeof(*vyu));
//   vzu=(double*)my_malloc(param->nAtom*sizeof(*vzu));
//   
//   xo=(double*)my_malloc(param->nAtom*sizeof(*xo));
//   yo=(double*)my_malloc(param->nAtom*sizeof(*yo));
//   zo=(double*)my_malloc(param->nAtom*sizeof(*zo));
//   
//   vxo=(double*)my_malloc(param->nAtom*sizeof(*vxo));
//   vyo=(double*)my_malloc(param->nAtom*sizeof(*vyo));
//   vzo=(double*)my_malloc(param->nAtom*sizeof(*vzo));
    
  
  #ifdef _OPENMP
  #pragma omp parallel for default(none) shared(xo,yo,zo,vxo,vyo,vzo,param,atom) private(i)
  #endif
  for(i=parallel->fAtom;i<parallel->lAtom;i++)
  {
    
// Store old coordinates and old velocities.

    xo[i]=x[i];
    yo[i]=y[i];
    zo[i]=z[i];
    
    vxo[i]=vx[i];
    vyo[i]=vy[i];
    vzo[i]=vz[i];
    
  }
  
  //Store initial box parameters
  
  cell0[0]=box->a1;
  cell0[1]=box->a2;
  cell0[2]=box->a3;
  cell0[3]=box->b1;
  cell0[4]=box->b2;
  cell0[5]=box->b3;
  cell0[6]=box->c1;
  cell0[7]=box->c2;
  cell0[8]=box->c3;
  
  volume=box->vol;
  
  if(param->nConst>0)
  {
//     dd=(DELTA*)my_malloc(param->nConst*sizeof(*dd));
    
    #ifdef _OPENMP
    #pragma omp parallel for default(none) shared(param,atom,dd,constList) private(i,ia,ib)
    #endif
    for(i=0;i<param->nConst;i++)
    {
      ia=constList[i].a;
      ib=constList[i].b;
      
      ddx[i]=x[ib]-x[ia];
      ddy[i]=y[ib]-y[ia];
      ddz[i]=z[ib]-z[ia];
    }
    
//     image_array(param->nConst,dd,simulCond,box);
    image_array(box,ddx,ddy,ddz,param->nConst);
    
//     xt=(double*)my_malloc(param->nAtom*sizeof(*xt));
//     yt=(double*)my_malloc(param->nAtom*sizeof(*yt));
//     zt=(double*)my_malloc(param->nAtom*sizeof(*zt));
    
  }
  
  // total mass and center of mass
  
  masst=0.;
  com[0]=0.;
  com[1]=0.;
  com[2]=0.;
  for(i=parallel->fAtom;i<parallel->lAtom;i++)
  {
    masst+=mass[i];
    com[0]+=mass[i]*x[i];
    com[1]+=mass[i]*y[i];
    com[2]+=mass[i]*z[i];
  }
  com[0]/=masst;
  com[1]/=masst;
  com[2]/=masst;
  
  //   Mass parameter for Nose-Hoover thermostat
  
  qmass=2.0*param->kinTemp0*X2(bath->tauT);
  pmass=2.0*param->kinTemp0*X2(bath->tauP);
  
  rts2=1./X2(param->timeStep);
  
  #ifdef _OPENMP
  #pragma omp parallel for default(none) shared(param,atom) private(i)
  #endif
  for(i=parallel->fAtom;i<parallel->lAtom;i++)
  { 
    vx[i]+=0.5*param->timeStep*fx[i]*rmass[i];
    vy[i]+=0.5*param->timeStep*fy[i]*rmass[i];
    vz[i]+=0.5*param->timeStep*fz[i]*rmass[i];
  }
  
  ener->kin=kinetic(param,vx,vy,vz,mass);
  
  gammb=(2.0*ener->kin - ener->virpot - virshake - 3.0*param->press0*volume)/pmass-
    bath->chiT*bath->chiP;
  gammc=bath->chiP+param->timeStep*gammb;
  gamma=0.5*(bath->chiP+gammc);
  
  lambdb=(2.0*(ener->kin - param->kinTemp0 ) + pmass*X2(bath->chiP)-
    rboltzui*param->temp0)/qmass;
  lambdc=bath->chiT+param->timeStep*lambdb;
  lambda=0.5*(bath->chiT+lambdc);
  
  if(param->nConst>0)
    nosecycle=4;
  else
    nosecycle=5;
  
  for(k=0;k<nosecycle;k++)
  {
    
// move atoms by leapfrog algorithm
    
    #ifdef _OPENMP
    #pragma omp parallel for default(none) shared(xt,yt,zt,xo,yo,zo,vxo,vyo,vzo,vxu,vyu,vzu,param,atom,lambda,gamma,com,gammc) private(i)
    #endif
    for(i=parallel->fAtom;i<parallel->lAtom;i++)
    {
      
// update velocities
      
      vxu[i]=vxo[i]+param->timeStep*(fx[i]*rmass[i]-vx[i]*(lambda+gamma));
      vyu[i]=vyo[i]+param->timeStep*(fy[i]*rmass[i]-vy[i]*(lambda+gamma));
      vzu[i]=vzo[i]+param->timeStep*(fz[i]*rmass[i]-vz[i]*(lambda+gamma));
      
// update positions
      
      x[i]=xo[i]+param->timeStep*(vxu[i]+gammc*(0.5*(x[i]+xo[i])-com[0]));
      y[i]=yo[i]+param->timeStep*(vyu[i]+gammc*(0.5*(y[i]+yo[i])-com[1]));
      z[i]=zo[i]+param->timeStep*(vzu[i]+gammc*(0.5*(z[i]+zo[i])-com[2]));
      
// Temporary storage of the uncorrected positions
      
      if(param->nConst>0)
      {
	xt[i]=x[i];
	yt[i]=y[i];
	zt[i]=z[i];
      }
      
    }
    
    if(param->nConst>0)
    {
// Apply constraint with Shake algorithm.
      
      cbrga=exp(3.*param->timeStep*gammc);
      cbrga=cbrt(cbrga);
      
//       scale_box(box,cbrga,cell0);
      scale_box(box,cell0,cbrga);
      
//       lf_shake(atom,simulCond,constList,dd,box,&virshakt,strest);
      lf_shake(param,box,constList,x,y,z,ddx,ddy,ddz,rmass,nAparallel->tConst,strest,&virshakt);
      
      virshake+=virshakt;
      for(i=0;i<6;i++)
	stress[i]+=strest[i];
      
      #ifdef _OPENMP
      #pragma omp parallel for default(none) shared(xt,yt,zt,vxu,vyu,vzu,param,atom,rts2) private(i)
      #endif
      for(i=parallel->fAtom;i<parallel->lAtom;i++)
      {
        
// Corrected velocities
      
	vxu[i]+=(x[i]-xt[i])*param->rTimeStep;
	vyu[i]+=(y[i]-yt[i])*param->rTimeStep;
	vzu[i]+=(z[i]-zt[i])*param->rTimeStep;
      
// Corrected Forces
      
	fx[i]+=(x[i]-xt[i])*mass[i]*rts2;
	fy[i]+=(y[i]-yt[i])*mass[i]*rts2;
	fz[i]+=(z[i]-zt[i])*mass[i]*rts2;
      
      }
    }
    
// calculate full timestep velocity

    #ifdef _OPENMP
    #pragma omp parallel for default(none) shared(vxo,vyo,vzo,vxu,vyu,vzu,param,atom) private(i)
    #endif
    for(i=parallel->fAtom;i<parallel->lAtom;i++)
    {
      
      vx[i]=0.5*(vxo[i]+vxu[i]);
      vy[i]=0.5*(vyo[i]+vyu[i]);
      vz[i]=0.5*(vzo[i]+vzu[i]);
      
    }
    
// calculate kinetic energy
    
    ener->kin=kinetic(param,vx,vy,vz,mass);
    
    gammb=(2.0*ener->kin-ener->virpot-virshake-3.0*param->press0*volume)/pmass-
      bath->chiT*bath->chiP;
    gammc=bath->chiP+param->timeStep*gammb;
    gamma=0.5*(bath->chiP+gammc);
    
    lambdb=(2.0*(ener->kin-param->kinTemp0)+pmass*X2(bath->chiP)-
      rboltzui*param->temp0)/qmass;
    lambdc=bath->chiT+param->timeStep*lambdb;
    lambda=0.5*(bath->chiT+lambdc);

  }
  
  ener->virshake=virshake;
  
//   stress_kinetic(atom,simulCond,stresk);
  stress_kinetic(param,vx,vy,vz,mass,stresk);
  
  box->stress1+=stress[0]+stresk[0];
  box->stress2+=stress[1]+stresk[1];
  box->stress3+=stress[2]+stresk[2];
  box->stress4+=stress[1]+stresk[1];
  box->stress5+=stress[3]+stresk[3];
  box->stress6+=stress[4]+stresk[4];
  box->stress7+=stress[2]+stresk[2];
  box->stress8+=stress[4]+stresk[4];
  box->stress9+=stress[5]+stresk[5];
  
  cbrga=exp(3.*param->timeStep*gammc);
  cbrga=cbrt(cbrga);
  
//   scale_box(box,cbrga,cell0);
  scale_box(box,cell0,cbrga);
  
  bath->chiT=lambdc;
  bath->chiP=gammc;
  
  ener->conint+=param->timeStep*lambda*(rboltzui*param->temp0+qmass/X2(bath->tauT));
  ener->consv=ener->conint+param->press0*volume+0.5*(qmass*X2(lambda)+pmass*X2(gamma));
  
// periodic boundary condition
  
//   image_update(atom,simulCond,box);
  image_update(param,box,x,y,z);
  
// updated velocity
  
  #ifdef _OPENMP
  #pragma omp parallel for default(none) shared(vxu,vyu,vzu,param,atom) private(i)
  #endif
  for(i=parallel->fAtom;i<parallel->lAtom;i++)
  {
    
    vx[i]=vxu[i];
    vy[i]=vyu[i];
    vz[i]=vzu[i];
    
  }
  
  vom[0]=0.;
  vom[1]=0.;
  vom[2]=0.;
  for(i=parallel->fAtom;i<parallel->lAtom;i++)
  {
    vom[0]+=mass[i]*vx[i];
    vom[1]+=mass[i]*vy[i];
    vom[2]+=mass[i]*vz[i];
  }
  vom[0]/=masst;
  vom[1]/=masst;
  vom[2]/=masst;
  
  for(i=parallel->fAtom;i<parallel->lAtom;i++)
  {
    vx[i]-=vom[0];
    vy[i]-=vom[1];
    vz[i]-=vom[2];
  }
  
// Free temporary arrays
  
//   free(vxu);
//   free(vyu);
//   free(vzu);
//   
//   free(xo);
//   free(yo);
//   free(zo);
//   
//   free(vxo);
//   free(vyo);
//   free(vzo);
//   
//   if(param->nConst>0)
//   {
//     free(dd);
//     
//     free(xt);
//     free(yt);
//     free(zt);
//   }
      
}

void vv_integrate(CTRL *ctrl,PARAM *param,ENERGY *ener,PBC *box,BATH *bath,
		  CONSTRAINT constList[],double *x,double *y,
		  double *z,double *vx,double *vy,double *vz,
		  double *fx,double *fy,double *fz,
		  double *mass,double *rmass,int *nAparallel->tConst,int stage)
{
  
#ifdef TIMER
  update_timer_begin(TIMER_INTEGRATE,__func__);
#endif
  
  switch (ctrl->ens)
  {
    case NVE:
      vv_nve(param,ener,box,constList,x,y,z,vx,vy,vz,fx,fy,fz,mass,rmass,nAparallel->tConst,stage);
      break;
    case NVT_B:
      vv_nvt_b(param,ener,box,bath,constList,x,y,z,vx,vy,vz,fx,fy,fz,mass,rmass,nAparallel->tConst,stage);
      break;
    case NPT_B:
      vv_npt_b(param,ener,box,bath,constList,x,y,z,vx,vy,vz,fx,fy,fz,mass,rmass,nAparallel->tConst,stage);
      break;
    case NVT_H:
      vv_nvt_h(param,ener,box,bath,constList,x,y,z,vx,vy,vz,fx,fy,fz,mass,rmass,nAparallel->tConst,stage);
      break;
    case NPT_H:
      vv_npt_h(param,ener,box,bath,constList,x,y,z,vx,vy,vz,fx,fy,fz,mass,rmass,nAparallel->tConst,stage);
      break;
    default:
      vv_nve(param,ener,box,constList,x,y,z,vx,vy,vz,fx,fy,fz,mass,rmass,nAparallel->tConst,stage);
      break;
  }
  
#ifdef TIMER
  update_timer_end(TIMER_INTEGRATE,__func__);
#endif
  
}

void vv_nve(PARAM *param,ENERGY *ener,PBC *box,CONSTRAINT constList[],
	    double *x,double *y,double *z,double *vx,double *vy,double *vz,
	    double *fx,double *fy,double *fz,
	    double *mass,double *rmass,int *nAparallel->tConst,int stage)
{
  int i,ia,ib;
  double virshake,stress[6]={0.},stresk[6]={0.};
//   DELTA *dd=NULL;
  
  if(param->nConst>0)
  {
//     dd=(DELTA*)my_malloc(param->nConst*sizeof(*dd));
    
    #ifdef _OPENMP
    #pragma omp parallel for default(none) shared(constList,dd,param,atom) private(i,ia,ib)
    #endif
    for(i=0;i<param->nConst;i++)
    {
      ia=constList[i].a;
      ib=constList[i].b;
      
      ddx[i]=x[ib]-x[ia];
      ddy[i]=y[ib]-y[ia];
      ddz[i]=z[ib]-z[ia];
    }
    
//     image_array(param->nConst,dd,simulCond,box);
    image_array(box,ddx,ddy,ddz,param->nConst);
    
  }

// move atoms by leapfrog algorithm
  
  #ifdef _OPENMP
  #pragma omp parallel for default(none) shared(param,atom) private(i)
  #endif
  for(i=parallel->fAtom;i<parallel->lAtom;i++)
  {
// update velocities
    
    vx[i]+=0.5*param->timeStep*fx[i]*rmass[i];
    vy[i]+=0.5*param->timeStep*fy[i]*rmass[i];
    vz[i]+=0.5*param->timeStep*fz[i]*rmass[i];
  }
  
  if(stage==1)
  {
    #ifdef _OPENMP
    #pragma omp parallel for default(none) shared(param,atom) private(i)
    #endif
    for(i=parallel->fAtom;i<parallel->lAtom;i++)
    {
// update positions
      
      x[i]+=param->timeStep*vx[i];
      y[i]+=param->timeStep*vy[i];
      z[i]+=param->timeStep*vz[i];
      
    }
    
    if(param->nConst>0)
    {
      
// Apply constraint with Shake algorithm.

//       vv_shake_r(atom,simulCond,constList,dd,box,&virshake,stress);
      vv_shake_r(param,box,constList,x,y,z,vx,vy,vz,ddx,ddy,ddz,rmass,nAparallel->tConst,stress,&virshake);
      ener->virshake=virshake;
      
    }
    
  }
  else
  {
// calculate kinetic energy

    if(param->nConst>0)
    {
      
// Apply constraint with Shake algorithm.

      vv_shake_v(param,constList,vx,vy,vz,ddx,ddy,ddz,rmass,nAparallel->tConst);
      
    }
  
    ener->kin=kinetic(param,vx,vy,vz,mass);
  
//     stress_kinetic(atom,simulCond,stresk);
    stress_kinetic(param,vx,vy,vz,mass,stresk);
    
    box->stress1+=stress[0]+stresk[0];
    box->stress2+=stress[1]+stresk[1];
    box->stress3+=stress[2]+stresk[2];
    box->stress4+=stress[1]+stresk[1];
    box->stress5+=stress[3]+stresk[3];
    box->stress6+=stress[4]+stresk[4];
    box->stress7+=stress[2]+stresk[2];
    box->stress8+=stress[4]+stresk[4];
    box->stress9+=stress[5]+stresk[5];
    
  }
  
  if(stage==2)
  {
    
// periodic boundary condition
    
//     image_update(atom,simulCond,box);
    image_update(param,box,x,y,z);
  }
  
//   if(param->nConst>0)
//   {
//     free(dd);
//   }
   
}

void vv_nvt_b(PARAM *param,ENERGY *ener,PBC *box,BATH *bath,CONSTRAINT constList[],
	      double *x,double *y,double *z,double *vx,double *vy,double *vz,
	      double *fx,double *fy,double *fz,
	      double *mass,double *rmass,int *nAparallel->tConst,int stage)
{
  int i,ia,ib;
  double lambda;
  double virshake,stress[6]={0.},stresk[6]={0.};
//   DELTA *dd=NULL;
  
  if(param->nConst>0)
  {
//     dd=(DELTA*)my_malloc(param->nConst*sizeof(*dd));
    
    #ifdef _OPENMP
    #pragma omp parallel for default(none) shared(constList,dd,param,atom) private(i,ia,ib)
    #endif
    for(i=0;i<param->nConst;i++)
    {
      ia=constList[i].a;
      ib=constList[i].b;
      
      ddx[i]=x[ib]-x[ia];
      ddy[i]=y[ib]-y[ia];
      ddz[i]=z[ib]-z[ia];
    }
    
//     image_array(param->nConst,dd,simulCond,box);
    image_array(box,ddx,ddy,ddz,param->nConst);
    
  }

// move atoms by leapfrog algorithm
  
  #ifdef _OPENMP
  #pragma omp parallel for default(none) shared(param,atom) private(i)
  #endif
  for(i=parallel->fAtom;i<parallel->lAtom;i++)
  {
// update velocities
    
    vx[i]+=0.5*param->timeStep*fx[i]*rmass[i];
    vy[i]+=0.5*param->timeStep*fy[i]*rmass[i];
    vz[i]+=0.5*param->timeStep*fz[i]*rmass[i];
  }
  
  if(stage==1)
  {
    #ifdef _OPENMP
    #pragma omp parallel for default(none) shared(param,atom) private(i)
    #endif
    for(i=parallel->fAtom;i<parallel->lAtom;i++)
    {
// update positions
      
      x[i]+=param->timeStep*vx[i];
      y[i]+=param->timeStep*vy[i];
      z[i]+=param->timeStep*vz[i];
      
    }
    
    if(param->nConst>0)
    {
      
// Apply constraint with Shake algorithm.

//       vv_shake_r(atom,simulCond,constList,dd,box,&virshake,stress);
      vv_shake_r(param,box,constList,x,y,z,vx,vy,vz,ddx,ddy,ddz,rmass,nAparallel->tConst,stress,&virshake);
      ener->virshake=virshake;
      
    }
    
  }
  else
  {
// calculate kinetic energy

    if(param->nConst>0)
    {
      
// Apply constraint with Shake algorithm.

//       vv_shake_v(atom,simulCond,constList,dd);
      vv_shake_v(param,constList,vx,vy,vz,ddx,ddy,ddz,rmass,nAparallel->tConst);
      
    }
  
    ener->kin=kinetic(param,vx,vy,vz,mass);
    
    lambda=sqrt(1.0+param->timeStep/bath->tauT*(param->kinTemp0/ener->kin-1.0));
    
    #ifdef _OPENMP
    #pragma omp parallel for default(none) shared(param,atom,lambda) private(i)
    #endif
    for(i=parallel->fAtom;i<parallel->lAtom;i++)
    {
      vx[i]*=lambda;
      vy[i]*=lambda;
      vz[i]*=lambda;
    }
    
    ener->kin*=X2(lambda);
  
//     stress_kinetic(atom,simulCond,stresk);
    stress_kinetic(param,vx,vy,vz,mass,stresk);
    
    box->stress1+=stress[0]+stresk[0];
    box->stress2+=stress[1]+stresk[1];
    box->stress3+=stress[2]+stresk[2];
    box->stress4+=stress[1]+stresk[1];
    box->stress5+=stress[3]+stresk[3];
    box->stress6+=stress[4]+stresk[4];
    box->stress7+=stress[2]+stresk[2];
    box->stress8+=stress[4]+stresk[4];
    box->stress9+=stress[5]+stresk[5];
    
  }
  
  if(stage==2)
  {
    
// periodic boundary condition
    
//     image_update(atom,simulCond,box);
    image_update(param,box,x,y,z);
  }
  
//   if(param->nConst>0)
//   {
//     free(dd);
//   }
   
}

void vv_npt_b(PARAM *param,ENERGY *ener,PBC *box,BATH *bath,CONSTRAINT constList[],
	      double *x,double *y,double *z,double *vx,double *vy,double *vz,
	      double *fx,double *fy,double *fz,
	      double *mass,double *rmass,int *nAparallel->tConst,int stage)
{
  int i,ia,ib,k,nosecycle;
  double lambda,gamma,cbrga,volume,pp;
//   double *xo=NULL,*yo=NULL,*zo=NULL;
//   double *vxo=NULL,*vyo=NULL,*vzo=NULL;
  double virshake,stress[6]={0.},stresk[6]={0.};
//   DELTA *dd=NULL;
  
//   xo=(double*)my_malloc(param->nAtom*sizeof(*xo));
//   yo=(double*)my_malloc(param->nAtom*sizeof(*yo));
//   zo=(double*)my_malloc(param->nAtom*sizeof(*zo));
//   
//   vxo=(double*)my_malloc(param->nAtom*sizeof(*vxo));
//   vyo=(double*)my_malloc(param->nAtom*sizeof(*vyo));
//   vzo=(double*)my_malloc(param->nAtom*sizeof(*vzo));
    
  volume=box->vol;
  
  if(param->nConst>0)
  {
//     dd=(DELTA*)my_malloc(param->nConst*sizeof(*dd));
    
    #ifdef _OPENMP
    #pragma omp parallel for default(none) shared(constList,dd,param,atom) private(i,ia,ib)
    #endif
    for(i=0;i<param->nConst;i++)
    {
      ia=constList[i].a;
      ib=constList[i].b;
      
      ddx[i]=x[ib]-x[ia];
      ddy[i]=y[ib]-y[ia];
      ddz[i]=z[ib]-z[ia];
    }
    
//     image_array(param->nConst,dd,simulCond,box);
    image_array(box,ddx,ddy,ddz,param->nConst);
    
  }
  
  if(stage==1)
  {
    
    ener->kin=kinetic(param,vx,vy,vz,mass);
    
    #ifdef _OPENMP
    #pragma omp parallel for default(none) shared(param,atom) private(i)
    #endif
    for(i=parallel->fAtom;i<parallel->lAtom;i++)
    {
//    update velocities
      
      vx[i]+=0.5*param->timeStep*fx[i]*rmass[i];
      vy[i]+=0.5*param->timeStep*fy[i]*rmass[i];
      vz[i]+=0.5*param->timeStep*fz[i]*rmass[i];
    }
    
    if(param->nConst>0)
    {
      #ifdef _OPENMP
      #pragma omp parallel for default(none) shared(xo,yo,zo,vxo,vyo,vzo,param,atom) private(i)
      #endif
      for(i=parallel->fAtom;i<parallel->lAtom;i++)
      {
	
    // Store old coordinates and old velocities.

	xo[i]=x[i];
	yo[i]=y[i];
	zo[i]=z[i];
	
	vxo[i]=vx[i];
	vyo[i]=vy[i];
	vzo[i]=vz[i];
	
      }
    }
    
    if( (stage==1) && (param->nConst>0) )
      nosecycle=2;
    else
      nosecycle=1;
    
    for(k=0;k<nosecycle;k++)
    {
      cbrga=1.;
      
      if(k==nosecycle-1)
      {
	pp=(2.*ener->kin-ener->virpot-virshake)/(3.*volume);
	gamma=1.+watercomp*param->timeStep*(pp-param->press0)/bath->tauP;
	cbrga=cbrt(gamma);
	
	vv_scale_box(box,cbrga);
	
      }
      
      #ifdef _OPENMP
      #pragma omp parallel for default(none) shared(param,atom,cbrga) private(i)
      #endif
      for(i=parallel->fAtom;i<parallel->lAtom;i++)
      {
  // update positions
	
	x[i]=cbrga*x[i]+param->timeStep*vx[i];
	y[i]=cbrga*y[i]+param->timeStep*vy[i];
	z[i]=cbrga*z[i]+param->timeStep*vz[i];
	
      }
      
      if(param->nConst>0)
      {
	
  // Apply constraint with Shake algorithm.

// 	vv_shake_r(atom,simulCond,constList,dd,box,&virshake,stress);
	vv_shake_r(param,box,constList,x,y,z,vx,vy,vz,ddx,ddy,ddz,rmass,nAparallel->tConst,stress,&virshake);
	ener->virshake=virshake;
	
      }
      
      if(k<nosecycle-1)
      {
	#ifdef _OPENMP
	#pragma omp parallel for default(none) shared(xo,yo,zo,vxo,vyo,vzo,param,atom) private(i)
	#endif
	for(i=parallel->fAtom;i<parallel->lAtom;i++)
	{
	  
      // Store old coordinates and old velocities.

	  x[i]=xo[i];
	  y[i]=yo[i];
	  z[i]=zo[i];
	  
	  vx[i]=vxo[i];
	  vy[i]=vyo[i];
	  vz[i]=vzo[i];
	  
	}
      }
    }
  }
  else
  {
    
    #ifdef _OPENMP
    #pragma omp parallel for default(none) shared(param,atom) private(i)
    #endif
    for(i=parallel->fAtom;i<parallel->lAtom;i++)
    {
//    update velocities
      
      vx[i]+=0.5*param->timeStep*fx[i]*rmass[i];
      vy[i]+=0.5*param->timeStep*fy[i]*rmass[i];
      vz[i]+=0.5*param->timeStep*fz[i]*rmass[i];
    }
    
// calculate kinetic energy

    ener->kin=kinetic(param,vx,vy,vz,mass);
    
    lambda=sqrt(1.0+param->timeStep/bath->tauT*(param->kinTemp0/ener->kin-1.0));
    
    for(i=parallel->fAtom;i<parallel->lAtom;i++)
    {
//    update velocities
      
      vx[i]*=lambda;
      vy[i]*=lambda;
      vz[i]*=lambda;
    }

    if(param->nConst>0)
    {
      
// Apply constraint with Shake algorithm.

//       vv_shake_v(atom,simulCond,constList,dd);
      vv_shake_v(param,constList,vx,vy,vz,ddx,ddy,ddz,rmass,nAparallel->tConst);
      
    }
    
    ener->kin=kinetic(param,vx,vy,vz,mass);
  
//     stress_kinetic(atom,simulCond,stresk);
    stress_kinetic(param,vx,vy,vz,mass,stresk);
    
    box->stress1+=stress[0]+stresk[0];
    box->stress2+=stress[1]+stresk[1];
    box->stress3+=stress[2]+stresk[2];
    box->stress4+=stress[1]+stresk[1];
    box->stress5+=stress[3]+stresk[3];
    box->stress6+=stress[4]+stresk[4];
    box->stress7+=stress[2]+stresk[2];
    box->stress8+=stress[4]+stresk[4];
    box->stress9+=stress[5]+stresk[5];
    
  }
  
  if(stage==2)
  {
    
// periodic boundary condition
    
//     image_update(atom,simulCond,box);
    image_update(param,box,x,y,z);
  }
  
//   free(xo);
//   free(yo);
//   free(zo);
//   
//   free(vxo);
//   free(vyo);
//   free(vzo);
//   
//   if(param->nConst>0)
//   {
//     free(dd);
//   }
   
}

void vv_nvt_h(PARAM *param,ENERGY *ener,PBC *box,BATH *bath,CONSTRAINT constList[],
	      double *x,double *y,double *z,double *vx,double *vy,double *vz,
	      double *fx,double *fy,double *fz,
	      double *mass,double *rmass,int *nAparallel->tConst,int stage)
{
  int i,ia,ib;
  double lambda,qmass;
  double virshake,stress[6]={0.},stresk[6]={0.};
//   DELTA *dd=NULL;
  
  if(param->nConst>0)
  {
//     dd=(DELTA*)my_malloc(param->nConst*sizeof(*dd));
    
    #ifdef _OPENMP
    #pragma omp parallel for default(none) shared(constList,dd,param,atom) private(i,ia,ib)
    #endif
    for(i=0;i<param->nConst;i++)
    {
      ia=constList[i].a;
      ib=constList[i].b;
      
      ddx[i]=x[ib]-x[ia];
      ddy[i]=y[ib]-y[ia];
      ddz[i]=z[ib]-z[ia];
    }
    
//     image_array(param->nConst,dd,simulCond,box);
    image_array(box,ddx,ddy,ddz,param->nConst);
    
  }
  
  qmass=2.0*param->kinTemp0*X2(bath->tauT);
  
  if(stage==1)
  {
    
    ener->kin=kinetic(param,vx,vy,vz,mass);
    
    bath->chiT+=0.5*param->timeStep*(ener->kin-param->kinTemp0)/qmass;
    
    lambda=exp(-0.5*param->timeStep*bath->chiT);
    
    #ifdef _OPENMP
    #pragma omp parallel for default(none) shared(param,atom,lambda) private(i)
    #endif
    for(i=parallel->fAtom;i<parallel->lAtom;i++)
    {
  // scale velocities
      
      vx[i]*=lambda;
      vy[i]*=lambda;
      vz[i]*=lambda;
      
  // update velocities
      
      vx[i]+=0.5*param->timeStep*fx[i]*rmass[i];
      vy[i]+=0.5*param->timeStep*fy[i]*rmass[i];
      vz[i]+=0.5*param->timeStep*fz[i]*rmass[i];
    }
    
    ener->kin*=X2(lambda);
    
    ener->conint+=0.5*param->timeStep*bath->chiT*qmass/X2(bath->tauT);
    
    bath->chiT+=0.5*param->timeStep*(ener->kin-param->kinTemp0)/qmass;
    
    #ifdef _OPENMP
    #pragma omp parallel for default(none) shared(param,atom) private(i)
    #endif
    for(i=parallel->fAtom;i<parallel->lAtom;i++)
    {
// update positions
      
      x[i]+=param->timeStep*vx[i];
      y[i]+=param->timeStep*vy[i];
      z[i]+=param->timeStep*vz[i];
      
    }
    
    if(param->nConst>0)
    {
      
// Apply constraint with Shake algorithm.

//       vv_shake_r(atom,simulCond,constList,dd,box,&virshake,stress);
      vv_shake_r(param,box,constList,x,y,z,vx,vy,vz,ddx,ddy,ddz,rmass,nAparallel->tConst,stress,&virshake);
      ener->virshake=virshake;
      
    }
    
  }
  else
  {
// calculate kinetic energy

    #ifdef _OPENMP
    #pragma omp parallel for default(none) shared(param,atom) private(i)
    #endif
    for(i=parallel->fAtom;i<parallel->lAtom;i++)
    {
  // update velocities
      
      vx[i]+=0.5*param->timeStep*fx[i]*rmass[i];
      vy[i]+=0.5*param->timeStep*fy[i]*rmass[i];
      vz[i]+=0.5*param->timeStep*fz[i]*rmass[i];
    }

    if(param->nConst>0)
    {
      
// Apply constraint with Shake algorithm.

//       vv_shake_v(atom,simulCond,constList,dd);
      vv_shake_v(param,constList,vx,vy,vz,ddx,ddy,ddz,rmass,nAparallel->tConst);
      
    }
  
    ener->kin=kinetic(param,vx,vy,vz,mass);
    
    bath->chiT+=0.5*param->timeStep*(ener->kin-param->kinTemp0)/qmass;
    
    lambda=exp(-0.5*param->timeStep*bath->chiT);
    
    #ifdef _OPENMP
    #pragma omp parallel for default(none) shared(param,atom,lambda) private(i)
    #endif
    for(i=parallel->fAtom;i<parallel->lAtom;i++)
    {
      vx[i]*=lambda;
      vy[i]*=lambda;
      vz[i]*=lambda;
    }
    
    ener->kin*=X2(lambda);
    
    ener->conint+=0.5*param->timeStep*bath->chiT*qmass/X2(bath->tauT);
    
    bath->chiT+=0.5*param->timeStep*(ener->kin-param->kinTemp0)/qmass;
    
    ener->consv=ener->conint+0.5*qmass*X2(bath->chiT);
    
//     stress_kinetic(atom,simulCond,stresk);
    stress_kinetic(param,vx,vy,vz,mass,stresk);
    
    box->stress1+=stress[0]+stresk[0];
    box->stress2+=stress[1]+stresk[1];
    box->stress3+=stress[2]+stresk[2];
    box->stress4+=stress[1]+stresk[1];
    box->stress5+=stress[3]+stresk[3];
    box->stress6+=stress[4]+stresk[4];
    box->stress7+=stress[2]+stresk[2];
    box->stress8+=stress[4]+stresk[4];
    box->stress9+=stress[5]+stresk[5];
    
  }
  
  if(stage==2)
  {
    
// periodic boundary condition
    
//     image_update(atom,simulCond,box);
    image_update(param,box,x,y,z);
  }
  
//   if(param->nConst>0)
//   {
//     free(dd);
//   }
   
}

void vv_npt_h(PARAM *param,ENERGY *ener,PBC *box,BATH *bath,CONSTRAINT constList[],
	      double *x,double *y,double *z,double *vx,double *vy,double *vz,
	      double *fx,double *fy,double *fz,
	      double *mass,double *rmass,int *nAparallel->tConst,int stage)
{
  int i,ia,ib,k,kk,nosecycle,hoovercycle=5;
  double hts,chts,cqts;
  double cons0,lambda,lambda0,qmass;
  double gamma,gamma0,pmass,cbrga,scale;
//   double *xo=NULL,*yo=NULL,*zo=NULL;
//   double *vxo=NULL,*vyo=NULL,*vzo=NULL;
  double volume,volume0,cell0[9],masst=0.,com[3]={0.},vom[3]={0.};
  double virshake,stress[6]={0.},stresk[6]={0.};
//   DELTA *dd=NULL;
  
//   xo=(double*)my_malloc(param->nAtom*sizeof(*xo));
//   yo=(double*)my_malloc(param->nAtom*sizeof(*yo));
//   zo=(double*)my_malloc(param->nAtom*sizeof(*zo));
//   
//   vxo=(double*)my_malloc(param->nAtom*sizeof(*vxo));
//   vyo=(double*)my_malloc(param->nAtom*sizeof(*vyo));
//   vzo=(double*)my_malloc(param->nAtom*sizeof(*vzo));
  
  hts=0.5*param->timeStep;
  chts=hts/(double)hoovercycle;
  cqts=0.25*param->timeStep/(double)hoovercycle;
  
  //Store initial box parameters
  
  cell0[0]=box->a1;
  cell0[1]=box->a2;
  cell0[2]=box->a3;
  cell0[3]=box->b1;
  cell0[4]=box->b2;
  cell0[5]=box->b3;
  cell0[6]=box->c1;
  cell0[7]=box->c2;
  cell0[8]=box->c3;
  
  volume=box->vol;
  volume0=box->vol;
  
  masst=0.;
  for(i=parallel->fAtom;i<parallel->lAtom;i++)
    masst+=mass[i];
  
  qmass=2.0*param->kinTemp0*X2(bath->tauT);
  pmass=2.0*param->kinTemp0*X2(bath->tauP);
  
  if(param->nConst>0)
  {
//     dd=(DELTA*)my_malloc(param->nConst*sizeof(*dd));
    
    #ifdef _OPENMP
    #pragma omp parallel for default(none) shared(param,constList,dd,atom) private(i,ia,ib)
    #endif
    for(i=0;i<param->nConst;i++)
    {
      ia=constList[i].a;
      ib=constList[i].b;
      
      ddx[i]=x[ib]-x[ia];
      ddy[i]=y[ib]-y[ia];
      ddz[i]=z[ib]-z[ia];
    }
    
//     image_array(param->nConst,dd,simulCond,box);
    image_array(box,ddx,ddy,ddz,param->nConst);
    
  }
  
  if(stage==1)
  {
    
    lambda0=bath->chiT;
    gamma0=bath->chiP;
    cons0=ener->conint;
    
    #ifdef _OPENMP
    #pragma omp parallel for default(none) shared(param,xo,yo,zo,vxo,vyo,vzo,atom) private(i)
    #endif
    for(i=parallel->fAtom;i<parallel->lAtom;i++)
    {
      
  // Store old coordinates and old velocities.

      xo[i]=x[i];
      yo[i]=y[i];
      zo[i]=z[i];
      
      vxo[i]=vx[i];
      vyo[i]=vy[i];
      vzo[i]=vz[i];
      
    }
    
    if( (stage==1) && (param->nConst>0) )
      nosecycle=2;
    else
      nosecycle=1;
    
    for(k=0;k<nosecycle;k++)
    {
      
      for(kk=0;kk<hoovercycle;kk++)
      {
	
      // apply nvt
	ener->kin=kinetic(param,vx,vy,vz,mass);
	
	bath->chiT+=0.5*cqts*(2.0*(ener->kin-param->kinTemp0)+
	  pmass*X2(bath->chiP)-rboltzui*param->temp0)/qmass;
	
	lambda=exp(-cqts*bath->chiT);
	
	#ifdef _OPENMP
	#pragma omp parallel for default(none) shared(param,atom,lambda) private(i)
	#endif
	for(i=parallel->fAtom;i<parallel->lAtom;i++)
	{
      // scale velocities
	  
	  vx[i]*=lambda;
	  vy[i]*=lambda;
	  vz[i]*=lambda;
	  
	}
	
	ener->kin*=X2(lambda);
	
	ener->conint+=cqts*bath->chiT*(rboltzui*param->temp0+qmass/X2(bath->tauT));
	
	bath->chiT+=0.5*cqts*(2.0*(ener->kin-param->kinTemp0)+
	  pmass*X2(bath->chiP)-rboltzui*param->temp0)/qmass;
	  
      // apply npt
	  
	bath->chiP+=0.5*chts*(((2.0*ener->kin-ener->virpot-virshake)-
	  3.0*param->press0*volume)/pmass-bath->chiP*bath->chiT);
	
	gamma=exp(-chts*bath->chiP);
	
	for(i=parallel->fAtom;i<parallel->lAtom;i++)
	{
      // scale velocities
	  
	  vx[i]*=gamma;
	  vy[i]*=gamma;
	  vz[i]*=gamma;
	  
	}
	
	ener->kin*=X2(gamma);
	
	volume*=exp(3.0*chts*bath->chiP);
	
	bath->chiP+=0.5*chts*(((2.0*ener->kin-ener->virpot-virshake)-
	  3.0*param->press0*volume)/pmass-bath->chiP*bath->chiT);
	
      // apply nvt
	
	ener->kin=kinetic(param,vx,vy,vz,mass);
	
	bath->chiT+=0.5*cqts*(2.0*(ener->kin-param->kinTemp0)+
	  pmass*X2(bath->chiP)-rboltzui*param->temp0)/qmass;
	
	lambda=exp(-cqts*bath->chiT);
	
	#ifdef _OPENMP
	#pragma omp parallel for default(none) shared(param,atom,lambda) private(i)
	#endif
	for(i=parallel->fAtom;i<parallel->lAtom;i++)
	{
      // scale velocities
	  
	  vx[i]*=lambda;
	  vy[i]*=lambda;
	  vz[i]*=lambda;
	  
	}
	
	ener->kin*=X2(lambda);
	
	ener->conint+=cqts*bath->chiT*(rboltzui*param->temp0+qmass/X2(bath->tauT));
	
	bath->chiT+=0.5*cqts*(2.0*(ener->kin-param->kinTemp0)+
	  pmass*X2(bath->chiP)-rboltzui*param->temp0)/qmass;
      }
      
      scale=cbrt(volume/volume0);
//       scale_box(box,scale,cell0);
      scale_box(box,cell0,scale);
      
      #ifdef _OPENMP
      #pragma omp parallel for default(none) shared(param,atom,lambda) private(i)
      #endif
      for(i=parallel->fAtom;i<parallel->lAtom;i++)
      {
    // update velocities
	
	vx[i]+=0.5*param->timeStep*fx[i]*rmass[i];
	vy[i]+=0.5*param->timeStep*fy[i]*rmass[i];
	vz[i]+=0.5*param->timeStep*fz[i]*rmass[i];
      }
      
      com[0]=0.;
      com[1]=0.;
      com[2]=0.;
      for(i=parallel->fAtom;i<parallel->lAtom;i++)
      {
	com[0]+=mass[i]*x[i];
	com[1]+=mass[i]*y[i];
	com[2]+=mass[i]*z[i];
      }
      com[0]/=masst;
      com[1]/=masst;
      com[2]/=masst;
      
      cbrga=exp(param->timeStep*bath->chiP);
      
      #ifdef _OPENMP
      #pragma omp parallel for default(none) shared(param,atom,com,cbrga) private(i)
      #endif
      for(i=parallel->fAtom;i<parallel->lAtom;i++)
      {
  // update positions
	
	x[i]=cbrga*(x[i]-com[0])+param->timeStep*vx[i]+com[0];
	y[i]=cbrga*(y[i]-com[1])+param->timeStep*vy[i]+com[1];
	z[i]=cbrga*(z[i]-com[2])+param->timeStep*vz[i]+com[2];
	
      }
    
      if(param->nConst>0)
      {
	
  // Apply constraint with Shake algorithm.

// 	vv_shake_r(atom,simulCond,constList,dd,box,&virshake,stress);
	vv_shake_r(param,box,constList,x,y,z,vx,vy,vz,ddx,ddy,ddz,rmass,nAparallel->tConst,stress,&virshake);
	ener->virshake=virshake;
	
      }
      
      if(k<nosecycle-1)
      {
	volume=volume0;
	bath->chiT=lambda0;
	bath->chiP=gamma0;
	ener->conint=cons0;
	
	#ifdef _OPENMP
	#pragma omp parallel for default(none) shared(xo,yo,zo,vxo,vyo,vzo,atom,param) private(i)
	#endif
	for(i=parallel->fAtom;i<parallel->lAtom;i++)
	{
	  
      // Store old coordinates and old velocities.

	  xo[i]=x[i];
	  yo[i]=y[i];
	  zo[i]=z[i];
	  
	  vxo[i]=vx[i];
	  vyo[i]=vy[i];
	  vzo[i]=vz[i];
	  
	}
	
      }
    }
  }
  else
  {
// calculate kinetic energy

    #ifdef _OPENMP
    #pragma omp parallel for default(none) shared(param,atom) private(i)
    #endif
    for(i=parallel->fAtom;i<parallel->lAtom;i++)
    {
  // update velocities
      
      vx[i]+=0.5*param->timeStep*fx[i]*rmass[i];
      vy[i]+=0.5*param->timeStep*fy[i]*rmass[i];
      vz[i]+=0.5*param->timeStep*fz[i]*rmass[i];
    }

    if(param->nConst>0)
    {
      
// Apply constraint with Shake algorithm.

//       vv_shake_v(atom,simulCond,constList,dd);
      vv_shake_v(param,constList,vx,vy,vz,ddx,ddy,ddz,rmass,nAparallel->tConst);
      
    }
    
    for(kk=0;kk<hoovercycle;kk++)
    {
      
    // apply nvt
      ener->kin=kinetic(param,vx,vy,vz,mass);
      
      bath->chiT+=0.5*cqts*(2.0*(ener->kin-param->kinTemp0)+
	pmass*X2(bath->chiP)-rboltzui*param->temp0)/qmass;
      
      lambda=exp(-cqts*bath->chiT);
      
      #ifdef _OPENMP
      #pragma omp parallel for default(none) shared(param,atom,lambda) private(i)
      #endif
      for(i=parallel->fAtom;i<parallel->lAtom;i++)
      {
    // scale velocities
	
	vx[i]*=lambda;
	vy[i]*=lambda;
	vz[i]*=lambda;
	
      }
      
      ener->kin*=X2(lambda);
      
      ener->conint+=cqts*bath->chiT*(rboltzui*param->temp0+qmass/X2(bath->tauT));
      
      bath->chiT+=0.5*cqts*(2.0*(ener->kin-param->kinTemp0)+
	pmass*X2(bath->chiP)-rboltzui*param->temp0)/qmass;
	
    // apply npt
	
      bath->chiP+=0.5*chts*(((2.0*ener->kin-ener->virpot-virshake)-
	3.0*param->press0*volume)/pmass-bath->chiP*bath->chiT);
      
      gamma=exp(-chts*bath->chiP);
      
      for(i=parallel->fAtom;i<parallel->lAtom;i++)
      {
    // scale velocities
	
	vx[i]*=gamma;
	vy[i]*=gamma;
	vz[i]*=gamma;
	
      }
      
      ener->kin*=X2(gamma);
      
      volume*=exp(3.0*chts*bath->chiP);
      
      bath->chiP+=0.5*chts*(((2.0*ener->kin-ener->virpot-virshake)-
	3.0*param->press0*volume)/pmass-bath->chiP*bath->chiT);
      
    // apply nvt
      
      ener->kin=kinetic(param,vx,vy,vz,mass);
      
      bath->chiT+=0.5*cqts*(2.0*(ener->kin-param->kinTemp0)+
	pmass*X2(bath->chiP)-rboltzui*param->temp0)/qmass;
      
      lambda=exp(-cqts*bath->chiT);
      
      #ifdef _OPENMP
      #pragma omp parallel for default(none) shared(param,atom,lambda) private(i)
      #endif
      for(i=parallel->fAtom;i<parallel->lAtom;i++)
      {
    // scale velocities
	
	vx[i]*=lambda;
	vy[i]*=lambda;
	vz[i]*=lambda;
	
      }
      
      ener->kin*=X2(lambda);
      
      ener->conint+=cqts*bath->chiT*(rboltzui*param->temp0+qmass/X2(bath->tauT));
      
      bath->chiT+=0.5*cqts*(2.0*(ener->kin-param->kinTemp0)+
	pmass*X2(bath->chiP)-rboltzui*param->temp0)/qmass;
    }
    
    vom[0]=0.;
    vom[1]=0.;
    vom[2]=0.;
    for(i=parallel->fAtom;i<parallel->lAtom;i++)
    {
      vom[0]+=mass[i]*vx[i];
      vom[1]+=mass[i]*vy[i];
      vom[2]+=mass[i]*vz[i];
    }
    vom[0]/=masst;
    vom[1]/=masst;
    vom[2]/=masst;
    
    for(i=parallel->fAtom;i<parallel->lAtom;i++)
    {
      vx[i]-=vom[0];
      vy[i]-=vom[1];
      vz[i]-=vom[2];
    }
    
    scale=cbrt(volume/volume0);
//     scale_box(box,scale,cell0);
    scale_box(box,cell0,scale);
    
    ener->consv=ener->conint+param->press0*volume+0.5*(qmass*X2(bath->chiT)+pmass*X2(bath->chiP));
    
    ener->kin=kinetic(param,vx,vy,vz,mass);
    
//     stress_kinetic(atom,simulCond,stresk);
    stress_kinetic(param,vx,vy,vz,mass,stresk);
    
    box->stress1+=stress[0]+stresk[0];
    box->stress2+=stress[1]+stresk[1];
    box->stress3+=stress[2]+stresk[2];
    box->stress4+=stress[1]+stresk[1];
    box->stress5+=stress[3]+stresk[3];
    box->stress6+=stress[4]+stresk[4];
    box->stress7+=stress[2]+stresk[2];
    box->stress8+=stress[4]+stresk[4];
    box->stress9+=stress[5]+stresk[5];
    
  }
  
  if(stage==2)
  {
    
// periodic boundary condition
    
//     image_update(atom,simulCond,box);
    image_update(param,box,x,y,z);
  }
  
//   free(xo);
//   free(yo);
//   free(zo);
//   
//   free(vxo);
//   free(vyo);
//   free(vzo);
//   
//   if(param->nConst>0)
//   {
//     free(dd);
//   }
   
}