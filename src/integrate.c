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

void lf_integrate(ATOM atom[], ENERGY *ener, SIMULPARAMS *simulCond,CONSTRAINT *constList,PBC *box)
{
  switch (simulCond->ens)
  {
    case 0:
      lf_nve(atom,ener,simulCond,constList,box);
      break;
    case 1:
      lf_nvt_b(atom,ener,simulCond,constList,box);
      break;
    case 2:
      lf_npt_b(atom,ener,simulCond,constList,box);
      break;
    case 3:
      lf_nvt_h(atom,ener,simulCond,constList,box);
      break;
    case 4:
      lf_npt_h(atom,ener,simulCond,constList,box);
      break;
    default:
      lf_nve(atom,ener,simulCond,constList,box);
      break;
  }
}

void lf_nve(ATOM atom[], ENERGY *ener, SIMULPARAMS *simulCond,CONSTRAINT *constList,PBC *box)
{
  int i,ia,ib;
  double *xo=NULL,*yo=NULL,*zo=NULL,*vxu=NULL,*vyu=NULL,*vzu=NULL;
  double virshake=0.,stress[6]={0.},stresk[6]={0.};
  DELTA *dd=NULL;
  
  vxu=(double*)malloc(simulCond->natom*sizeof(*vxu));
  vyu=(double*)malloc(simulCond->natom*sizeof(*vyu));
  vzu=(double*)malloc(simulCond->natom*sizeof(*vzu));
  
  if(simulCond->nconst>0)
  {
    dd=(DELTA*)malloc(simulCond->nconst*sizeof(*dd));
    
    #ifdef _OPENMP
    #pragma omp parallel for default(none) shared(simulCond,atom,constList,dd) private(i,ia,ib)
    #endif
    for(i=0;i<simulCond->nconst;i++)
    {
      ia=constList[i].a;
      ib=constList[i].b;
      
      dd[i].x=atom[ib].x-atom[ia].x;
      dd[i].y=atom[ib].y-atom[ia].y;
      dd[i].z=atom[ib].z-atom[ia].z;
    }
    
    image_array(simulCond->nconst,dd,simulCond,box);
    
    xo=(double*)malloc(simulCond->natom*sizeof(*xo));
    yo=(double*)malloc(simulCond->natom*sizeof(*yo));
    zo=(double*)malloc(simulCond->natom*sizeof(*zo));
    
    #ifdef _OPENMP
    #pragma omp parallel for default(none) shared(simulCond,xo,yo,zo,atom) private(i)
    #endif
    for(i=0;i<simulCond->natom;i++)
    {
      
// Store old coordinates.

      xo[i]=atom[i].x;
      yo[i]=atom[i].y;
      zo[i]=atom[i].z;
      
    }
    
  }

// move atoms by leapfrog algorithm
  #ifdef _OPENMP
  #pragma omp parallel for default(none) shared(simulCond,vxu,vyu,vzu,atom) private(i)
  #endif
  for(i=0;i<simulCond->natom;i++)
  {
    
// update velocities
    
    vxu[i]=atom[i].vx+simulCond->timeStep*atom[i].fx/atom[i].m;
    vyu[i]=atom[i].vy+simulCond->timeStep*atom[i].fy/atom[i].m;
    vzu[i]=atom[i].vz+simulCond->timeStep*atom[i].fz/atom[i].m;
    
// update positions
    
    atom[i].x+=simulCond->timeStep*vxu[i];
    atom[i].y+=simulCond->timeStep*vyu[i];
    atom[i].z+=simulCond->timeStep*vzu[i];
    
  }
  
  if(simulCond->nconst>0)
  {
// Apply constraint with Shake algorithm.

    lf_shake(atom,simulCond,constList,dd,box,&virshake,stress);
    
    #ifdef _OPENMP
    #pragma omp parallel for default(none) shared(simulCond,vxu,vyu,vzu,xo,yo,zo,atom) private(i)
    #endif
    for(i=0;i<simulCond->natom;i++)
    {
        
// Corrected velocities
    
      vxu[i]=(atom[i].x-xo[i])/simulCond->timeStep;
      vyu[i]=(atom[i].y-yo[i])/simulCond->timeStep;
      vzu[i]=(atom[i].z-zo[i])/simulCond->timeStep;
    
// Corrected Forces
    
      atom[i].fx=(vxu[i]-atom[i].vx)*atom[i].m/simulCond->timeStep;
      atom[i].fy=(vyu[i]-atom[i].vy)*atom[i].m/simulCond->timeStep;
      atom[i].fz=(vzu[i]-atom[i].vz)*atom[i].m/simulCond->timeStep;
    
    }
  }
  
// calculate full timestep velocity
  #ifdef _OPENMP
  #pragma omp parallel for default(none) shared(simulCond,vxu,vyu,vzu,atom) private(i)
  #endif
  for(i=0;i<simulCond->natom;i++)
  {
    
    atom[i].vx=0.5*(atom[i].vx+vxu[i]);
    atom[i].vy=0.5*(atom[i].vy+vyu[i]);
    atom[i].vz=0.5*(atom[i].vz+vzu[i]);
    
  }
  
// calculate kinetic energy
  
  ener->kin=kinetic(atom,simulCond);
  
  ener->virshake=virshake;
  
  stress_kinetic(atom,simulCond,stresk);
  
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
  
  image_update(atom,simulCond,box);
  
// updated velocity
  #ifdef _OPENMP
  #pragma omp parallel for default(none) shared(simulCond,vxu,vyu,vzu,atom) private(i)
  #endif
  for(i=0;i<simulCond->natom;i++)
  {
    
    atom[i].vx=vxu[i];
    atom[i].vy=vyu[i];
    atom[i].vz=vzu[i];
    
  }
  
  free(vxu);
  free(vyu);
  free(vzu);
  
  if(simulCond->nconst>0)
  {
    free(xo);
    free(yo);
    free(zo);
    free(dd);
  }
      
}

void lf_nvt_b(ATOM atom[], ENERGY *ener, SIMULPARAMS *simulCond,CONSTRAINT *constList,PBC *box)
{
  int i,k,ia,ib,bercycle;
  double lambda,ts2;
  double *xo=NULL,*yo=NULL,*zo=NULL;
  double *vxo=NULL,*vyo=NULL,*vzo=NULL;
  double *xt=NULL,*yt=NULL,*zt=NULL;
  double *vxu=NULL,*vyu=NULL,*vzu=NULL;
  double virshake=0.,virshakt=0.,stress[6]={0.},strest[6]={0.},stresk[6]={0.};
  DELTA *dd=NULL;
  
  vxu=(double*)malloc(simulCond->natom*sizeof(*vxu));
  vyu=(double*)malloc(simulCond->natom*sizeof(*vyu));
  vzu=(double*)malloc(simulCond->natom*sizeof(*vzu));
  
  xo=(double*)malloc(simulCond->natom*sizeof(*xo));
  yo=(double*)malloc(simulCond->natom*sizeof(*yo));
  zo=(double*)malloc(simulCond->natom*sizeof(*zo));
  
  vxo=(double*)malloc(simulCond->natom*sizeof(*vxo));
  vyo=(double*)malloc(simulCond->natom*sizeof(*vyo));
  vzo=(double*)malloc(simulCond->natom*sizeof(*vzo));
    
  #ifdef _OPENMP
  #pragma omp parallel for default(none) shared(simulCond,xo,yo,zo,vxo,vyo,vzo,atom) private(i)
  #endif
  for(i=0;i<simulCond->natom;i++)
  {
    
// Store old coordinates and old velocities.

    xo[i]=atom[i].x;
    yo[i]=atom[i].y;
    zo[i]=atom[i].z;
    
    vxo[i]=atom[i].vx;
    vyo[i]=atom[i].vy;
    vzo[i]=atom[i].vz;
    
  }
  
  if(simulCond->nconst>0)
  {
    dd=(DELTA*)malloc(simulCond->nconst*sizeof(*dd));
    
    #ifdef _OPENMP
    #pragma omp parallel for default(none) shared(simulCond,atom,constList,dd) private(i,ia,ib)
    #endif
    for(i=0;i<simulCond->nconst;i++)
    {
      ia=constList[i].a;
      ib=constList[i].b;
      
      dd[i].x=atom[ib].x-atom[ia].x;
      dd[i].y=atom[ib].y-atom[ia].y;
      dd[i].z=atom[ib].z-atom[ia].z;
    }
    
    image_array(simulCond->nconst,dd,simulCond,box);
    
    xt=(double*)malloc(simulCond->natom*sizeof(*xt));
    yt=(double*)malloc(simulCond->natom*sizeof(*yt));
    zt=(double*)malloc(simulCond->natom*sizeof(*zt));
    
  }
  
  ts2=X2(simulCond->timeStep);
  
  #ifdef _OPENMP
  #pragma omp parallel for default(none) shared(simulCond,atom) private(i)
  #endif
  for(i=0;i<simulCond->natom;i++)
  { 
    atom[i].vx+=0.5*simulCond->timeStep*atom[i].fx/atom[i].m;
    atom[i].vy+=0.5*simulCond->timeStep*atom[i].fy/atom[i].m;
    atom[i].vz+=0.5*simulCond->timeStep*atom[i].fz/atom[i].m;
  }
  
  ener->kin=kinetic(atom,simulCond);
  
  if(simulCond->nconst>0)
    bercycle=2;
  else
    bercycle=3;
  
  for(k=0;k<bercycle;k++)
  {
   
    lambda=sqrt(1.0+simulCond->timeStep/simulCond->taut*(simulCond->kintemp0/ener->kin-1.0));
    
// move atoms by leapfrog algorithm
    
    #ifdef _OPENMP
    #pragma omp parallel for default(none) shared(simulCond,atom,xo,yo,zo,xt,yt,zt,vxo,vyo,vzo,vxu,vyu,vzu,lambda) private(i)
    #endif
    for(i=0;i<simulCond->natom;i++)
    {
      
// update velocities
      
      vxu[i]=(vxo[i]+simulCond->timeStep*atom[i].fx/atom[i].m)*lambda;
      vyu[i]=(vyo[i]+simulCond->timeStep*atom[i].fy/atom[i].m)*lambda;
      vzu[i]=(vzo[i]+simulCond->timeStep*atom[i].fz/atom[i].m)*lambda;
      
// update positions
      
      atom[i].x=xo[i]+simulCond->timeStep*vxu[i];
      atom[i].y=yo[i]+simulCond->timeStep*vyu[i];
      atom[i].z=zo[i]+simulCond->timeStep*vzu[i];
      
// Temporary storage of the uncorrected positions
      
      if(simulCond->nconst>0)
      {
	xt[i]=atom[i].x;
	yt[i]=atom[i].y;
	zt[i]=atom[i].z;
      }
      
    }
    
    if( (simulCond->nconst>0) && (k==0) )
    {
// Apply constraint with Shake algorithm.
      
      lf_shake(atom,simulCond,constList,dd,box,&virshakt,strest);
      
      virshake+=virshakt;
      for(i=0;i<6;i++)
	stress[i]+=strest[i];
      
      #ifdef _OPENMP
      #pragma omp parallel for default(none) shared(simulCond,atom,xt,yt,zt,vxu,vyu,vzu,ts2) private(i)
      #endif
      for(i=0;i<simulCond->natom;i++)
      {
        
// Corrected velocities
      
	vxu[i]+=(atom[i].x-xt[i])/simulCond->timeStep;
	vyu[i]+=(atom[i].y-yt[i])/simulCond->timeStep;
	vzu[i]+=(atom[i].z-zt[i])/simulCond->timeStep;
      
// Corrected Forces
      
	atom[i].fx+=(atom[i].x-xt[i])*atom[i].m/ts2;
	atom[i].fy+=(atom[i].y-yt[i])*atom[i].m/ts2;
	atom[i].fz+=(atom[i].z-zt[i])*atom[i].m/ts2;
      
      }
    }
    
// calculate full timestep velocity
    #ifdef _OPENMP
    #pragma omp parallel for default(none) shared(simulCond,vxo,vyo,vzo,vxu,vyu,vzu,atom) private(i)
    #endif
    for(i=0;i<simulCond->natom;i++)
    {
      
      atom[i].vx=0.5*(vxo[i]+vxu[i]);
      atom[i].vy=0.5*(vyo[i]+vyu[i]);
      atom[i].vz=0.5*(vzo[i]+vzu[i]);
      
    }
    
// calculate kinetic energy
    
    ener->kin=kinetic(atom,simulCond);
    
  }
  
  ener->virshake=virshake;
  
  stress_kinetic(atom,simulCond,stresk);
  
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
  
  image_update(atom,simulCond,box);
  
// updated velocity
  
  #ifdef _OPENMP
  #pragma omp parallel for default(none) shared(simulCond,vxu,vyu,vzu,atom) private(i)
  #endif
  for(i=0;i<simulCond->natom;i++)
  {
    
    atom[i].vx=vxu[i];
    atom[i].vy=vyu[i];
    atom[i].vz=vzu[i];
    
  }
  
// Free temporary arrays
  
  free(vxu);
  free(vyu);
  free(vzu);
  
  free(xo);
  free(yo);
  free(zo);
  
  free(vxo);
  free(vyo);
  free(vzo);
  
  if(simulCond->nconst>0)
  {
    free(dd);
    
    free(xt);
    free(yt);
    free(zt);
  }
      
}

void lf_npt_b(ATOM atom[], ENERGY *ener, SIMULPARAMS *simulCond,CONSTRAINT *constList,PBC *box)
{
  int i,k,ia,ib,bercycle;
  double lambda,gamma,cbrga,pp,ts2;
  double *xo=NULL,*yo=NULL,*zo=NULL;
  double *vxo=NULL,*vyo=NULL,*vzo=NULL;
  double *xt=NULL,*yt=NULL,*zt=NULL;
  double *vxu=NULL,*vyu=NULL,*vzu=NULL;
  double volume,cell0[9];
  double virshake=0.,virshakt=0.,stress[6]={0.},strest[6]={0.},stresk[6]={0.};
  DELTA *dd=NULL;
  
  vxu=(double*)malloc(simulCond->natom*sizeof(*vxu));
  vyu=(double*)malloc(simulCond->natom*sizeof(*vyu));
  vzu=(double*)malloc(simulCond->natom*sizeof(*vzu));
  
  xo=(double*)malloc(simulCond->natom*sizeof(*xo));
  yo=(double*)malloc(simulCond->natom*sizeof(*yo));
  zo=(double*)malloc(simulCond->natom*sizeof(*zo));
  
  vxo=(double*)malloc(simulCond->natom*sizeof(*vxo));
  vyo=(double*)malloc(simulCond->natom*sizeof(*vyo));
  vzo=(double*)malloc(simulCond->natom*sizeof(*vzo));
    
  #ifdef _OPENMP
  #pragma omp parallel for default(none) shared(simulCond,xo,yo,zo,vxo,vyo,vzo,atom) private(i)
  #endif
  for(i=0;i<simulCond->natom;i++)
  {
    
// Store old coordinates and old velocities.

    xo[i]=atom[i].x;
    yo[i]=atom[i].y;
    zo[i]=atom[i].z;
    
    vxo[i]=atom[i].vx;
    vyo[i]=atom[i].vy;
    vzo[i]=atom[i].vz;
    
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
  
  if(simulCond->nconst>0)
  {
    dd=(DELTA*)malloc(simulCond->nconst*sizeof(*dd));
    
    #ifdef _OPENMP
    #pragma omp parallel for default(none) shared(simulCond,atom,constList,dd) private(i,ia,ib)
    #endif
    for(i=0;i<simulCond->nconst;i++)
    {
      ia=constList[i].a;
      ib=constList[i].b;
      
      dd[i].x=atom[ib].x-atom[ia].x;
      dd[i].y=atom[ib].y-atom[ia].y;
      dd[i].z=atom[ib].z-atom[ia].z;
    }
    
    image_array(simulCond->nconst,dd,simulCond,box);
    
    xt=(double*)malloc(simulCond->natom*sizeof(*xt));
    yt=(double*)malloc(simulCond->natom*sizeof(*yt));
    zt=(double*)malloc(simulCond->natom*sizeof(*zt));
    
  }
  
  ts2=X2(simulCond->timeStep);
  
  #ifdef _OPENMP
  #pragma omp parallel for default(none) shared(simulCond,atom) private(i)
  #endif
  for(i=0;i<simulCond->natom;i++)
  { 
    atom[i].vx+=0.5*simulCond->timeStep*atom[i].fx/atom[i].m;
    atom[i].vy+=0.5*simulCond->timeStep*atom[i].fy/atom[i].m;
    atom[i].vz+=0.5*simulCond->timeStep*atom[i].fz/atom[i].m;
  }
  
  ener->kin=kinetic(atom,simulCond);
  
  pp=(2.*ener->kin-ener->virpot-virshake)/(3.*volume);
  gamma=1.+watercomp*simulCond->timeStep*(pp-simulCond->press)/simulCond->taup;
  cbrga=cbrt(gamma);
  
  lambda=sqrt(1.0+simulCond->timeStep/simulCond->taut*(simulCond->kintemp0/ener->kin-1.0));
    
  if(simulCond->nconst>0)
    bercycle=4;
  else
    bercycle=5;
  
  for(k=0;k<bercycle;k++)
  {
    
// move atoms by leapfrog algorithm
    
    #ifdef _OPENMP
    #pragma omp parallel for default(none) shared(simulCond,atom,xo,yo,zo,xt,yt,zt,vxo,vyo,vzo,vxu,vyu,vzu,lambda) private(i)
    #endif
    for(i=0;i<simulCond->natom;i++)
    {
      
// update velocities
      
      vxu[i]=(vxo[i]+simulCond->timeStep*atom[i].fx/atom[i].m)*lambda;
      vyu[i]=(vyo[i]+simulCond->timeStep*atom[i].fy/atom[i].m)*lambda;
      vzu[i]=(vzo[i]+simulCond->timeStep*atom[i].fz/atom[i].m)*lambda;
      
// update positions
      
      atom[i].x=cbrga*xo[i]+simulCond->timeStep*vxu[i];
      atom[i].y=cbrga*yo[i]+simulCond->timeStep*vyu[i];
      atom[i].z=cbrga*zo[i]+simulCond->timeStep*vzu[i];
      
// Temporary storage of the uncorrected positions
      
      if(simulCond->nconst>0)
      {
	xt[i]=atom[i].x;
	yt[i]=atom[i].y;
	zt[i]=atom[i].z;
      }
      
    }
    
    if( (simulCond->nconst>0) && (k==0) )
    {
      
      scale_box(box,cbrga,cell0);
      
// Apply constraint with Shake algorithm.
      
      lf_shake(atom,simulCond,constList,dd,box,&virshakt,strest);
      
      virshake+=virshakt;
      for(i=0;i<6;i++)
	stress[i]+=strest[i];
      
      #ifdef _OPENMP
      #pragma omp parallel for default(none) shared(simulCond,atom,xt,yt,zt,vxu,vyu,vzu,ts2) private(i)
      #endif
      for(i=0;i<simulCond->natom;i++)
      {
        
// Corrected velocities
      
	vxu[i]+=(atom[i].x-xt[i])/simulCond->timeStep;
	vyu[i]+=(atom[i].y-yt[i])/simulCond->timeStep;
	vzu[i]+=(atom[i].z-zt[i])/simulCond->timeStep;
      
// Corrected Forces
      
	atom[i].fx+=(atom[i].x-xt[i])*atom[i].m/ts2;
	atom[i].fy+=(atom[i].y-yt[i])*atom[i].m/ts2;
	atom[i].fz+=(atom[i].z-zt[i])*atom[i].m/ts2;
      
      }
    }
    
// calculate full timestep velocity
    #ifdef _OPENMP
    #pragma omp parallel for default(none) shared(simulCond,vxo,vyo,vzo,vxu,vyu,vzu,atom) private(i)
    #endif
    for(i=0;i<simulCond->natom;i++)
    {
      
      atom[i].vx=0.5*(vxo[i]+vxu[i]);
      atom[i].vy=0.5*(vyo[i]+vyu[i]);
      atom[i].vz=0.5*(vzo[i]+vzu[i]);
      
    }
    
// calculate kinetic energy
    
    ener->kin=kinetic(atom,simulCond);
    
    pp=(2.*ener->kin-ener->virpot-virshake)/(3.*volume);
    gamma=1.+watercomp*simulCond->timeStep*(pp-simulCond->press)/simulCond->taup;
    cbrga=cbrt(gamma);
    
    lambda=sqrt(1.0+simulCond->timeStep/simulCond->taut*(simulCond->kintemp0/ener->kin-1.0));
    
  }
  
  ener->virshake=virshake;
  
  stress_kinetic(atom,simulCond,stresk);
  
  box->stress1+=stress[0]+stresk[0];
  box->stress2+=stress[1]+stresk[1];
  box->stress3+=stress[2]+stresk[2];
  box->stress4+=stress[1]+stresk[1];
  box->stress5+=stress[3]+stresk[3];
  box->stress6+=stress[4]+stresk[4];
  box->stress7+=stress[2]+stresk[2];
  box->stress8+=stress[4]+stresk[4];
  box->stress9+=stress[5]+stresk[5];
  
  scale_box(box,cbrga,cell0);
  
// periodic boundary condition
  
  image_update(atom,simulCond,box);
  
// updated velocity
  
  #ifdef _OPENMP
  #pragma omp parallel for default(none) shared(simulCond,vxu,vyu,vzu,atom) private(i)
  #endif
  for(i=0;i<simulCond->natom;i++)
  {
    
    atom[i].vx=vxu[i];
    atom[i].vy=vyu[i];
    atom[i].vz=vzu[i];
    
  }
  
// Free temporary arrays
  
  free(vxu);
  free(vyu);
  free(vzu);
  
  free(xo);
  free(yo);
  free(zo);
  
  free(vxo);
  free(vyo);
  free(vzo);
  
  if(simulCond->nconst>0)
  {
    free(dd);
    
    free(xt);
    free(yt);
    free(zt);
  }
      
}

void lf_nvt_h(ATOM atom[], ENERGY *ener, SIMULPARAMS *simulCond,CONSTRAINT *constList,PBC *box)
{
  int i,k,ia,ib,nosecycle;
  double lambda,lambdb,lambdc,ts2,qmass;
  double *xo=NULL,*yo=NULL,*zo=NULL;
  double *vxo=NULL,*vyo=NULL,*vzo=NULL;
  double *xt=NULL,*yt=NULL,*zt=NULL;
  double *vxu=NULL,*vyu=NULL,*vzu=NULL;
  double virshake=0.,virshakt=0.,stress[6]={0.},strest[6]={0.},stresk[6]={0.};
  DELTA *dd=NULL;
  
  vxu=(double*)malloc(simulCond->natom*sizeof(*vxu));
  vyu=(double*)malloc(simulCond->natom*sizeof(*vyu));
  vzu=(double*)malloc(simulCond->natom*sizeof(*vzu));
  
  xo=(double*)malloc(simulCond->natom*sizeof(*xo));
  yo=(double*)malloc(simulCond->natom*sizeof(*yo));
  zo=(double*)malloc(simulCond->natom*sizeof(*zo));
  
  vxo=(double*)malloc(simulCond->natom*sizeof(*vxo));
  vyo=(double*)malloc(simulCond->natom*sizeof(*vyo));
  vzo=(double*)malloc(simulCond->natom*sizeof(*vzo));
    
  
  #ifdef _OPENMP
  #pragma omp parallel for default(none) shared(simulCond,xo,yo,zo,vxo,vyo,vzo,atom) private(i)
  #endif
  for(i=0;i<simulCond->natom;i++)
  {
    
// Store old coordinates and old velocities.

    xo[i]=atom[i].x;
    yo[i]=atom[i].y;
    zo[i]=atom[i].z;
    
    vxo[i]=atom[i].vx;
    vyo[i]=atom[i].vy;
    vzo[i]=atom[i].vz;
    
  }
  
  if(simulCond->nconst>0)
  {
    dd=(DELTA*)malloc(simulCond->nconst*sizeof(*dd));
    
    #ifdef _OPENMP
    #pragma omp parallel for default(none) shared(simulCond,atom,dd,constList) private(i,ia,ib)
    #endif
    for(i=0;i<simulCond->nconst;i++)
    {
      ia=constList[i].a;
      ib=constList[i].b;
      
      dd[i].x=atom[ib].x-atom[ia].x;
      dd[i].y=atom[ib].y-atom[ia].y;
      dd[i].z=atom[ib].z-atom[ia].z;
    }
    
    image_array(simulCond->nconst,dd,simulCond,box);
    
    xt=(double*)malloc(simulCond->natom*sizeof(*xt));
    yt=(double*)malloc(simulCond->natom*sizeof(*yt));
    zt=(double*)malloc(simulCond->natom*sizeof(*zt));
    
  }
  
  //   Mass parameter for Nose-Hoover thermostat
  
  qmass=2.0*simulCond->kintemp0*X2(simulCond->taut);
  
  ts2=X2(simulCond->timeStep);
  
  #ifdef _OPENMP
  #pragma omp parallel for default(none) shared(simulCond,atom) private(i)
  #endif
  for(i=0;i<simulCond->natom;i++)
  { 
    atom[i].vx+=0.5*simulCond->timeStep*atom[i].fx/atom[i].m;
    atom[i].vy+=0.5*simulCond->timeStep*atom[i].fy/atom[i].m;
    atom[i].vz+=0.5*simulCond->timeStep*atom[i].fz/atom[i].m;
  }
  
  ener->kin=kinetic(atom,simulCond);
  
  lambdb=2.0*(ener->kin-simulCond->kintemp0)/qmass;
  lambdc=simulCond->lambdat+simulCond->timeStep*lambdb;
  lambda=0.5*(simulCond->lambdat+lambdc);
  
  if(simulCond->nconst>0)
    nosecycle=3;
  else
    nosecycle=4;
  
  for(k=0;k<nosecycle;k++)
  {
    
// move atoms by leapfrog algorithm
    
    #ifdef _OPENMP
    #pragma omp parallel for default(none) shared(simulCond,xt,yt,zt,xo,yo,zo,vxo,vyo,vzo,vxu,vyu,vzu,atom,lambda) private(i)
    #endif
    for(i=0;i<simulCond->natom;i++)
    {
      
// update velocities
      
      vxu[i]=vxo[i]+simulCond->timeStep*(atom[i].fx/atom[i].m-atom[i].vx*lambda);
      vyu[i]=vyo[i]+simulCond->timeStep*(atom[i].fy/atom[i].m-atom[i].vy*lambda);
      vzu[i]=vzo[i]+simulCond->timeStep*(atom[i].fz/atom[i].m-atom[i].vz*lambda);
      
// update positions
      
      atom[i].x=xo[i]+simulCond->timeStep*vxu[i];
      atom[i].y=yo[i]+simulCond->timeStep*vyu[i];
      atom[i].z=zo[i]+simulCond->timeStep*vzu[i];
      
// Temporary storage of the uncorrected positions
      
      if(simulCond->nconst>0)
      {
	xt[i]=atom[i].x;
	yt[i]=atom[i].y;
	zt[i]=atom[i].z;
      }
      
    }
    
    if( (simulCond->nconst>0) && (k==0) )
    {
// Apply constraint with Shake algorithm.
      
      lf_shake(atom,simulCond,constList,dd,box,&virshakt,strest);
      
      virshake+=virshakt;
      for(i=0;i<6;i++)
	stress[i]+=strest[i];
      
      #ifdef _OPENMP
      #pragma omp parallel for default(none) shared(simulCond,xt,yt,zt,vxu,vyu,vzu,atom,ts2) private(i)
      #endif
      for(i=0;i<simulCond->natom;i++)
      {
        
// Corrected velocities
      
	vxu[i]+=(atom[i].x-xt[i])/simulCond->timeStep;
	vyu[i]+=(atom[i].y-yt[i])/simulCond->timeStep;
	vzu[i]+=(atom[i].z-zt[i])/simulCond->timeStep;
      
// Corrected Forces
      
	atom[i].fx+=(atom[i].x-xt[i])*atom[i].m/ts2;
	atom[i].fy+=(atom[i].y-yt[i])*atom[i].m/ts2;
	atom[i].fz+=(atom[i].z-zt[i])*atom[i].m/ts2;
      
      }
    }
    
// calculate full timestep velocity

    #ifdef _OPENMP
    #pragma omp parallel for default(none) shared(simulCond,vxo,vyo,vzo,vxu,vyu,vzu,atom) private(i)
    #endif
    for(i=0;i<simulCond->natom;i++)
    {
      
      atom[i].vx=0.5*(vxo[i]+vxu[i]);
      atom[i].vy=0.5*(vyo[i]+vyu[i]);
      atom[i].vz=0.5*(vzo[i]+vzu[i]);
      
    }
    
// calculate kinetic energy
    
    ener->kin=kinetic(atom,simulCond);
    
    lambdb=2.0*(ener->kin-simulCond->kintemp0)/qmass;
    lambdc=simulCond->lambdat+simulCond->timeStep*lambdb;
    lambda=0.5*(simulCond->lambdat+lambdc);
    
  }
  
  ener->virshake=virshake;
  
  stress_kinetic(atom,simulCond,stresk);
  
  box->stress1+=stress[0]+stresk[0];
  box->stress2+=stress[1]+stresk[1];
  box->stress3+=stress[2]+stresk[2];
  box->stress4+=stress[1]+stresk[1];
  box->stress5+=stress[3]+stresk[3];
  box->stress6+=stress[4]+stresk[4];
  box->stress7+=stress[2]+stresk[2];
  box->stress8+=stress[4]+stresk[4];
  box->stress9+=stress[5]+stresk[5];
  
  simulCond->lambdat=lambdc;
  
  ener->conint+=simulCond->timeStep*lambda*qmass/X2(simulCond->taut);
  ener->consv=ener->conint+0.5*qmass*X2(lambda);
  
// periodic boundary condition
  
  image_update(atom,simulCond,box);
  
// updated velocity
  
  #ifdef _OPENMP
  #pragma omp parallel for default(none) shared(simulCond,vxu,vyu,vzu,atom) private(i)
  #endif
  for(i=0;i<simulCond->natom;i++)
  {
    
    atom[i].vx=vxu[i];
    atom[i].vy=vyu[i];
    atom[i].vz=vzu[i];
    
  }
  
// Free temporary arrays
  
  free(vxu);
  free(vyu);
  free(vzu);
  
  free(xo);
  free(yo);
  free(zo);
  
  free(vxo);
  free(vyo);
  free(vzo);
  
  if(simulCond->nconst>0)
  {
    free(dd);
    
    free(xt);
    free(yt);
    free(zt);
  }
      
}

void lf_npt_h(ATOM atom[], ENERGY *ener, SIMULPARAMS *simulCond,CONSTRAINT *constList,PBC *box)
{
  int i,k,ia,ib,nosecycle;
  double lambda,lambdb,lambdc,ts2,qmass;
  double gamma,gammb,gammc,pmass;
  double *xo=NULL,*yo=NULL,*zo=NULL;
  double *vxo=NULL,*vyo=NULL,*vzo=NULL;
  double *xt=NULL,*yt=NULL,*zt=NULL;
  double *vxu=NULL,*vyu=NULL,*vzu=NULL;
  double volume,masst=0.,cell0[9],com[3]={0.},vom[3]={0.};
  double virshake=0.,virshakt=0.,stress[6]={0.},strest[6]={0.},stresk[6]={0.};
  DELTA *dd=NULL;
  
  vxu=(double*)malloc(simulCond->natom*sizeof(*vxu));
  vyu=(double*)malloc(simulCond->natom*sizeof(*vyu));
  vzu=(double*)malloc(simulCond->natom*sizeof(*vzu));
  
  xo=(double*)malloc(simulCond->natom*sizeof(*xo));
  yo=(double*)malloc(simulCond->natom*sizeof(*yo));
  zo=(double*)malloc(simulCond->natom*sizeof(*zo));
  
  vxo=(double*)malloc(simulCond->natom*sizeof(*vxo));
  vyo=(double*)malloc(simulCond->natom*sizeof(*vyo));
  vzo=(double*)malloc(simulCond->natom*sizeof(*vzo));
    
  
  #ifdef _OPENMP
  #pragma omp parallel for default(none) shared(simulCond,xo,yo,zo,vxo,vyo,vzo,atom) private(i)
  #endif
  for(i=0;i<simulCond->natom;i++)
  {
    
// Store old coordinates and old velocities.

    xo[i]=atom[i].x;
    yo[i]=atom[i].y;
    zo[i]=atom[i].z;
    
    vxo[i]=atom[i].vx;
    vyo[i]=atom[i].vy;
    vzo[i]=atom[i].vz;
    
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
  
  if(simulCond->nconst>0)
  {
    dd=(DELTA*)malloc(simulCond->nconst*sizeof(*dd));
    
    #ifdef _OPENMP
    #pragma omp parallel for default(none) shared(simulCond,atom,dd,constList) private(i,ia,ib)
    #endif
    for(i=0;i<simulCond->nconst;i++)
    {
      ia=constList[i].a;
      ib=constList[i].b;
      
      dd[i].x=atom[ib].x-atom[ia].x;
      dd[i].y=atom[ib].y-atom[ia].y;
      dd[i].z=atom[ib].z-atom[ia].z;
    }
    
    image_array(simulCond->nconst,dd,simulCond,box);
    
    xt=(double*)malloc(simulCond->natom*sizeof(*xt));
    yt=(double*)malloc(simulCond->natom*sizeof(*yt));
    zt=(double*)malloc(simulCond->natom*sizeof(*zt));
    
  }
  
  // total mass and center of mass
  
  masst=0.;
  com[0]=0.;
  com[1]=0.;
  com[2]=0.;
  for(i=0;i<simulCond->natom;i++)
  {
    masst+=atom[i].m;
    com[0]+=atom[i].m*atom[i].x;
    com[1]+=atom[i].m*atom[i].y;
    com[2]+=atom[i].m*atom[i].z;
  }
  com[0]/=masst;
  com[1]/=masst;
  com[2]/=masst;
  
  //   Mass parameter for Nose-Hoover thermostat
  
  qmass=2.0*simulCond->kintemp0*X2(simulCond->taut);
  pmass=2.0*simulCond->kintemp0*X2(simulCond->taup);
  
  ts2=X2(simulCond->timeStep);
  
  #ifdef _OPENMP
  #pragma omp parallel for default(none) shared(simulCond,atom) private(i)
  #endif
  for(i=0;i<simulCond->natom;i++)
  { 
    atom[i].vx+=0.5*simulCond->timeStep*atom[i].fx/atom[i].m;
    atom[i].vy+=0.5*simulCond->timeStep*atom[i].fy/atom[i].m;
    atom[i].vz+=0.5*simulCond->timeStep*atom[i].fz/atom[i].m;
  }
  
  ener->kin=kinetic(atom,simulCond);
  
  gammb=(2.0*ener->kin - ener->virpot - virshake - 3.0*simulCond->press*volume)/pmass-
    simulCond->lambdat*simulCond->gammap;
  gammc=simulCond->gammap+simulCond->timeStep*gammb;
  gamma=0.5*(simulCond->gammap+gammc)
  
  lambdb=(2.0*(ener->kin - simulCond->kintemp0 ) + pmass*X2(simulCond->gammap)-
    rboltzui*simulCond->temp)/qmass;
  lambdc=simulCond->lambdat+simulCond->timeStep*lambdb;
  lambda=0.5*(simulCond->lambdat+lambdc);
  
  if(simulCond->nconst>0)
    nosecycle=4;
  else
    nosecycle=5;
  
  for(k=0;k<nosecycle;k++)
  {
    
// move atoms by leapfrog algorithm
    
    #ifdef _OPENMP
    #pragma omp parallel for default(none) shared(simulCond,xt,yt,zt,xo,yo,zo,vxo,vyo,vzo,vxu,vyu,vzu,atom,lambda) private(i)
    #endif
    for(i=0;i<simulCond->natom;i++)
    {
      
// update velocities
      
      vxu[i]=vxo[i]+simulCond->timeStep*(atom[i].fx/atom[i].m-atom[i].vx*(lambda+gamma));
      vyu[i]=vyo[i]+simulCond->timeStep*(atom[i].fy/atom[i].m-atom[i].vy*(lambda+gamma));
      vzu[i]=vzo[i]+simulCond->timeStep*(atom[i].fz/atom[i].m-atom[i].vz*(lambda+gamma));
      
// update positions
      
      atom[i].x=xo[i]+simulCond->timeStep*(vxu[i]+gammc*(0.5*(atom[i].x+xo[i])-com[0]));
      atom[i].y=yo[i]+simulCond->timeStep*(vyu[i]+gammc*(0.5*(atom[i].y+yo[i])-com[1]));
      atom[i].z=zo[i]+simulCond->timeStep*(vzu[i]+gammc*(0.5*(atom[i].z+zo[i])-com[2]));
      
// Temporary storage of the uncorrected positions
      
      if(simulCond->nconst>0)
      {
	xt[i]=atom[i].x;
	yt[i]=atom[i].y;
	zt[i]=atom[i].z;
      }
      
    }
    
    if(simulCond->nconst>0)
    {
// Apply constraint with Shake algorithm.
      
      cbrga=exp(3.*simulCond->timeStep*gammc);
      cbrga=cbrt(cbrga);
      
      scale_box(box,cbrga,cell0);
      
      lf_shake(atom,simulCond,constList,dd,box,&virshakt,strest);
      
      virshake+=virshakt;
      for(i=0;i<6;i++)
	stress[i]+=strest[i];
      
      #ifdef _OPENMP
      #pragma omp parallel for default(none) shared(simulCond,xt,yt,zt,vxu,vyu,vzu,atom,ts2) private(i)
      #endif
      for(i=0;i<simulCond->natom;i++)
      {
        
// Corrected velocities
      
	vxu[i]+=(atom[i].x-xt[i])/simulCond->timeStep;
	vyu[i]+=(atom[i].y-yt[i])/simulCond->timeStep;
	vzu[i]+=(atom[i].z-zt[i])/simulCond->timeStep;
      
// Corrected Forces
      
	atom[i].fx+=(atom[i].x-xt[i])*atom[i].m/ts2;
	atom[i].fy+=(atom[i].y-yt[i])*atom[i].m/ts2;
	atom[i].fz+=(atom[i].z-zt[i])*atom[i].m/ts2;
      
      }
    }
    
// calculate full timestep velocity

    #ifdef _OPENMP
    #pragma omp parallel for default(none) shared(simulCond,vxo,vyo,vzo,vxu,vyu,vzu,atom) private(i)
    #endif
    for(i=0;i<simulCond->natom;i++)
    {
      
      atom[i].vx=0.5*(vxo[i]+vxu[i]);
      atom[i].vy=0.5*(vyo[i]+vyu[i]);
      atom[i].vz=0.5*(vzo[i]+vzu[i]);
      
    }
    
// calculate kinetic energy
    
    ener->kin=kinetic(atom,simulCond);
    
    gammb=(2.0*ener->kin-ener->virpot-virshake-3.0*simulCond->press*volume)/pmass-
      simulCond->lambdat*simulCond->gammap;
    gammc=simulCond->gammap+simulCond->timeStep*gammb;
    gamma=0.5*(simulCond->gammap+gammc)
    
    lambdb=(2.0*(ener->kin-simulCond->kintemp0)+pmass*X2(simulCond->gammap)-
      rboltzui*simulCond->temp)/qmass;
    lambdc=simulCond->lambdat+simulCond->timeStep*lambdb;
    lambda=0.5*(simulCond->lambdat+lambdc);

  }
  
  ener->virshake=virshake;
  
  stress_kinetic(atom,simulCond,stresk);
  
  box->stress1+=stress[0]+stresk[0];
  box->stress2+=stress[1]+stresk[1];
  box->stress3+=stress[2]+stresk[2];
  box->stress4+=stress[1]+stresk[1];
  box->stress5+=stress[3]+stresk[3];
  box->stress6+=stress[4]+stresk[4];
  box->stress7+=stress[2]+stresk[2];
  box->stress8+=stress[4]+stresk[4];
  box->stress9+=stress[5]+stresk[5];
  
  cbrga=exp(3.*simulCond->timeStep*gammc);
  cbrga=cbrt(cbrga);
  
  scale_box(box,cbrga,cell0);
  
  simulCond->lambdat=lambdc;
  simulCond->gammap=gammc;
  
  ener->conint+=simulCond->timeStep*lambda*(rboltzui*simulCond->temp+qmass/X2(simulCond->taut));
  ener->consv=ener->conint+simulCond->press*volume+0.5*(qmass*X2(lambda)+pmass*X2(gamma));
  
// periodic boundary condition
  
  image_update(atom,simulCond,box);
  
// updated velocity
  
  #ifdef _OPENMP
  #pragma omp parallel for default(none) shared(simulCond,vxu,vyu,vzu,atom) private(i)
  #endif
  for(i=0;i<simulCond->natom;i++)
  {
    
    atom[i].vx=vxu[i];
    atom[i].vy=vyu[i];
    atom[i].vz=vzu[i];
    
  }
  
  vom[0]=0.;
  vom[1]=0.;
  vom[2]=0.;
  for(i=0;i<simulCond->natom;i++)
  {
    vom[0]+=atom[i].m*atom[i].vx;
    vom[1]+=atom[i].m*atom[i].vy;
    vom[2]+=atom[i].m*atom[i].vz;
  }
  vom[0]/=masst;
  vom[1]/=masst;
  vom[2]/=masst;
  
  for(i=0;i<simulCond->natom;i++)
  {
    atom[i].vx-=vom[0];
    atom[i].vy-=vom[1];
    atom[i].vz-=vom[2];
  }
  
// Free temporary arrays
  
  free(vxu);
  free(vyu);
  free(vzu);
  
  free(xo);
  free(yo);
  free(zo);
  
  free(vxo);
  free(vyo);
  free(vzo);
  
  if(simulCond->nconst>0)
  {
    free(dd);
    
    free(xt);
    free(yt);
    free(zt);
  }
      
}

void vv_integrate(ATOM atom[], ENERGY *ener, SIMULPARAMS *simulCond,CONSTRAINT *constList,PBC *box,int stage)
{
  switch (simulCond->ens)
  {
    case 0:
      vv_nve(atom,ener,simulCond,constList,box,stage);
      break;
    case 1:
      vv_nvt_b(atom,ener,simulCond,constList,box,stage);
      break;
    case 2:
      vv_npt_b(atom,ener,simulCond,constList,box,stage);
      break;
    case 3:
      vv_nvt_h(atom,ener,simulCond,constList,box,stage);
      break;
    case 4:
      vv_npt_h(atom,ener,simulCond,constList,box,stage);
      break;
    default:
      vv_nve(atom,ener,simulCond,constList,box,stage);
      break;
  }
}

void vv_nve(ATOM atom[], ENERGY *ener, SIMULPARAMS *simulCond,CONSTRAINT *constList,PBC *box,int stage)
{
  int i,ia,ib;
  double virshake,stress[6]={0.},stresk[6]={0.};
  DELTA *dd=NULL;
  
  if(simulCond->nconst>0)
  {
    dd=(DELTA*)malloc(simulCond->nconst*sizeof(*dd));
    
    #ifdef _OPENMP
    #pragma omp parallel for default(none) shared(simulCond,constList,dd,atom) private(i,ia,ib)
    #endif
    for(i=0;i<simulCond->nconst;i++)
    {
      ia=constList[i].a;
      ib=constList[i].b;
      
      dd[i].x=atom[ib].x-atom[ia].x;
      dd[i].y=atom[ib].y-atom[ia].y;
      dd[i].z=atom[ib].z-atom[ia].z;
    }
    
    image_array(simulCond->nconst,dd,simulCond,box);
    
  }

// move atoms by leapfrog algorithm
  
  #ifdef _OPENMP
  #pragma omp parallel for default(none) shared(simulCond,atom) private(i)
  #endif
  for(i=0;i<simulCond->natom;i++)
  {
// update velocities
    
    atom[i].vx+=0.5*simulCond->timeStep*atom[i].fx/atom[i].m;
    atom[i].vy+=0.5*simulCond->timeStep*atom[i].fy/atom[i].m;
    atom[i].vz+=0.5*simulCond->timeStep*atom[i].fz/atom[i].m;
  }
  
  if(stage==1)
  {
    #ifdef _OPENMP
    #pragma omp parallel for default(none) shared(simulCond,atom) private(i)
    #endif
    for(i=0;i<simulCond->natom;i++)
    {
// update positions
      
      atom[i].x+=simulCond->timeStep*atom[i].vx;
      atom[i].y+=simulCond->timeStep*atom[i].vy;
      atom[i].z+=simulCond->timeStep*atom[i].vz;
      
    }
    
    if(simulCond->nconst>0)
    {
      
// Apply constraint with Shake algorithm.

      vv_shake_r(atom,simulCond,constList,dd,box,&virshake,stress);
      ener->virshake=virshake;
      
    }
    
  }
  else
  {
// calculate kinetic energy

    if(simulCond->nconst>0)
    {
      
// Apply constraint with Shake algorithm.

      vv_shake_v(atom,simulCond,constList,dd);
      
    }
  
    ener->kin=kinetic(atom,simulCond);
  
    stress_kinetic(atom,simulCond,stresk);
    
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
    
    image_update(atom,simulCond,box);
  }
  
  if(simulCond->nconst>0)
  {
    free(dd);
  }
   
}

void vv_nvt_b(ATOM atom[], ENERGY *ener, SIMULPARAMS *simulCond,CONSTRAINT *constList,PBC *box,int stage)
{
  int i,ia,ib;
  double lambda;
  double virshake,stress[6]={0.},stresk[6]={0.};
  DELTA *dd=NULL;
  
  if(simulCond->nconst>0)
  {
    dd=(DELTA*)malloc(simulCond->nconst*sizeof(*dd));
    
    #ifdef _OPENMP
    #pragma omp parallel for default(none) shared(simulCond,constList,dd,atom) private(i,ia,ib)
    #endif
    for(i=0;i<simulCond->nconst;i++)
    {
      ia=constList[i].a;
      ib=constList[i].b;
      
      dd[i].x=atom[ib].x-atom[ia].x;
      dd[i].y=atom[ib].y-atom[ia].y;
      dd[i].z=atom[ib].z-atom[ia].z;
    }
    
    image_array(simulCond->nconst,dd,simulCond,box);
    
  }

// move atoms by leapfrog algorithm
  
  #ifdef _OPENMP
  #pragma omp parallel for default(none) shared(simulCond,atom) private(i)
  #endif
  for(i=0;i<simulCond->natom;i++)
  {
// update velocities
    
    atom[i].vx+=0.5*simulCond->timeStep*atom[i].fx/atom[i].m;
    atom[i].vy+=0.5*simulCond->timeStep*atom[i].fy/atom[i].m;
    atom[i].vz+=0.5*simulCond->timeStep*atom[i].fz/atom[i].m;
  }
  
  if(stage==1)
  {
    #ifdef _OPENMP
    #pragma omp parallel for default(none) shared(simulCond,atom) private(i)
    #endif
    for(i=0;i<simulCond->natom;i++)
    {
// update positions
      
      atom[i].x+=simulCond->timeStep*atom[i].vx;
      atom[i].y+=simulCond->timeStep*atom[i].vy;
      atom[i].z+=simulCond->timeStep*atom[i].vz;
      
    }
    
    if(simulCond->nconst>0)
    {
      
// Apply constraint with Shake algorithm.

      vv_shake_r(atom,simulCond,constList,dd,box,&virshake,stress);
      ener->virshake=virshake;
      
    }
    
  }
  else
  {
// calculate kinetic energy

    if(simulCond->nconst>0)
    {
      
// Apply constraint with Shake algorithm.

      vv_shake_v(atom,simulCond,constList,dd);
      
    }
  
    ener->kin=kinetic(atom,simulCond);
    
    lambda=sqrt(1.0+simulCond->timeStep/simulCond->taut*(simulCond->kintemp0/ener->kin-1.0));
    
    #ifdef _OPENMP
    #pragma omp parallel for default(none) shared(simulCond,atom,lambda) private(i)
    #endif
    for(i=0;i<simulCond->natom;i++)
    {
      atom[i].vx*=lambda;
      atom[i].vy*=lambda;
      atom[i].vz*=lambda;
    }
    
    ener->kin*=X2(lambda);
  
    stress_kinetic(atom,simulCond,stresk);
    
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
    
    image_update(atom,simulCond,box);
  }
  
  if(simulCond->nconst>0)
  {
    free(dd);
  }
   
}

void vv_npt_b(ATOM atom[], ENERGY *ener, SIMULPARAMS *simulCond,CONSTRAINT *constList,PBC *box,int stage)
{
  int i,ia,ib,nosecycle;
  double lambda,gamma,cbrga,pp;
  double *xo=NULL,*yo=NULL,*zo=NULL;
  double *vxo=NULL,*vyo=NULL,*vzo=NULL;
  double virshake,stress[6]={0.},stresk[6]={0.};
  DELTA *dd=NULL;
  
  xo=(double*)malloc(simulCond->natom*sizeof(*xo));
  yo=(double*)malloc(simulCond->natom*sizeof(*yo));
  zo=(double*)malloc(simulCond->natom*sizeof(*zo));
  
  vxo=(double*)malloc(simulCond->natom*sizeof(*vxo));
  vyo=(double*)malloc(simulCond->natom*sizeof(*vyo));
  vzo=(double*)malloc(simulCond->natom*sizeof(*vzo));
  
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
  
  if(simulCond->nconst>0)
  {
    dd=(DELTA*)malloc(simulCond->nconst*sizeof(*dd));
    
    #ifdef _OPENMP
    #pragma omp parallel for default(none) shared(simulCond,constList,dd,atom) private(i,ia,ib)
    #endif
    for(i=0;i<simulCond->nconst;i++)
    {
      ia=constList[i].a;
      ib=constList[i].b;
      
      dd[i].x=atom[ib].x-atom[ia].x;
      dd[i].y=atom[ib].y-atom[ia].y;
      dd[i].z=atom[ib].z-atom[ia].z;
    }
    
    image_array(simulCond->nconst,dd,simulCond,box);
    
  }
  
  if(stage==1)
  {
    
    ener->kin=kinetic(atom,simulCond);
    
    #ifdef _OPENMP
    #pragma omp parallel for default(none) shared(simulCond,atom) private(i)
    #endif
    for(i=0;i<simulCond->natom;i++)
    {
//    update velocities
      
      atom[i].vx+=0.5*simulCond->timeStep*atom[i].fx/atom[i].m;
      atom[i].vy+=0.5*simulCond->timeStep*atom[i].fy/atom[i].m;
      atom[i].vz+=0.5*simulCond->timeStep*atom[i].fz/atom[i].m;
    }
    
    if(simulCond->nconst>0)
    {
      #ifdef _OPENMP
      #pragma omp parallel for default(none) shared(simulCond,xo,yo,zo,vxo,vyo,vzo,atom) private(i)
      #endif
      for(i=0;i<simulCond->natom;i++)
      {
	
    // Store old coordinates and old velocities.

	xo[i]=atom[i].x;
	yo[i]=atom[i].y;
	zo[i]=atom[i].z;
	
	vxo[i]=atom[i].vx;
	vyo[i]=atom[i].vy;
	vzo[i]=atom[i].vz;
	
      }
    }
    
    if( (stage==1) && (simulCond->nconst>0) )
      nosecycle=2;
    else
      nosecycle=1;
    
    for(k=0;k<nosecycle;k++)
    {
      cbrga=1.;
      
      if(k==nosecycle-1)
      {
	pp=(2.*ener->kin-ener->virpot-virshake)/(3.*volume);
	gamma=1.+watercomp*simulCond->timeStep*(pp-simulCond->press)/simulCond->taup;
	cbrga=cbrt(gamma);
	
	vv_scale_box(box,cbrga);
	
      }
      
      #ifdef _OPENMP
      #pragma omp parallel for default(none) shared(simulCond,atom) private(i)
      #endif
      for(i=0;i<simulCond->natom;i++)
      {
  // update positions
	
	atom[i].x=cbrga*atom[i].x+simulCond->timeStep*atom[i].vx;
	atom[i].y=cbrga*atom[i].y+simulCond->timeStep*atom[i].vy;
	atom[i].z=cbrga*atom[i].z+simulCond->timeStep*atom[i].vz;
	
      }
      
      if(simulCond->nconst>0)
      {
	
  // Apply constraint with Shake algorithm.

	vv_shake_r(atom,simulCond,constList,dd,box,&virshake,stress);
	ener->virshake=virshake;
	
      }
      
      if(k<nosecycle-1)
      {
	#ifdef _OPENMP
	#pragma omp parallel for default(none) shared(simulCond,xo,yo,zo,vxo,vyo,vzo,atom) private(i)
	#endif
	for(i=0;i<simulCond->natom;i++)
	{
	  
      // Store old coordinates and old velocities.

	  atom[i].x=xo[i];
	  atom[i].y=yo[i];
	  atom[i].z=zo[i];
	  
	  atom[i].vx=vxo[i];
	  atom[i].vy=vyo[i];
	  atom[i].vz=vzo[i];
	  
	}
      }
    }
  }
  else
  {
    
    #ifdef _OPENMP
    #pragma omp parallel for default(none) shared(simulCond,atom) private(i)
    #endif
    for(i=0;i<simulCond->natom;i++)
    {
//    update velocities
      
      atom[i].vx+=0.5*simulCond->timeStep*atom[i].fx/atom[i].m;
      atom[i].vy+=0.5*simulCond->timeStep*atom[i].fy/atom[i].m;
      atom[i].vz+=0.5*simulCond->timeStep*atom[i].fz/atom[i].m;
    }
    
// calculate kinetic energy

    ener->kin=kinetic(atom,simulCond);
    
    lambda=sqrt(1.0+simulCond->timeStep/simulCond->taut*(simulCond->kintemp0/ener->kin-1.0));
    
    for(i=0;i<simulCond->natom;i++)
    {
//    update velocities
      
      atom[i].vx*=lambda;
      atom[i].vy*=lambda;
      atom[i].vz*=lambda;
    }

    if(simulCond->nconst>0)
    {
      
// Apply constraint with Shake algorithm.

      vv_shake_v(atom,simulCond,constList,dd);
      
    }
    
    ener->kin=kinetic(atom,simulCond);
  
    stress_kinetic(atom,simulCond,stresk);
    
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
    
    image_update(atom,simulCond,box);
  }
  
  free(xo);
  free(yo);
  free(zo);
  
  free(vxo);
  free(vyo);
  free(vzo)
  
  if(simulCond->nconst>0)
  {
    free(dd);
  }
   
}

void vv_nvt_h(ATOM atom[], ENERGY *ener, SIMULPARAMS *simulCond,CONSTRAINT *constList,PBC *box,int stage)
{
  int i,ia,ib;
  double lambda,qmass;
  double virshake,stress[6]={0.},stresk[6]={0.};
  DELTA *dd=NULL;
  
  if(simulCond->nconst>0)
  {
    dd=(DELTA*)malloc(simulCond->nconst*sizeof(*dd));
    
    #ifdef _OPENMP
    #pragma omp parallel for default(none) shared(simulCond,constList,dd,atom) private(i,ia,ib)
    #endif
    for(i=0;i<simulCond->nconst;i++)
    {
      ia=constList[i].a;
      ib=constList[i].b;
      
      dd[i].x=atom[ib].x-atom[ia].x;
      dd[i].y=atom[ib].y-atom[ia].y;
      dd[i].z=atom[ib].z-atom[ia].z;
    }
    
    image_array(simulCond->nconst,dd,simulCond,box);
    
  }
  
  qmass=2.0*simulCond->kintemp0*X2(simulCond->taut);
  
  if(stage==1)
  {
    
    ener->kin=kinetic(atom,simulCond);
    
    simulCond->lambdat+=0.5*simulCond->timeStep*(ener->kin-simulCond->kintemp0)/qmass;
    
    lambda=exp(-0.5*simulCond->timeStep*simulCond->lambdat);
    
    #ifdef _OPENMP
    #pragma omp parallel for default(none) shared(simulCond,atom,lambda) private(i)
    #endif
    for(i=0;i<simulCond->natom;i++)
    {
  // scale velocities
      
      atom[i].vx*=lambda;
      atom[i].vy*=lambda;
      atom[i].vz*=lambda;
      
  // update velocities
      
      atom[i].vx+=0.5*simulCond->timeStep*atom[i].fx/atom[i].m;
      atom[i].vy+=0.5*simulCond->timeStep*atom[i].fy/atom[i].m;
      atom[i].vz+=0.5*simulCond->timeStep*atom[i].fz/atom[i].m;
    }
    
    ener->kin*=X2(lambda);
    
    ener->conint+=0.5*simulCond->timeStep*simulCond->lambdat*qmass/X2(simulCond->taut);
    
    simulCond->lambdat+=0.5*simulCond->timeStep*(ener->kin-simulCond->kintemp0)/qmass;
    
    #ifdef _OPENMP
    #pragma omp parallel for default(none) shared(simulCond,atom) private(i)
    #endif
    for(i=0;i<simulCond->natom;i++)
    {
// update positions
      
      atom[i].x+=simulCond->timeStep*atom[i].vx;
      atom[i].y+=simulCond->timeStep*atom[i].vy;
      atom[i].z+=simulCond->timeStep*atom[i].vz;
      
    }
    
    if(simulCond->nconst>0)
    {
      
// Apply constraint with Shake algorithm.

      vv_shake_r(atom,simulCond,constList,dd,box,&virshake,stress);
      ener->virshake=virshake;
      
    }
    
  }
  else
  {
// calculate kinetic energy

    #ifdef _OPENMP
    #pragma omp parallel for default(none) shared(simulCond,atom) private(i)
    #endif
    for(i=0;i<simulCond->natom;i++)
    {
  // update velocities
      
      atom[i].vx+=0.5*simulCond->timeStep*atom[i].fx/atom[i].m;
      atom[i].vy+=0.5*simulCond->timeStep*atom[i].fy/atom[i].m;
      atom[i].vz+=0.5*simulCond->timeStep*atom[i].fz/atom[i].m;
    }

    if(simulCond->nconst>0)
    {
      
// Apply constraint with Shake algorithm.

      vv_shake_v(atom,simulCond,constList,dd);
      
    }
  
    ener->kin=kinetic(atom,simulCond);
    
    simulCond->lambdat+=0.5*simulCond->timeStep*(ener->kin-simulCond->kintemp0)/qmass;
    
    lambda=exp(-0.5*simulCond->timeStep*simulCond->lambdat);
    
    #ifdef _OPENMP
    #pragma omp parallel for default(none) shared(simulCond,atom,lambda) private(i)
    #endif
    for(i=0;i<simulCond->natom;i++)
    {
      atom[i].vx*=lambda;
      atom[i].vy*=lambda;
      atom[i].vz*=lambda;
    }
    
    ener->kin*=X2(lambda);
    
    ener->conint+=0.5*simulCond->timeStep*simulCond->lambdat*qmass/X2(simulCond->taut);
    
    simulCond->lambdat+=0.5*simulCond->timeStep*(ener->kin-simulCond->kintemp0)/qmass;
    
    ener->consv=ener->conint+0.5*qmass*X2(simulCond->lambdat);
    
    stress_kinetic(atom,simulCond,stresk);
    
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
    
    image_update(atom,simulCond,box);
  }
  
  if(simulCond->nconst>0)
  {
    free(dd);
  }
   
}

void vv_npt_h(ATOM atom[], ENERGY *ener, SIMULPARAMS *simulCond,CONSTRAINT *constList,PBC *box,int stage)
{
  int i,ia,ib,nosecycle,hoovercycle=5;
  double hts,chts,cqts;
  double cons0,lambda,lambda0,qmass;
  double gamma,gamma0,pmass,scale;
  double *xo=NULL,*yo=NULL,*zo=NULL;
  double *vxo=NULL,*vyo=NULL,*vzo=NULL;
  double volume,volume0,cell0[9],masst=0.,com[3]={0.},vom[3]={0.};
  double virshake,stress[6]={0.},stresk[6]={0.};
  DELTA *dd=NULL;
  
  xo=(double*)malloc(simulCond->natom*sizeof(*xo));
  yo=(double*)malloc(simulCond->natom*sizeof(*yo));
  zo=(double*)malloc(simulCond->natom*sizeof(*zo));
  
  vxo=(double*)malloc(simulCond->natom*sizeof(*vxo));
  vyo=(double*)malloc(simulCond->natom*sizeof(*vyo));
  vzo=(double*)malloc(simulCond->natom*sizeof(*vzo));
  
  hts=0.5*simulCond->timeStep;
  chts=hts/(double)hoovercycle;
  cqts=0.25*simulCond->timeStep/(double)hoovercycle;
  
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
  for(i=0;i<simulCond->natom;i++)
    masst+=atom[i].m;
  
  qmass=2.0*simulCond->kintemp0*X2(simulCond->taut);
  pmass=2.0*simulCond->kintemp0*X2(simulCond->taup);
  
  if(simulCond->nconst>0)
  {
    dd=(DELTA*)malloc(simulCond->nconst*sizeof(*dd));
    
    #ifdef _OPENMP
    #pragma omp parallel for default(none) shared(simulCond,constList,dd,atom) private(i,ia,ib)
    #endif
    for(i=0;i<simulCond->nconst;i++)
    {
      ia=constList[i].a;
      ib=constList[i].b;
      
      dd[i].x=atom[ib].x-atom[ia].x;
      dd[i].y=atom[ib].y-atom[ia].y;
      dd[i].z=atom[ib].z-atom[ia].z;
    }
    
    image_array(simulCond->nconst,dd,simulCond,box);
    
  }
  
  if(stage==1)
  {
    
    if(simulCond->nconst>0)
    {
      
      lambda0=simulCond->lambdat;
      gamma0=simulCond->gammap;
      cons0=ener->conint;
      
      #ifdef _OPENMP
      #pragma omp parallel for default(none) shared(simulCond,xo,yo,zo,vxo,vyo,vzo,atom) private(i)
      #endif
      for(i=0;i<simulCond->natom;i++)
      {
	
    // Store old coordinates and old velocities.

	xo[i]=atom[i].x;
	yo[i]=atom[i].y;
	zo[i]=atom[i].z;
	
	vxo[i]=atom[i].vx;
	vyo[i]=atom[i].vy;
	vzo[i]=atom[i].vz;
	
      }
      
    }
    
    if( (stage==1) && (simulCond->nconst>0) )
      nosecycle=2;
    else
      nosecycle=1;
    
    for(k=0;k<nosecycle;k++)
    {
      
      for(kk=0;kk<hoovercycle;kk++)
      {
	
      // apply nvt
	ener->kin=kinetic(atom,simulCond);
	
	simulCond->lambdat+=0.5*cqst*(2.0*(ener->kin-simulCond->kintemp0)+
	  pmass*X2(simulCond->gammap)-rboltzui*simulCond->temp)/qmass;
	
	lambda=exp(-cqst*simulCond->lambdat);
	
	#ifdef _OPENMP
	#pragma omp parallel for default(none) shared(simulCond,atom,lambda) private(i)
	#endif
	for(i=0;i<simulCond->natom;i++)
	{
      // scale velocities
	  
	  atom[i].vx*=lambda;
	  atom[i].vy*=lambda;
	  atom[i].vz*=lambda;
	  
	}
	
	ener->kin*=X2(lambda);
	
	ener->conint+=cqst*simulCond->lambdat*(rboltzui*simulCond->temp+qmass/X2(simulCond->taut));
	
	simulCond->lambdat+=0.5*cqst*(2.0*(ener->kin-simulCond->kintemp0)+
	  pmass*X2(simulCond->gammap)-rboltzui*simulCond->temp)/qmass;
	  
      // apply npt
	  
	simulCond->gammap+=0.5*chst*(((2.0*ener->kin-ener->virpot-virshake)-
	  3.0*simulCond->press*volume)/pmass-simulCond->gammap*simulCond->lambdat);
	
	gamma=exp(-chst*simulCond->gammap);
	
	for(i=0;i<simulCond->natom;i++)
	{
      // scale velocities
	  
	  atom[i].vx*=gamma;
	  atom[i].vy*=gamma;
	  atom[i].vz*=gamma;
	  
	}
	
	ener->kin*=X2(gamma);
	
	volume*=exp(3.0*chst*simulCond->gammap);
	
	simulCond->gammap+=0.5*chst*(((2.0*ener->kin-ener->virpot-virshake)-
	  3.0*simulCond->press*volume)/pmass-simulCond->gammap*simulCond->lambdat);
	
      // apply nvt
	
	ener->kin=kinetic(atom,simulCond);
	
	simulCond->lambdat+=0.5*cqst*(2.0*(ener->kin-simulCond->kintemp0)+
	  pmass*X2(simulCond->gammap)-rboltzui*simulCond->temp)/qmass;
	
	lambda=exp(-cqst*simulCond->lambdat);
	
	#ifdef _OPENMP
	#pragma omp parallel for default(none) shared(simulCond,atom,lambda) private(i)
	#endif
	for(i=0;i<simulCond->natom;i++)
	{
      // scale velocities
	  
	  atom[i].vx*=lambda;
	  atom[i].vy*=lambda;
	  atom[i].vz*=lambda;
	  
	}
	
	ener->kin*=X2(lambda);
	
	ener->conint+=cqst*simulCond->lambdat*(rboltzui*simulCond->temp+qmass/X2(simulCond->taut));
	
	simulCond->lambdat+=0.5*cqst*(2.0*(ener->kin-simulCond->kintemp0)+
	  pmass*X2(simulCond->gammap)-rboltzui*simulCond->temp)/qmass;
      }
      
      scale=cbrt(volume/volume0);
      scale_box(box,scale,cell0);
      
      #ifdef _OPENMP
      #pragma omp parallel for default(none) shared(simulCond,atom,lambda) private(i)
      #endif
      for(i=0;i<simulCond->natom;i++)
      {
    // update velocities
	
	atom[i].vx+=0.5*simulCond->timeStep*atom[i].fx/atom[i].m;
	atom[i].vy+=0.5*simulCond->timeStep*atom[i].fy/atom[i].m;
	atom[i].vz+=0.5*simulCond->timeStep*atom[i].fz/atom[i].m;
      }
      
      com[0]=0.;
      com[1]=0.;
      com[2]=0.;
      for(i=0;i<simulCond->natom;i++)
      {
	com[0]+=atom[i].m*atom[i].x;
	com[1]+=atom[i].m*atom[i].y;
	com[2]+=atom[i].m*atom[i].z;
      }
      com[0]/=masst;
      com[1]/=masst;
      com[2]/=masst;
      
      cbrga=exp(simulCond->timeStep*simulCond->gammap);
      
      #ifdef _OPENMP
      #pragma omp parallel for default(none) shared(simulCond,atom) private(i)
      #endif
      for(i=0;i<simulCond->natom;i++)
      {
  // update positions
	
	atom[i].x=cbrga*(atom[i].x-com[0])+simulCond->timeStep*atom[i].vx+com[0];
	atom[i].y=cbrga*(atom[i].y-com[1])+simulCond->timeStep*atom[i].vy+com[1];
	atom[i].z=cbrga*(atom[i].z-com[2])+simulCond->timeStep*atom[i].vz+com[2];
	
      }
    
      if(simulCond->nconst>0)
      {
	
  // Apply constraint with Shake algorithm.

	vv_shake_r(atom,simulCond,constList,dd,box,&virshake,stress);
	ener->virshake=virshake;
	
      }
      
      if(k<nosecycle-1)
      {
	volume=volume0;
	simulCond->lambdat=lambda0;
	simulCond->gammap=gamma0;
	ener->conint=cons0;
	
	#ifdef _OPENMP
	#pragma omp parallel for default(none) shared(simulCond,xo,yo,zo,vxo,vyo,vzo,atom) private(i)
	#endif
	for(i=0;i<simulCond->natom;i++)
	{
	  
      // Store old coordinates and old velocities.

	  xo[i]=atom[i].x;
	  yo[i]=atom[i].y;
	  zo[i]=atom[i].z;
	  
	  vxo[i]=atom[i].vx;
	  vyo[i]=atom[i].vy;
	  vzo[i]=atom[i].vz;
	  
	}
	
      }
    }
  }
  else
  {
// calculate kinetic energy

    #ifdef _OPENMP
    #pragma omp parallel for default(none) shared(simulCond,atom) private(i)
    #endif
    for(i=0;i<simulCond->natom;i++)
    {
  // update velocities
      
      atom[i].vx+=0.5*simulCond->timeStep*atom[i].fx/atom[i].m;
      atom[i].vy+=0.5*simulCond->timeStep*atom[i].fy/atom[i].m;
      atom[i].vz+=0.5*simulCond->timeStep*atom[i].fz/atom[i].m;
    }

    if(simulCond->nconst>0)
    {
      
// Apply constraint with Shake algorithm.

      vv_shake_v(atom,simulCond,constList,dd);
      
    }
    
    for(kk=0;kk<hoovercycle;kk++)
    {
      
    // apply nvt
      ener->kin=kinetic(atom,simulCond);
      
      simulCond->lambdat+=0.5*cqst*(2.0*(ener->kin-simulCond->kintemp0)+
	pmass*X2(simulCond->gammap)-rboltzui*simulCond->temp)/qmass;
      
      lambda=exp(-cqst*simulCond->lambdat);
      
      #ifdef _OPENMP
      #pragma omp parallel for default(none) shared(simulCond,atom,lambda) private(i)
      #endif
      for(i=0;i<simulCond->natom;i++)
      {
    // scale velocities
	
	atom[i].vx*=lambda;
	atom[i].vy*=lambda;
	atom[i].vz*=lambda;
	
      }
      
      ener->kin*=X2(lambda);
      
      ener->conint+=cqst*simulCond->lambdat*(rboltzui*simulCond->temp+qmass/X2(simulCond->taut));
      
      simulCond->lambdat+=0.5*cqst*(2.0*(ener->kin-simulCond->kintemp0)+
	pmass*X2(simulCond->gammap)-rboltzui*simulCond->temp)/qmass;
	
    // apply npt
	
      simulCond->gammap+=0.5*chst*(((2.0*ener->kin-ener->virpot-virshake)-
	3.0*simulCond->press*volume)/pmass-simulCond->gammap*simulCond->lambdat);
      
      gamma=exp(-chst*simulCond->gammap);
      
      for(i=0;i<simulCond->natom;i++)
      {
    // scale velocities
	
	atom[i].vx*=gamma;
	atom[i].vy*=gamma;
	atom[i].vz*=gamma;
	
      }
      
      ener->kin*=X2(gamma);
      
      volume*=exp(3.0*chst*simulCond->gammap);
      
      simulCond->gammap+=0.5*chst*(((2.0*ener->kin-ener->virpot-virshake)-
	3.0*simulCond->press*volume)/pmass-simulCond->gammap*simulCond->lambdat);
      
    // apply nvt
      
      ener->kin=kinetic(atom,simulCond);
      
      simulCond->lambdat+=0.5*cqst*(2.0*(ener->kin-simulCond->kintemp0)+
	pmass*X2(simulCond->gammap)-rboltzui*simulCond->temp)/qmass;
      
      lambda=exp(-cqst*simulCond->lambdat);
      
      #ifdef _OPENMP
      #pragma omp parallel for default(none) shared(simulCond,atom,lambda) private(i)
      #endif
      for(i=0;i<simulCond->natom;i++)
      {
    // scale velocities
	
	atom[i].vx*=lambda;
	atom[i].vy*=lambda;
	atom[i].vz*=lambda;
	
      }
      
      ener->kin*=X2(lambda);
      
      ener->conint+=cqst*simulCond->lambdat*(rboltzui*simulCond->temp+qmass/X2(simulCond->taut));
      
      simulCond->lambdat+=0.5*cqst*(2.0*(ener->kin-simulCond->kintemp0)+
	pmass*X2(simulCond->gammap)-rboltzui*simulCond->temp)/qmass;
    }
    
    vom[0]=0.;
    vom[1]=0.;
    vom[2]=0.;
    for(i=0;i<simulCond->natom;i++)
    {
      vom[0]+=atom[i].m*atom[i].vx;
      vom[1]+=atom[i].m*atom[i].vy;
      vom[2]+=atom[i].m*atom[i].vz;
    }
    vom[0]/=masst;
    vom[1]/=masst;
    vom[2]/=masst;
    
    for(i=0;i<simulCond->natom;i++)
    {
      atom[i].vx-=vom[0];
      atom[i].vy-=vom[1];
      atom[i].vz-=vom[2];
    }
    
    scale=cbrt(volume/volume0);
    scale_box(box,scale,cell0);
    
    ener->consv=ener->conint+simulCond->press*volume+0.5*(qmass*X2(simulCond->lambdat)+pmass*X2(simulCond->gammat));
    
    ener->kin=kinetic(atom,simulCond);
    
    stress_kinetic(atom,simulCond,stresk);
    
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
    
    image_update(atom,simulCond,box);
  }
  
  free(xo);
  free(yo);
  free(zo);
  
  free(vxo);
  free(vyo);
  free(vzo)
  
  if(simulCond->nconst>0)
  {
    free(dd);
  }
   
}