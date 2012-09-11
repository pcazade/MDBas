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
    case 3:
      lf_nvt_h(atom,ener,simulCond,constList,box);
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

    lf_shake(atom,simulCond,constList,dd,box);
    
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
      
      lf_shake(atom,simulCond,constList,dd,box);
      
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
      
      lf_shake(atom,simulCond,constList,dd,box);
      
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
    case 3:
      vv_nvt_h(atom,ener,simulCond,constList,box,stage);
      break;
    default:
      vv_nve(atom,ener,simulCond,constList,box,stage);
      break;
  }
}

void vv_nve(ATOM atom[], ENERGY *ener, SIMULPARAMS *simulCond,CONSTRAINT *constList,PBC *box,int stage)
{
  int i,ia,ib;
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

      vv_shake_r(atom,simulCond,constList,dd,box);
      
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

      vv_shake_r(atom,simulCond,constList,dd,box);
      
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

void vv_nvt_h(ATOM atom[], ENERGY *ener, SIMULPARAMS *simulCond,CONSTRAINT *constList,PBC *box,int stage)
{
  int i,ia,ib;
  double lambda,qmass;
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

      vv_shake_r(atom,simulCond,constList,dd,box);
      
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
