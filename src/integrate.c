#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "global.h"
#include "utils.h"
#include "shake.h"
#include "integrate.h"

void lf_integrate(ATOM *atom, ENERGYFORCE *enerFor, SIMULPARAMS *simulCond,CONSTRAINT *constList)
{
  switch (simulCond->ens)
  {
    case 0:
      lf_nve(atom,enerFor,simulCond,constList);
      break;
    case 1:
      lf_nvt_b(atom,enerFor,simulCond,constList);
      break;
    default:
      lf_nve(atom,enerFor,simulCond,constList);
      break;
  }
}

void lf_nve(ATOM *atom, ENERGYFORCE *enerFor, SIMULPARAMS *simulCond,CONSTRAINT *constList)
{
  int i,ia,ib;
  double *xo=NULL,*yo=NULL,*zo=NULL,*vxu=NULL,*vyu=NULL,*vzu=NULL;
  DELTA *dd=NULL;
  
  vxu=(double*)malloc(atom->natom*sizeof(*vxu));
  vyu=(double*)malloc(atom->natom*sizeof(*vyu));
  vzu=(double*)malloc(atom->natom*sizeof(*vzu));
  
  if(simulCond->nconst>0)
  {
    dd=(DELTA*)malloc(simulCond->nconst*sizeof(*dd));
    
    for(i=0;i<simulCond->nconst;i++)
    {
      ia=constList[i].a;
      ib=constList[i].b;
      
      dd[i].x=atom->x[ib]-atom->x[ia];
      dd[i].y=atom->y[ib]-atom->y[ia];
      dd[i].z=atom->z[ib]-atom->z[ia];
    }
    
    image_array(simulCond->nconst,dd,simulCond);
    
    xo=(double*)malloc(atom->natom*sizeof(*xo));
    yo=(double*)malloc(atom->natom*sizeof(*yo));
    zo=(double*)malloc(atom->natom*sizeof(*zo));
    
    for(i=0;i<atom->natom;i++)
    {
      
// Store old coordinates.

      xo[i]=atom->x[i];
      yo[i]=atom->y[i];
      zo[i]=atom->z[i];
      
    }
    
  }

// move atoms by leapfrog algorithm
  
  for(i=0;i<atom->natom;i++)
  {
    
// update velocities
    
    vxu[i]=atom->vx[i]+simulCond->timeStep*atom->fx[i]/atom->m[i];
    vyu[i]=atom->vy[i]+simulCond->timeStep*atom->fy[i]/atom->m[i];
    vzu[i]=atom->vz[i]+simulCond->timeStep*atom->fz[i]/atom->m[i];
    
// update positions
    
    atom->x[i]+=simulCond->timeStep*vxu[i];
    atom->y[i]+=simulCond->timeStep*vyu[i];
    atom->z[i]+=simulCond->timeStep*vzu[i];
    
  }
  
  if(simulCond->nconst>0)
  {
// Apply constraint with Shake algorithm.

    lf_shake(atom,simulCond,constList,dd);
    for(i=0;i<atom->natom;i++)
    {
        
// Corrected velocities
    
      vxu[i]=(atom->x[i]-xo[i])/simulCond->timeStep;
      vyu[i]=(atom->y[i]-yo[i])/simulCond->timeStep;
      vzu[i]=(atom->z[i]-zo[i])/simulCond->timeStep;
    
// Corrected Forces
    
      atom->fx[i]=(vxu[i]-atom->vx[i])*atom->m[i]/simulCond->timeStep;
      atom->fy[i]=(vyu[i]-atom->vy[i])*atom->m[i]/simulCond->timeStep;
      atom->fz[i]=(vzu[i]-atom->vz[i])*atom->m[i]/simulCond->timeStep;
    
    }
  }
  
// calculate full timestep velocity

  for(i=0;i<atom->natom;i++)
  {
    
    atom->vx[i]=0.5*(atom->vx[i]+vxu[i]);
    atom->vy[i]=0.5*(atom->vy[i]+vyu[i]);
    atom->vz[i]=0.5*(atom->vz[i]+vzu[i]);
    
  }
  
// calculate kinetic energy
  
  enerFor->energyKin=kinetic(atom);
  
// periodic boundary condition
  
  image_update(atom,simulCond);
  
// updated velocity
  
  for(i=0;i<atom->natom;i++)
  {
    
    atom->vx[i]=vxu[i];
    atom->vy[i]=vyu[i];
    atom->vz[i]=vzu[i];
    
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

void lf_nvt_b(ATOM *atom, ENERGYFORCE *enerFor, SIMULPARAMS *simulCond,CONSTRAINT *constList)
{
  int i,k,ia,ib,bercycle;
  double lambda,ts2;
  double *xo=NULL,*yo=NULL,*zo=NULL;
  double *vxo=NULL,*vyo=NULL,*vzo=NULL;
  double *xt=NULL,*yt=NULL,*zt=NULL;
  double *vxu=NULL,*vyu=NULL,*vzu=NULL;
  DELTA *dd=NULL;
  
  vxu=(double*)malloc(atom->natom*sizeof(*vxu));
  vyu=(double*)malloc(atom->natom*sizeof(*vyu));
  vzu=(double*)malloc(atom->natom*sizeof(*vzu));
  
  xo=(double*)malloc(atom->natom*sizeof(*xo));
  yo=(double*)malloc(atom->natom*sizeof(*yo));
  zo=(double*)malloc(atom->natom*sizeof(*zo));
  
  vxo=(double*)malloc(atom->natom*sizeof(*vxo));
  vyo=(double*)malloc(atom->natom*sizeof(*vyo));
  vzo=(double*)malloc(atom->natom*sizeof(*vzo));
    
  for(i=0;i<atom->natom;i++)
  {
    
// Store old coordinates and old velocities.

    xo[i]=atom->x[i];
    yo[i]=atom->y[i];
    zo[i]=atom->z[i];
    
    vxo[i]=atom->vx[i];
    vyo[i]=atom->vy[i];
    vzo[i]=atom->vz[i];
    
  }
  
  if(simulCond->nconst>0)
  {
    dd=(DELTA*)malloc(simulCond->nconst*sizeof(*dd));
    
    for(i=0;i<simulCond->nconst;i++)
    {
      ia=constList[i].a;
      ib=constList[i].b;
      
      dd[i].x=atom->x[ib]-atom->x[ia];
      dd[i].y=atom->y[ib]-atom->y[ia];
      dd[i].z=atom->z[ib]-atom->z[ia];
    }
    
    image_array(simulCond->nconst,dd,simulCond);
    
    xt=(double*)malloc(atom->natom*sizeof(*xt));
    yt=(double*)malloc(atom->natom*sizeof(*yt));
    zt=(double*)malloc(atom->natom*sizeof(*zt));
    
  }
  
  ts2=X2(simulCond->timeStep);
  
  for(i=0;i<atom->natom;i++)
  { 
    atom->vx[i]+=0.5*simulCond->timeStep*atom->fx[i]/atom->m[i];
    atom->vy[i]+=0.5*simulCond->timeStep*atom->fy[i]/atom->m[i];
    atom->vz[i]+=0.5*simulCond->timeStep*atom->fz[i]/atom->m[i];
  }
  
  enerFor->energyKin=kinetic(atom);
  
  if(simulCond->nconst>0)
    bercycle=2;
  else
    bercycle=3;
  
  for(k=0;k<bercycle;k++)
  {
   
    lambda=sqrt(1.0+simulCond->timeStep/simulCond->taut*(simulCond->kintemp0/enerFor->energyKin-1.0));
    
// move atoms by leapfrog algorithm
    
    for(i=0;i<atom->natom;i++)
    {
      
// update velocities
      
      vxu[i]=(vxo[i]+simulCond->timeStep*atom->fx[i]/atom->m[i])*lambda;
      vyu[i]=(vyo[i]+simulCond->timeStep*atom->fy[i]/atom->m[i])*lambda;
      vzu[i]=(vzo[i]+simulCond->timeStep*atom->fz[i]/atom->m[i])*lambda;
      
// update positions
      
      atom->x[i]=xo[i]+simulCond->timeStep*vxu[i];
      atom->y[i]=yo[i]+simulCond->timeStep*vyu[i];
      atom->z[i]=zo[i]+simulCond->timeStep*vzu[i];
      
// Temporary storage of the uncorrected positions
      
      if(simulCond->nconst>0)
      {
	xt[i]=atom->x[i];
	yt[i]=atom->y[i];
	zt[i]=atom->z[i];
      }
      
    }
    
    if(simulCond->nconst>0)
    {
// Apply constraint with Shake algorithm.
      
      lf_shake(atom,simulCond,constList,dd);
      for(i=0;i<atom->natom;i++)
      {
        
// Corrected velocities
      
	vxu[i]+=(atom->x[i]-xt[i])/simulCond->timeStep;
	vyu[i]+=(atom->y[i]-yt[i])/simulCond->timeStep;
	vzu[i]+=(atom->z[i]-zt[i])/simulCond->timeStep;
      
// Corrected Forces
      
	atom->fx[i]+=(atom->x[i]-xt[i])*atom->m[i]/ts2;
	atom->fy[i]+=(atom->y[i]-yt[i])*atom->m[i]/ts2;
	atom->fz[i]+=(atom->z[i]-zt[i])*atom->m[i]/ts2;
      
      }
    }
    
// calculate full timestep velocity

    for(i=0;i<atom->natom;i++)
    {
      
      atom->vx[i]=0.5*(vxo[i]+vxu[i]);
      atom->vy[i]=0.5*(vyo[i]+vyu[i]);
      atom->vz[i]=0.5*(vzo[i]+vzu[i]);
      
    }
    
// calculate kinetic energy
    
    enerFor->energyKin=kinetic(atom);
    
  }
  
// periodic boundary condition
  
  image_update(atom,simulCond);
  
// updated velocity
  
  for(i=0;i<atom->natom;i++)
  {
    
    atom->vx[i]=vxu[i];
    atom->vy[i]=vyu[i];
    atom->vz[i]=vzu[i];
    
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


void vv_integrate(ATOM *atom, ENERGYFORCE *enerFor, SIMULPARAMS *simulCond,CONSTRAINT *constList,int stage)
{
  switch (simulCond->ens)
  {
    case 0:
      vv_nve(atom,enerFor,simulCond,stage);
      break;
    default:
      vv_nve(atom,enerFor,simulCond,stage);
      break;
  }
}

void vv_nve(ATOM *atom, ENERGYFORCE *enerFor, SIMULPARAMS *simulCond,int stage)
{
  int i;

// move atoms by leapfrog algorithm
  
  for(i=0;i<atom->natom;i++)
  {
// update velocities
    
    atom->vx[i]+=0.5*simulCond->timeStep*atom->fx[i]/atom->m[i];
    atom->vy[i]+=0.5*simulCond->timeStep*atom->fy[i]/atom->m[i];
    atom->vz[i]+=0.5*simulCond->timeStep*atom->fz[i]/atom->m[i];
  }
  
  if(stage==1)
  {
    for(i=0;i<atom->natom;i++)
    {
// update positions
      
      atom->x[i]+=simulCond->timeStep*atom->vx[i];
      atom->y[i]+=simulCond->timeStep*atom->vy[i];
      atom->z[i]+=simulCond->timeStep*atom->vz[i];
      
    }
  }
  else
  {
// calculate kinetic energy
  
    enerFor->energyKin=kinetic(atom);
  }
  
  if(stage==2)
  {
    
// periodic boundary condition
    
    image_update(atom,simulCond);
  }
      
}
