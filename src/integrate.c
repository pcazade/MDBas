#include <stdio.h>
#include <stdlib.h>
#include "global.h"
#include "utils.h"
#include "shake.h"

void lf_nve(ATOM *atom, ENERGYFORCE *enerFor, SIMULPARAMS *simulCond,CONSTRAINT *constList)
{
  int i,ia,ib;
  double *xo,*yo,*zo,*vxu,*vyu,*vzu;
  DELTA *dd;
  
  xo=(double*)malloc(atom->natom*sizeof(*xo));
  yo=(double*)malloc(atom->natom*sizeof(*yo));
  zo=(double*)malloc(atom->natom*sizeof(*zo));
  
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
  }

// move atoms by leapfrog algorithm
  
  for(i=0;i<atom->natom;i++)
  {
    
// Store old coordinates.

    xo[i]=atom->x[i];
    yo[i]=atom->y[i];
    zo[i]=atom->z[i];
    
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
