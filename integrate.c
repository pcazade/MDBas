#include <stdio.h>
#include <stdlib.h>
#include "global.h"
#include "utils.h"

void lf_nve(ATOM *atom, ENERGYFORCE *enerFor, SIMULPARAMS *simulCond)
{
  int i;
  double *vxu,*vyu,*vzu;
  
  vxu=(double*)malloc(atom->natom*sizeof(*vxu));
  vyu=(double*)malloc(atom->natom*sizeof(*vyu));
  vzu=(double*)malloc(atom->natom*sizeof(*vzu));

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
