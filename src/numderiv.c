#include <stdlib.h>

#include "global.h"
#include "energy.h"

void numforce(ATOM atom[],DELTA *nForce,FORCEFIELD *ff,ENERGY *ener,SIMULPARAMS *simulCond,PBC *box,int npoints,double h)
{
  int i;
  double coord;
  
  if(npoints==2)
  {
    
    double a,b;
    
    for(i=0;i<simulCond->natom;i++)
    {
      
      coord=atom[i].x;
      
      atom[i].x=coord-h;
      energy(atom,ff,ener,simulCond,box);
      a=ener->pot;
      
      atom[i].x=coord+h;
      energy(atom,ff,ener,simulCond,box);
      b=ener->pot;
      
      nForce[i].x=-(b-a)/(2.*h);
      
      atom[i].x=coord;
      
      /***************************************************/
      
      coord=atom[i].y;
      
      atom[i].y=coord-h;
      energy(atom,ff,ener,simulCond,box);
      a=ener->pot;
      
      atom[i].y=coord+h;
      energy(atom,ff,ener,simulCond,box);
      b=ener->pot;
      
      nForce[i].y=-(b-a)/(2.*h);
      
      atom[i].y=coord;
      
      /***************************************************/
      
      coord=atom[i].z;
      
      atom[i].z=coord-h;
      energy(atom,ff,ener,simulCond,box);
      a=ener->pot;
      
      atom[i].z=coord+h;
      energy(atom,ff,ener,simulCond,box);
      b=ener->pot;
      
      nForce[i].z=-(b-a)/(2.*h);
      
      atom[i].z=coord;
    }
  }
  else if(npoints==4)
  {
    
    double a,b,c,d;
    
    for(i=0;i<simulCond->natom;i++)
    {
      
      coord=atom[i].x;
      
      atom[i].x=coord-(2.*h);
      energy(atom,ff,ener,simulCond,box);
      a=ener->pot;
      
      atom[i].x=coord-h;
      energy(atom,ff,ener,simulCond,box);
      b=ener->pot;
      
      atom[i].x=coord+h;
      energy(atom,ff,ener,simulCond,box);
      c=ener->pot;
      
      atom[i].x=coord+(2.*h);
      energy(atom,ff,ener,simulCond,box);
      d=ener->pot;
      
      nForce[i].x=-(8.*(c-b)-(d-a))/(12.*h);
      
      atom[i].x=coord;
      
      /***************************************************/
      
      coord=atom[i].y;
      
      atom[i].y=coord-(2.*h);
      energy(atom,ff,ener,simulCond,box);
      a=ener->pot;
      
      atom[i].y=coord-h;
      energy(atom,ff,ener,simulCond,box);
      b=ener->pot;
      
      atom[i].y=coord+h;
      energy(atom,ff,ener,simulCond,box);
      c=ener->pot;
      
      atom[i].y=coord+(2.*h);
      energy(atom,ff,ener,simulCond,box);
      d=ener->pot;
      
      nForce[i].y=-(8.*(c-b)-(d-a))/(12.*h);
      
      atom[i].y=coord;
      
      /***************************************************/
      
      coord=atom[i].z;
      
      atom[i].z=coord-(2.*h);
      energy(atom,ff,ener,simulCond,box);
      a=ener->pot;
      
      atom[i].z=coord-h;
      energy(atom,ff,ener,simulCond,box);
      b=ener->pot;
      
      atom[i].z=coord+h;
      energy(atom,ff,ener,simulCond,box);
      c=ener->pot;
      
      atom[i].z=coord+(2.*h);
      energy(atom,ff,ener,simulCond,box);
      d=ener->pot;
      
      nForce[i].z=-(8.*(c-b)-(d-a))/(12.*h);
      
      atom[i].z=coord;
      
    }
  }
  else if(npoints==6)
  {
    
    double a,b,c,d,e,f;
    
    for(i=0;i<simulCond->natom;i++)
    {
      
      coord=atom[i].x;
      
      atom[i].x=coord-(3.*h);
      energy(atom,ff,ener,simulCond,box);
      a=ener->pot;
      
      atom[i].x=coord-(2.*h);
      energy(atom,ff,ener,simulCond,box);
      b=ener->pot;
      
      atom[i].x=coord-h;
      energy(atom,ff,ener,simulCond,box);
      c=ener->pot;
      
      atom[i].x=coord+h;
      energy(atom,ff,ener,simulCond,box);
      d=ener->pot;
      
      atom[i].x=coord+(2.*h);
      energy(atom,ff,ener,simulCond,box);
      e=ener->pot;
      
      atom[i].x=coord+(3.*h);
      energy(atom,ff,ener,simulCond,box);
      f=ener->pot;
      
      nForce[i].x=-(45.*(d-c)-9.*(e-b)+(f-a))/(60.*h);
      
      atom[i].x=coord;
      
      /***************************************************/
      
      coord=atom[i].y;
      
      atom[i].y=coord-(3.*h);
      energy(atom,ff,ener,simulCond,box);
      a=ener->pot;
      
      atom[i].y=coord-(2.*h);
      energy(atom,ff,ener,simulCond,box);
      b=ener->pot;
      
      atom[i].y=coord-h;
      energy(atom,ff,ener,simulCond,box);
      c=ener->pot;
      
      atom[i].y=coord+h;
      energy(atom,ff,ener,simulCond,box);
      d=ener->pot;
      
      atom[i].y=coord+(2.*h);
      energy(atom,ff,ener,simulCond,box);
      e=ener->pot;
      
      atom[i].y=coord+(3.*h);
      energy(atom,ff,ener,simulCond,box);
      f=ener->pot;
      
      nForce[i].y=-(45.*(d-c)-9.*(e-b)+(f-a))/(60.*h);
      
      atom[i].y=coord;
      
      /***************************************************/
      
      coord=atom[i].z;
      
      atom[i].z=coord-(3.*h);
      energy(atom,ff,ener,simulCond,box);
      a=ener->pot;
      
      atom[i].z=coord-(2.*h);
      energy(atom,ff,ener,simulCond,box);
      b=ener->pot;
      
      atom[i].z=coord-h;
      energy(atom,ff,ener,simulCond,box);
      c=ener->pot;
      
      atom[i].z=coord+h;
      energy(atom,ff,ener,simulCond,box);
      d=ener->pot;
      
      atom[i].z=coord+(2.*h);
      energy(atom,ff,ener,simulCond,box);
      e=ener->pot;
      
      atom[i].z=coord+(3.*h);
      energy(atom,ff,ener,simulCond,box);
      f=ener->pot;
      
      nForce[i].z=-(45.*(d-c)-9.*(e-b)+(f-a))/(60.*h);
      
      atom[i].z=coord;
    }
  }
  
}
