#include <stdlib.h>
#include "global.h"
#include "energy.h"

void numforce(ATOM *atom,FORCEFIELD *ff,ENERGYFORCE *enerFor,SIMULPARAMS *simulCond,int npoints,double h)
{
  int i;
  double coord;
  
  atom->numFx=(double*)malloc(atom->natom*sizeof(*(atom->numFx)));
  atom->numFy=(double*)malloc(atom->natom*sizeof(*(atom->numFy)));
  atom->numFz=(double*)malloc(atom->natom*sizeof(*(atom->numFz)));
  
  if(npoints==2)
  {
    
    double a,b;
    
    for(i=0;i<atom->natom;i++)
    {
      
      coord=atom->x[i];
      
      atom->x[i]=coord-h;
      energy(atom,ff,enerFor,simulCond);
      a=enerFor->energyPot;
      
      atom->x[i]=coord+h;
      energy(atom,ff,enerFor,simulCond);
      b=enerFor->energyPot;
      
      atom->numFx[i]=-(b-a)/(2.*h);
      
      atom->x[i]=coord;
      
      /***************************************************/
      
      coord=atom->y[i];
      
      atom->y[i]=coord-h;
      energy(atom,ff,enerFor,simulCond);
      a=enerFor->energyPot;
      
      atom->y[i]=coord+h;
      energy(atom,ff,enerFor,simulCond);
      b=enerFor->energyPot;
      
      atom->numFy[i]=-(b-a)/(2.*h);
      
      atom->y[i]=coord;
      
      /***************************************************/
      
      coord=atom->z[i];
      
      atom->z[i]=coord-h;
      energy(atom,ff,enerFor,simulCond);
      a=enerFor->energyPot;
      
      atom->z[i]=coord+h;
      energy(atom,ff,enerFor,simulCond);
      b=enerFor->energyPot;
      
      atom->numFz[i]=-(b-a)/(2.*h);
      
      atom->z[i]=coord;
    }
  }
  else if(npoints==4)
  {
    
    double a,b,c,d;
    
    for(i=0;i<atom->natom;i++)
    {
      
      coord=atom->x[i];
      
      atom->x[i]=coord-(2.*h);
      energy(atom,ff,enerFor,simulCond);
      a=enerFor->energyPot;
      
      atom->x[i]=coord-h;
      energy(atom,ff,enerFor,simulCond);
      b=enerFor->energyPot;
      
      atom->x[i]=coord+h;
      energy(atom,ff,enerFor,simulCond);
      c=enerFor->energyPot;
      
      atom->x[i]=coord+(2.*h);
      energy(atom,ff,enerFor,simulCond);
      d=enerFor->energyPot;
      
      atom->numFx[i]=-(8.*(c-b)-(d-a))/(12.*h);
      
      atom->x[i]=coord;
      
      /***************************************************/
      
      coord=atom->y[i];
      
      atom->y[i]=coord-(2.*h);
      energy(atom,ff,enerFor,simulCond);
      a=enerFor->energyPot;
      
      atom->y[i]=coord-h;
      energy(atom,ff,enerFor,simulCond);
      b=enerFor->energyPot;
      
      atom->y[i]=coord+h;
      energy(atom,ff,enerFor,simulCond);
      c=enerFor->energyPot;
      
      atom->y[i]=coord+(2.*h);
      energy(atom,ff,enerFor,simulCond);
      d=enerFor->energyPot;
      
      atom->numFy[i]=-(8.*(c-b)-(d-a))/(12.*h);
      
      atom->y[i]=coord;
      
      /***************************************************/
      
      coord=atom->z[i];
      
      atom->z[i]=coord-(2.*h);
      energy(atom,ff,enerFor,simulCond);
      a=enerFor->energyPot;
      
      atom->z[i]=coord-h;
      energy(atom,ff,enerFor,simulCond);
      b=enerFor->energyPot;
      
      atom->z[i]=coord+h;
      energy(atom,ff,enerFor,simulCond);
      c=enerFor->energyPot;
      
      atom->z[i]=coord+(2.*h);
      energy(atom,ff,enerFor,simulCond);
      d=enerFor->energyPot;
      
      atom->numFz[i]=-(8.*(c-b)-(d-a))/(12.*h);
      
      atom->z[i]=coord;
      
    }
  }
  else if(npoints==6)
  {
    
    double a,b,c,d,e,f;
    
    for(i=0;i<atom->natom;i++)
    {
      
      coord=atom->x[i];
      
      atom->x[i]=coord-(3.*h);
      energy(atom,ff,enerFor,simulCond);
      a=enerFor->energyPot;
      
      atom->x[i]=coord-(2.*h);
      energy(atom,ff,enerFor,simulCond);
      b=enerFor->energyPot;
      
      atom->x[i]=coord-h;
      energy(atom,ff,enerFor,simulCond);
      c=enerFor->energyPot;
      
      atom->x[i]=coord+h;
      energy(atom,ff,enerFor,simulCond);
      d=enerFor->energyPot;
      
      atom->x[i]=coord+(2.*h);
      energy(atom,ff,enerFor,simulCond);
      e=enerFor->energyPot;
      
      atom->x[i]=coord+(3.*h);
      energy(atom,ff,enerFor,simulCond);
      f=enerFor->energyPot;
      
      atom->numFx[i]=-(45.*(d-c)-9.*(e-b)+(f-a))/(60.*h);
      
      atom->x[i]=coord;
      
      /***************************************************/
      
      coord=atom->y[i];
      
      atom->y[i]=coord-(3.*h);
      energy(atom,ff,enerFor,simulCond);
      a=enerFor->energyPot;
      
      atom->y[i]=coord-(2.*h);
      energy(atom,ff,enerFor,simulCond);
      b=enerFor->energyPot;
      
      atom->y[i]=coord-h;
      energy(atom,ff,enerFor,simulCond);
      c=enerFor->energyPot;
      
      atom->y[i]=coord+h;
      energy(atom,ff,enerFor,simulCond);
      d=enerFor->energyPot;
      
      atom->y[i]=coord+(2.*h);
      energy(atom,ff,enerFor,simulCond);
      e=enerFor->energyPot;
      
      atom->y[i]=coord+(3.*h);
      energy(atom,ff,enerFor,simulCond);
      f=enerFor->energyPot;
      
      atom->numFy[i]=-(45.*(d-c)-9.*(e-b)+(f-a))/(60.*h);
      
      atom->y[i]=coord;
      
      /***************************************************/
      
      coord=atom->z[i];
      
      atom->z[i]=coord-(3.*h);
      energy(atom,ff,enerFor,simulCond);
      a=enerFor->energyPot;
      
      atom->z[i]=coord-(2.*h);
      energy(atom,ff,enerFor,simulCond);
      b=enerFor->energyPot;
      
      atom->z[i]=coord-h;
      energy(atom,ff,enerFor,simulCond);
      c=enerFor->energyPot;
      
      atom->z[i]=coord+h;
      energy(atom,ff,enerFor,simulCond);
      d=enerFor->energyPot;
      
      atom->z[i]=coord+(2.*h);
      energy(atom,ff,enerFor,simulCond);
      e=enerFor->energyPot;
      
      atom->z[i]=coord+(3.*h);
      energy(atom,ff,enerFor,simulCond);
      f=enerFor->energyPot;
      
      atom->numFz[i]=-(45.*(d-c)-9.*(e-b)+(f-a))/(60.*h);
      
      atom->z[i]=coord;
    }
  }
  
}
