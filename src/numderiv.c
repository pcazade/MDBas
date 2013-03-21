/**
 * \file numderiv.c
 * \brief Contains functions for estimating derivatives numerically (test purposes only).
 * \author Pierre-Andre Cazade and Florent Hedin
 * \version alpha-branch
 * \date 2012
 */

#include <stdlib.h>

#include "global.h"
#include "energy.h"

void numforce(CTRL *ctrl,PARAM *param,ENERGY *ener,PBC *box,NEIGH *neigh,
	      BOND bond[],BOND ub[],ANGLE angle[],DIHE dihe[],DIHE impr[],
	      DELTA nForce[],double x[],double y[],double z[],int npoints,double h)
{
  int i;
  double coord;
  
  if(npoints==2)
  {
    
    double a,b;
    
    for(i=0;i<param->nAtom;i++)
    {
      
      coord=x[i];
      
      x[i]=coord-h;
      //energy(ctrl,param,ener,box,neigh,atom,bond,ub,angle,dihe,impr,vdw);
      a=ener->pot;
      
      x[i]=coord+h;
      //energy(ctrl,param,ener,box,neigh,atom,bond,ub,angle,dihe,impr,vdw);
      b=ener->pot;
      
      nForce[i].x=-(b-a)/(2.*h);
      
      x[i]=coord;
      
      /***************************************************/
      
      coord=y[i];
      
      y[i]=coord-h;
      //energy(ctrl,param,ener,box,neigh,atom,bond,ub,angle,dihe,impr,vdw);
      a=ener->pot;
      
      y[i]=coord+h;
      //energy(ctrl,param,ener,box,neigh,atom,bond,ub,angle,dihe,impr,vdw);
      b=ener->pot;
      
      nForce[i].y=-(b-a)/(2.*h);
      
      y[i]=coord;
      
      /***************************************************/
      
      coord=z[i];
      
      z[i]=coord-h;
      //energy(ctrl,param,ener,box,neigh,atom,bond,ub,angle,dihe,impr,vdw);
      a=ener->pot;
      
      z[i]=coord+h;
      //energy(ctrl,param,ener,box,neigh,atom,bond,ub,angle,dihe,impr,vdw);
      b=ener->pot;
      
      nForce[i].z=-(b-a)/(2.*h);
      
      z[i]=coord;
    }
  }
  else if(npoints==4)
  {
    
    double a,b,c,d;
    
    for(i=0;i<param->nAtom;i++)
    {
      
      coord=x[i];
      
      x[i]=coord-(2.*h);
      //energy(ctrl,param,ener,box,neigh,atom,bond,ub,angle,dihe,impr,vdw);
      a=ener->pot;
      
      x[i]=coord-h;
      //energy(ctrl,param,ener,box,neigh,atom,bond,ub,angle,dihe,impr,vdw);
      b=ener->pot;
      
      x[i]=coord+h;
      //energy(ctrl,param,ener,box,neigh,atom,bond,ub,angle,dihe,impr,vdw);
      c=ener->pot;
      
      x[i]=coord+(2.*h);
      //energy(ctrl,param,ener,box,neigh,atom,bond,ub,angle,dihe,impr,vdw);
      d=ener->pot;
      
      nForce[i].x=-(8.*(c-b)-(d-a))/(12.*h);
      
      x[i]=coord;
      
      /***************************************************/
      
      coord=y[i];
      
      y[i]=coord-(2.*h);
      //energy(ctrl,param,ener,box,neigh,atom,bond,ub,angle,dihe,impr,vdw);
      a=ener->pot;
      
      y[i]=coord-h;
      //energy(ctrl,param,ener,box,neigh,atom,bond,ub,angle,dihe,impr,vdw);
      b=ener->pot;
      
      y[i]=coord+h;
      //energy(ctrl,param,ener,box,neigh,atom,bond,ub,angle,dihe,impr,vdw);
      c=ener->pot;
      
      y[i]=coord+(2.*h);
      //energy(ctrl,param,ener,box,neigh,atom,bond,ub,angle,dihe,impr,vdw);
      d=ener->pot;
      
      nForce[i].y=-(8.*(c-b)-(d-a))/(12.*h);
      
      y[i]=coord;
      
      /***************************************************/
      
      coord=z[i];
      
      z[i]=coord-(2.*h);
      //energy(ctrl,param,ener,box,neigh,atom,bond,ub,angle,dihe,impr,vdw);
      a=ener->pot;
      
      z[i]=coord-h;
      //energy(ctrl,param,ener,box,neigh,atom,bond,ub,angle,dihe,impr,vdw);
      b=ener->pot;
      
      z[i]=coord+h;
      //energy(ctrl,param,ener,box,neigh,atom,bond,ub,angle,dihe,impr,vdw);
      c=ener->pot;
      
      z[i]=coord+(2.*h);
      //energy(ctrl,param,ener,box,neigh,atom,bond,ub,angle,dihe,impr,vdw);
      d=ener->pot;
      
      nForce[i].z=-(8.*(c-b)-(d-a))/(12.*h);
      
      z[i]=coord;
      
    }
  }
  else if(npoints==6)
  {
    
    double a,b,c,d,e,f;
    
    for(i=0;i<param->nAtom;i++)
    {
      
      coord=x[i];
      
      x[i]=coord-(3.*h);
      //energy(ctrl,param,ener,box,neigh,atom,bond,ub,angle,dihe,impr,vdw);
      a=ener->pot;
      
      x[i]=coord-(2.*h);
      //energy(ctrl,param,ener,box,neigh,atom,bond,ub,angle,dihe,impr,vdw);
      b=ener->pot;
      
      x[i]=coord-h;
      //energy(ctrl,param,ener,box,neigh,atom,bond,ub,angle,dihe,impr,vdw);
      c=ener->pot;
      
      x[i]=coord+h;
      //energy(ctrl,param,ener,box,neigh,atom,bond,ub,angle,dihe,impr,vdw);
      d=ener->pot;
      
      x[i]=coord+(2.*h);
      //energy(ctrl,param,ener,box,neigh,atom,bond,ub,angle,dihe,impr,vdw);
      e=ener->pot;
      
      x[i]=coord+(3.*h);
      //energy(ctrl,param,ener,box,neigh,atom,bond,ub,angle,dihe,impr,vdw);
      f=ener->pot;
      
      nForce[i].x=-(45.*(d-c)-9.*(e-b)+(f-a))/(60.*h);
      
      x[i]=coord;
      
      /***************************************************/
      
      coord=y[i];
      
      y[i]=coord-(3.*h);
      //energy(ctrl,param,ener,box,neigh,atom,bond,ub,angle,dihe,impr,vdw);
      a=ener->pot;
      
      y[i]=coord-(2.*h);
      //energy(ctrl,param,ener,box,neigh,atom,bond,ub,angle,dihe,impr,vdw);
      b=ener->pot;
      
      y[i]=coord-h;
      //energy(ctrl,param,ener,box,neigh,atom,bond,ub,angle,dihe,impr,vdw);
      c=ener->pot;
      
      y[i]=coord+h;
      //energy(ctrl,param,ener,box,neigh,atom,bond,ub,angle,dihe,impr,vdw);
      d=ener->pot;
      
      y[i]=coord+(2.*h);
      //energy(ctrl,param,ener,box,neigh,atom,bond,ub,angle,dihe,impr,vdw);
      e=ener->pot;
      
      y[i]=coord+(3.*h);
      //energy(ctrl,param,ener,box,neigh,atom,bond,ub,angle,dihe,impr,vdw);
      f=ener->pot;
      
      nForce[i].y=-(45.*(d-c)-9.*(e-b)+(f-a))/(60.*h);
      
      y[i]=coord;
      
      /***************************************************/
      
      coord=z[i];
      
      z[i]=coord-(3.*h);
      //energy(ctrl,param,ener,box,neigh,atom,bond,ub,angle,dihe,impr,vdw);
      a=ener->pot;
      
      z[i]=coord-(2.*h);
      //energy(ctrl,param,ener,box,neigh,atom,bond,ub,angle,dihe,impr,vdw);
      b=ener->pot;
      
      z[i]=coord-h;
      //energy(ctrl,param,ener,box,neigh,atom,bond,ub,angle,dihe,impr,vdw);
      c=ener->pot;
      
      z[i]=coord+h;
      //energy(ctrl,param,ener,box,neigh,atom,bond,ub,angle,dihe,impr,vdw);
      d=ener->pot;
      
      z[i]=coord+(2.*h);
      //energy(ctrl,param,ener,box,neigh,atom,bond,ub,angle,dihe,impr,vdw);
      e=ener->pot;
      
      z[i]=coord+(3.*h);
      //energy(ctrl,param,ener,box,neigh,atom,bond,ub,angle,dihe,impr,vdw);
      f=ener->pot;
      
      nForce[i].z=-(45.*(d-c)-9.*(e-b)+(f-a))/(60.*h);
      
      z[i]=coord;
    }
  }
  
}
