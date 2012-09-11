#include "global.h"
#include "utils.h"

double vdw_none(ATOM atom[],FORCEFIELD *ff,SIMULPARAMS *simulCond,PBC *box,
		int i, int j, double r, double *dvdw)
{
  *dvdw=0.;
  return 0.;
}

void vdw_full(ATOM atom[],FORCEFIELD *ff,ENERGY *ener,SIMULPARAMS *simulCond,PBC *box)
{
    
  int i,j,k,exclude;;
  double vdw=0.,pvdw,dvdw;
  double r,fx,fy,fz,fxi,fyi,fzi;
  double delta[3];
  
  for(i=0;i<simulCond->natom-1;i++)
  {
    fxi=0.;
    fyi=0.;
    fzi=0.;
    
    for(j=i+1;j<simulCond->natom;j++)
    {
      
      exclude=0;
      for (k=0;k<simulCond->excludeNum[i];k++)
      {
	if(simulCond->excludeAtom[i][k]==j)
	{
	  exclude=1;
	  break;
	}
      }
      
      if(!exclude)
      {
	
	r=distance(i,j,atom,delta,simulCond,box);
	
	pvdw=4.*ff->parmVdw[i][0]*ff->parmVdw[j][0]*(X12((ff->parmVdw[i][1]+ff->parmVdw[j][1])/r)-
	  X6((ff->parmVdw[i][1]+ff->parmVdw[j][1])/r));
	
	dvdw=24.*ff->parmVdw[i][0]*ff->parmVdw[j][0]/r*(X6((ff->parmVdw[i][1]+ff->parmVdw[j][1])/r)-
	  2.*X12((ff->parmVdw[i][1]+ff->parmVdw[j][1])/r));
	
	vdw+=pvdw;
	
	fx=dvdw*delta[0]/r;
	fy=dvdw*delta[1]/r;
	fz=dvdw*delta[2]/r;
	
	fxi+=fx;
	fyi+=fy;
	fzi+=fz;
	
	atom[j].fx+=-fx;
	atom[j].fy+=-fy;
	atom[j].fz+=-fz;
	
      }
      
    }
    atom[i].fx+=fxi;
    atom[i].fy+=fyi;
    atom[i].fz+=fzi;
  }
  ener->vdw+=vdw;
}

double vdw_switch(ATOM atom[],FORCEFIELD *ff,SIMULPARAMS *simulCond,PBC *box,
		int i, int j, double r, double *dvdw)
{
  
/**
 * \brief Switched van der Waals potential
 * vdwSwitch=elecPot(r)*switchFunc(r)
 * vdwPot=4*eps*((sig/r)**12-(sig/r)**6)
 * switchFunc=(rc**2+2r**2-3ro**2)*(rc**2-r**2)**2/(rc**2-ro**2)**3
 * dvdwSwitch=delecPot(r)*switchFunc(r)+elecPot(r)*dswitchFunc(r)
 * dvdwPot(r)=-elecPot(r)/r
 * dswitchFunc(r)=-12*r*(rc**2-r**2)*(ro**2-r**2)/(rc**2-ro**2)**3
 *
 */
  
  double vdw=0.,pvdw,dpvdw,switchFunc,dswitchFunc;
  
  if(r<=simulCond->cuton)
  {
    
    pvdw=4.*ff->parmVdw[i][0]*ff->parmVdw[j][0]*(X12((ff->parmVdw[i][1]+ff->parmVdw[j][1])/r)-
      X6((ff->parmVdw[i][1]+ff->parmVdw[j][1])/r));
    
    *dvdw=24.*ff->parmVdw[i][0]*ff->parmVdw[j][0]/r*(X6((ff->parmVdw[i][1]+ff->parmVdw[j][1])/r)-
      2.*X12((ff->parmVdw[i][1]+ff->parmVdw[j][1])/r));
    
    vdw=pvdw;
    
  }
  else
  {
    
    switchFunc=X2(X2(simulCond->cutoff)-X2(r))*
      (X2(simulCond->cutoff)+2.*X2(r)-3.*X2(simulCond->cuton))/
      X3(X2(simulCond->cutoff)-X2(simulCond->cuton));
      
    dswitchFunc=12.*r*(X2(simulCond->cutoff)-X2(r))*
      (X2(simulCond->cuton)-X2(r))/
      X3(X2(simulCond->cutoff)-X2(simulCond->cuton));
    
    pvdw=4.*ff->parmVdw[i][0]*ff->parmVdw[j][0]*(X12((ff->parmVdw[i][1]+ff->parmVdw[j][1])/r)-
      X6((ff->parmVdw[i][1]+ff->parmVdw[j][1])/r));
    
    dpvdw=24.*ff->parmVdw[i][0]*ff->parmVdw[j][0]/r*(X6((ff->parmVdw[i][1]+ff->parmVdw[j][1])/r)-
      2.*X12((ff->parmVdw[i][1]+ff->parmVdw[j][1])/r));
    
    vdw=pvdw*switchFunc;
    
    *dvdw=pvdw*dswitchFunc+dpvdw*switchFunc;
    
  }

  return vdw;
}

double vdw14_none(ATOM atom[],FORCEFIELD *ff,SIMULPARAMS *simulCond,PBC *box,
		  int i, int j, double r, double *dvdw)
{
  *dvdw=0;
  return 0.;
}

double vdw14_full(ATOM atom[],FORCEFIELD *ff,SIMULPARAMS *simulCond,PBC *box,
		int i, int j, double r, double *dvdw)
{
  
  double vdw=0.;
  
  vdw=4.*ff->scal14*ff->parmVdw[i][3]*ff->parmVdw[j][3]*(X12((ff->parmVdw[i][4]+ff->parmVdw[j][4])/r)-
    X6((ff->parmVdw[i][4]+ff->parmVdw[j][4])/r));
  
  *dvdw=24.*ff->scal14*ff->parmVdw[i][3]*ff->parmVdw[j][3]/r*(X6((ff->parmVdw[i][4]+ff->parmVdw[j][4])/r)-
    2.*X12((ff->parmVdw[i][4]+ff->parmVdw[j][4])/r));
  
  return vdw;
}

double vdw14_switch(ATOM atom[],FORCEFIELD *ff,SIMULPARAMS *simulCond,PBC *box,
		  int i, int j, double r, double *dvdw)
{
  
/*****************************************************************************
 * Switched van der Waals potential :
 * vdwSwitch=elecPot(r)*switchFunc(r)
 * vdwPot=4*eps*((sig/r)**12-(sig/r)**6)
 * switchFunc=(rc**2+2r**2-3ro**2)*(rc**2-r**2)**2/(rc**2-ro**2)**3
 * dvdwSwitch=delecPot(r)*switchFunc(r)+elecPot(r)*dswitchFunc(r)
 * dvdwPot(r)=-elecPot(r)/r
 * dswitchFunc(r)=-12*r*(rc**2-r**2)*(ro**2-r**2)/(rc**2-ro**2)**3
 ****************************************************************************/
  
  double vdw=0.,pvdw,dpvdw,switchFunc,dswitchFunc;
  
  if(r<=simulCond->cuton)
  {
    
    pvdw=4.*ff->scal14*ff->parmVdw[i][3]*ff->parmVdw[j][3]*(X12((ff->parmVdw[i][4]+ff->parmVdw[j][4])/r)-
      X6((ff->parmVdw[i][4]+ff->parmVdw[j][4])/r));
    
    *dvdw=24.*ff->scal14*ff->parmVdw[i][3]*ff->parmVdw[j][3]/r*(X6((ff->parmVdw[i][4]+ff->parmVdw[j][4])/r)-
      2.*X12((ff->parmVdw[i][4]+ff->parmVdw[j][4])/r));
    
    vdw=pvdw;
    
  }
  else if(r<=simulCond->cutoff)
  {
    
    switchFunc=X2(X2(simulCond->cutoff)-X2(r))*
      (X2(simulCond->cutoff)+2.*X2(r)-3.*X2(simulCond->cuton))/
      X3(X2(simulCond->cutoff)-X2(simulCond->cuton));
      
    dswitchFunc=12.*r*(X2(simulCond->cutoff)-X2(r))*
      (X2(simulCond->cuton)-X2(r))/
      X3(X2(simulCond->cutoff)-X2(simulCond->cuton));
    
    pvdw=4.*ff->scal14*ff->parmVdw[i][3]*ff->parmVdw[j][3]*(X12((ff->parmVdw[i][4]+ff->parmVdw[j][4])/r)-
      X6((ff->parmVdw[i][4]+ff->parmVdw[j][4])/r));
    
    dpvdw=24.*ff->scal14*ff->parmVdw[i][3]*ff->parmVdw[j][3]/r*(X6((ff->parmVdw[i][4]+ff->parmVdw[j][4])/r)-
      2.*X12((ff->parmVdw[i][4]+ff->parmVdw[j][4])/r));
    
    vdw=pvdw*switchFunc;
    
    *dvdw=pvdw*dswitchFunc+dpvdw*switchFunc;
    
  }
  
  return vdw;
}
