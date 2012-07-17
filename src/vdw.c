#include "global.h"
#include "utils.h"

void vdw_full(ATOM *atom,FORCEFIELD *ff,ENERGYFORCE *enerFor,SIMULPARAMS *simulCond)
{
    
  int i,j;
  double vdw=0.,pvdw,dvdw;
  double r,fxi,fyi,fzi;
  double delta[3];
  
  for(i=0;i<atom->natom-1;i++)
  {
    fxi=0.;
    fyi=0.;
    fzi=0.;
    
    for(j=i+1;j<atom->natom;j++)
    {
      
      r=distance(i,j,atom,delta,simulCond);
      
      pvdw=4.*ff->parmVdw[i][0]*ff->parmVdw[j][0]*(X12((ff->parmVdw[i][1]+ff->parmVdw[j][1])/r)-
	X6((ff->parmVdw[i][1]+ff->parmVdw[j][1])/r));
      
      dvdw=24.*ff->parmVdw[i][0]*ff->parmVdw[j][0]/r*(X6((ff->parmVdw[i][1]+ff->parmVdw[j][1])/r)-
	2.*X12((ff->parmVdw[i][1]+ff->parmVdw[j][1])/r));
      
      vdw+=pvdw;
      
      fxi+=dvdw*delta[0]/r;
      fyi+=dvdw*delta[1]/r;
      fzi+=dvdw*delta[2]/r;
      
      atom->fx[j]+=-dvdw*delta[0]/r;
      atom->fy[j]+=-dvdw*delta[1]/r;
      atom->fz[j]+=-dvdw*delta[2]/r;
      
    }
    atom->fx[i]+=fxi;
    atom->fy[i]+=fyi;
    atom->fz[i]+=fzi;
  }
  enerFor->energyVdw+=vdw;
}

void vdw_switch(ATOM *atom,FORCEFIELD *ff,ENERGYFORCE *enerFor,SIMULPARAMS *simulCond)
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
  
  int i,j,k,ipr;
  double vdw=0.,pvdw,dpvdw,dvdw,switchFunc,dswitchFunc;
  double r,fxi,fyi,fzi;
  double delta[3];
  
  ipr=0;
  for(i=0;i<atom->natom-1;i++)
  {
    fxi=0.;
    fyi=0.;
    fzi=0.;
    
    for(k=0;k<ff->verPair[i];k++)
    {
      j=ff->verList[ipr];
      ipr++;
      
      r=distance(i,j,atom,delta,simulCond);
      
      if(r<=simulCond->cuton)
      {
	pvdw=4.*ff->parmVdw[i][0]*ff->parmVdw[j][0]*(X12((ff->parmVdw[i][1]+ff->parmVdw[j][1])/r)-
	  X6((ff->parmVdw[i][1]+ff->parmVdw[j][1])/r));
	
	dvdw=24.*ff->parmVdw[i][0]*ff->parmVdw[j][0]/r*(X6((ff->parmVdw[i][1]+ff->parmVdw[j][1])/r)-
	  2.*X12((ff->parmVdw[i][1]+ff->parmVdw[j][1])/r));
	
	vdw+=pvdw;
	
	fxi+=dvdw*delta[0]/r;
	fyi+=dvdw*delta[1]/r;
	fzi+=dvdw*delta[2]/r;
	
	atom->fx[j]+=-dvdw*delta[0]/r;
	atom->fy[j]+=-dvdw*delta[1]/r;
	atom->fz[j]+=-dvdw*delta[2]/r;	
      }
      else if(r<=simulCond->cutoff)
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
	
	vdw+=pvdw*switchFunc;
	
	dvdw=pvdw*dswitchFunc+dpvdw*switchFunc;
	
	fxi+=dvdw*delta[0]/r;
	fyi+=dvdw*delta[1]/r;
	fzi+=dvdw*delta[2]/r;
	
	atom->fx[j]+=-dvdw*delta[0]/r;
	atom->fy[j]+=-dvdw*delta[1]/r;
	atom->fz[j]+=-dvdw*delta[2]/r;
      }
    }
    atom->fx[i]+=fxi;
    atom->fy[i]+=fyi;
    atom->fz[i]+=fzi;
  }
  enerFor->energyVdw+=vdw;
}

void vdw14_full(ATOM *atom,FORCEFIELD *ff,ENERGYFORCE *enerFor,SIMULPARAMS *simulCond)
{
    
  int i,j,k;
  double vdw=0.,pvdw,dvdw;
  double r;
  double delta[3];
  
  for(k=0;k<ff->npr14;k++)
  {
    i=ff->ver14[k][0];
    j=ff->ver14[k][1];
    
    r=distance(i,j,atom,delta,simulCond);
    
    pvdw=4.*ff->scal14*ff->parmVdw[i][3]*ff->parmVdw[j][3]*(X12((ff->parmVdw[i][4]+ff->parmVdw[j][4])/r)-
      X6((ff->parmVdw[i][4]+ff->parmVdw[j][4])/r));
    
    dvdw=24.*ff->scal14*ff->parmVdw[i][3]*ff->parmVdw[j][3]/r*(X6((ff->parmVdw[i][4]+ff->parmVdw[j][4])/r)-
      2.*X12((ff->parmVdw[i][4]+ff->parmVdw[j][4])/r));
    
    vdw+=pvdw;
    
    atom->fx[i]+=dvdw*delta[0]/r;
    atom->fy[i]+=dvdw*delta[1]/r;
    atom->fz[i]+=dvdw*delta[2]/r;
    
    atom->fx[j]+=-dvdw*delta[0]/r;
    atom->fy[j]+=-dvdw*delta[1]/r;
    atom->fz[j]+=-dvdw*delta[2]/r;
  }
  enerFor->energyVdw+=vdw;
}

void vdw14_switch(ATOM *atom,FORCEFIELD *ff,ENERGYFORCE *enerFor,SIMULPARAMS *simulCond)
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
  
  int i,j,k;
  double vdw=0.,pvdw,dpvdw,dvdw,switchFunc,dswitchFunc;
  double r;
  double delta[3];
  
  for(k=0;k<ff->npr14;k++)
  {
    i=ff->ver14[k][0];
    j=ff->ver14[k][1];
    
    r=distance(i,j,atom,delta,simulCond);
    
    if(r<=simulCond->cuton)
    {
      pvdw=4.*ff->scal14*ff->parmVdw[i][3]*ff->parmVdw[j][3]*(X12((ff->parmVdw[i][4]+ff->parmVdw[j][4])/r)-
	X6((ff->parmVdw[i][4]+ff->parmVdw[j][4])/r));
      
      dvdw=24.*ff->scal14*ff->parmVdw[i][3]*ff->parmVdw[j][3]/r*(X6((ff->parmVdw[i][4]+ff->parmVdw[j][4])/r)-
	2.*X12((ff->parmVdw[i][4]+ff->parmVdw[j][4])/r));
      
      vdw+=pvdw;
      
      atom->fx[i]+=dvdw*delta[0]/r;
      atom->fy[i]+=dvdw*delta[1]/r;
      atom->fz[i]+=dvdw*delta[2]/r;
      
      atom->fx[j]+=-dvdw*delta[0]/r;
      atom->fy[j]+=-dvdw*delta[1]/r;
      atom->fz[j]+=-dvdw*delta[2]/r;	
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
      
      vdw+=pvdw*switchFunc;
      
      dvdw=pvdw*dswitchFunc+dpvdw*switchFunc;
      
      atom->fx[i]+=dvdw*delta[0]/r;
      atom->fy[i]+=dvdw*delta[1]/r;
      atom->fz[i]+=dvdw*delta[2]/r;
      
      atom->fx[j]+=-dvdw*delta[0]/r;
      atom->fy[j]+=-dvdw*delta[1]/r;
      atom->fz[j]+=-dvdw*delta[2]/r;
    }     
  }
  enerFor->energyVdw+=vdw;
}
