#include "global.h"
#include "utils.h"

void vdw_full(ATOM *atom,FORCEFIELD *ff,ENERGY *ener,SIMULPARAMS *simulCond)
{
    
  int i,j;
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
      
      r=distance(i,j,atom,delta,simulCond);
      
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
    atom[i].fx+=fxi;
    atom[i].fy+=fyi;
    atom[i].fz+=fzi;
  }
  ener->vdw+=vdw;
}

void vdw_switch(ATOM *atom,FORCEFIELD *ff,ENERGY *ener,SIMULPARAMS *simulCond)
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
  double r,fx,fy,fz,fxi,fyi,fzi;
  double delta[3];
  
  #ifndef _OPENMP
  int ipr=0;
  #else
  #pragma omp parallel default(none) shared(atom,ff,ener,simulCond) private(i,j,k,pvdw,dpvdw,dvdw,switchFunc,dswitchFunc,r,fx,fy,fz,fxi,fyi,fzi,delta) reduction(+:vdw)
  {
    #pragma omp for schedule(dynamic) nowait
    #endif
    for(i=0;i<simulCond->natom-1;i++)
    {
      fxi=0.;
      fyi=0.;
      fzi=0.;
      
      for(k=0;k<ff->verPair[i];k++)
      {
	#ifndef _OPENMP
	j=ff->verList[ipr];
	ipr++;
	#else
	j=ff->verList[ ff->verCumSum[i] + k ];
	#endif
	
	r=distance(i,j,atom,delta,simulCond);
	
	if(r<=simulCond->cuton)
	{
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
	
	  #pragma omp atomic
	  atom[j].fx+=-fx;
	  #pragma omp atomic
	  atom[j].fy+=-fy;
	  #pragma omp atomic
	  atom[j].fz+=-fz;
		  
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
	  
	  fx=dvdw*delta[0]/r;
	  fy=dvdw*delta[1]/r;
	  fz=dvdw*delta[2]/r;
	
	  fxi+=fx;
	  fyi+=fy;
	  fzi+=fz;
	
	  #pragma omp atomic
	  atom[j].fx+=-fx;
	  #pragma omp atomic
	  atom[j].fy+=-fy;
	  #pragma omp atomic
	  atom[j].fz+=-fz;
	  
	}
      }
      
      #pragma omp atomic
      atom[i].fx+=fxi;
      #pragma omp atomic
      atom[i].fy+=fyi;
      #pragma omp atomic
      atom[i].fz+=fzi;

    }
  #ifdef _OPENMP
  }
  #endif
  ener->vdw+=vdw;
  
}

void vdw14_full(ATOM *atom,FORCEFIELD *ff,ENERGY *ener,SIMULPARAMS *simulCond)
{
    
  int i,j,k;
  double vdw=0.,pvdw,dvdw;
  double r,fx,fy,fz;
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
        
    fx=dvdw*delta[0]/r;
    fy=dvdw*delta[1]/r;
    fz=dvdw*delta[2]/r;
    
    atom[i].fx+=fx;
    atom[i].fy+=fy;
    atom[i].fz+=fz;
    
    atom[j].fx+=-fx;
    atom[j].fy+=-fy;
    atom[j].fz+=-fz;
  }
  ener->vdw+=vdw;
}

void vdw14_switch(ATOM *atom,FORCEFIELD *ff,ENERGY *ener,SIMULPARAMS *simulCond)
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
  double r,fx,fy,fz;
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
      
      fx=dvdw*delta[0]/r;
      fy=dvdw*delta[1]/r;
      fz=dvdw*delta[2]/r;
    
      atom[i].fx+=fx;
      atom[i].fy+=fy;
      atom[i].fz+=fz;
    
      atom[j].fx+=-fx;
      atom[j].fy+=-fy;
      atom[j].fz+=-fz;
      	
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
      
      fx=dvdw*delta[0]/r;
      fy=dvdw*delta[1]/r;
      fz=dvdw*delta[2]/r;
    
      atom[i].fx+=fx;
      atom[i].fy+=fy;
      atom[i].fz+=fz;
    
      atom[j].fx+=-fx;
      atom[j].fy+=-fy;
      atom[j].fz+=-fz;
      
    }     
  }
  ener->vdw+=vdw;
}
