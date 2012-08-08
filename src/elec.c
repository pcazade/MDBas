#include "global.h"
#include "utils.h"

void coulomb_full(ATOM *atom,FORCEFIELD *ff,ENERGY *ener,SIMULPARAMS *simulCond,PBC *box)
{
  
  int i,j;
  double elec=0.,pelec,delec;
  double r,fx,fy,fz,fxi,fyi,fzi;
  double delta[3];
  
  for(i=0;i<simulCond->natom-1;i++)
  {
    fxi=0.;
    fyi=0.;
    fzi=0.;
    
    for(j=i+1;j<simulCond->natom;j++)
    {
      
      r=distance(i,j,atom,delta,simulCond,box);
      
      pelec=simulCond->chargeConst*atom[i].q*atom[j].q/r;
      elec+=pelec;
      delec=-pelec/r;
      
      fx=delec*delta[0]/r;
      fy=delec*delta[1]/r;
      fz=delec*delta[2]/r;
      
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
  ener->elec+=elec;
  
}

void coulomb_shift1(ATOM *atom,FORCEFIELD *ff,ENERGY *ener,SIMULPARAMS *simulCond,PBC *box)
{
  
 /*****************************************************************************
 * Shifted electrostatic potential with the shift fuctional form 1:
 * elecShift=elecPot(r)*shiftFunc(r)
 * elecPot=cte*qi*qj/r
 * shiftFunc=1-2r/rc+r**2/rc**2
 * delecShift=delecPot(r)*shiftFunc(r)+elecPot(r)*dshiftFunc(r)
 * delecPot(r)=-elecPot(r)/r
 * dshiftFunc(r)=-2/rc+2r/rc**2
 ****************************************************************************/
  
  int i,j,k;
  double elec=0.,pelec,delec,shiftFunc,dshiftFunc;
  double r,fx,fy,fz,fxi,fyi,fzi;
  double delta[3];
  
  #ifndef _OPENMP
  int ipr=0;
  #else
  #pragma omp parallel default(none) shared(atom,ff,ener,simulCond) private(i,j,k,pelec,delec,shiftFunc,dshiftFunc,r,fx,fy,fz,fxi,fyi,fzi,delta) reduction(+:elec)
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
	
	r=distance(i,j,atom,delta,simulCond,box);
	
	if(r<=simulCond->cutoff)
	{
	  shiftFunc=1.-2.*r/simulCond->cutoff+X2(r/simulCond->cutoff);
	  dshiftFunc=-2./simulCond->cutoff+2.*r/X2(simulCond->cutoff);
	    
	  pelec=simulCond->chargeConst*atom[i].q*atom[j].q/r;
	  elec+=pelec*shiftFunc;
	  delec=pelec*(dshiftFunc-shiftFunc/r);
	  
	  fx=delec*delta[0]/r;
	  fy=delec*delta[1]/r;
	  fz=delec*delta[2]/r;
	  
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
  
  ener->elec+=elec;
  
}

void coulomb_shift2(ATOM *atom,FORCEFIELD *ff,ENERGY *ener,SIMULPARAMS *simulCond,PBC *box)
{
  
/*****************************************************************************
 * Shifted electrostatic potential with the shift fuctional form 2 :
 * elecShift=elecPot(r)*shiftFunc(r)
 * elecPot=cte*qi*qj/r
 * shiftFunc=1-2r**2/rc**2+r**4/rc**4
 * delecShift=delecPot(r)*shiftFunc(r)+elecPot(r)*dshiftFunc(r)
 * delecPot(r)=-elecPot(r)/r
 * dshiftFunc(r)=-4r/rc**2+4r**3/rc**4
 ****************************************************************************/
  
  int i,j,k;
  double elec=0.,pelec,delec,shiftFunc,dshiftFunc;
  double r,fx,fy,fz,fxi,fyi,fzi;
  double delta[3];
  

  #ifndef _OPENMP
  int ipr=0;
  #else
  #pragma omp parallel default(none) shared(atom,ff,ener,simulCond) private(i,j,k,pelec,delec,shiftFunc,dshiftFunc,r,fx,fy,fz,fxi,fyi,fzi,delta) reduction(+:elec)
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
	
	r=distance(i,j,atom,delta,simulCond,box);
	
	if(r<=simulCond->cutoff)
	{
	  shiftFunc=1.-2.*X2(r/simulCond->cutoff)+X4(r/simulCond->cutoff);
	  dshiftFunc=-4.*r/X2(simulCond->cutoff)+4.*X3(r)/X4(simulCond->cutoff);
	  
	  pelec=simulCond->chargeConst*atom[i].q*atom[j].q/r;
	  elec+=pelec*shiftFunc;
	  delec=pelec*(dshiftFunc-shiftFunc/r);
	  
	  fx=delec*delta[0]/r;
	  fy=delec*delta[1]/r;
	  fz=delec*delta[2]/r;
	  
	  fxi+=fx;
	  fyi+=fy;
	  fzi+=fz;
	  
	  #pragma omp atomic
	  atom[j].fx+=-fx;
	  #pragma omp atomic
	  atom[j].fy+=-fy;
	  #pragma omp atomic
	  atom[j].fz+=-fz;

	  
	}	//END of if
      }	//END of inner for
      
      #pragma omp atomic
      atom[i].fx+=fxi;
      #pragma omp atomic
      atom[i].fy+=fyi;
      #pragma omp atomic
      atom[i].fz+=fzi;
      
    } //END of main for
  
  #ifdef _OPENMP
  }
  #endif
  
  ener->elec+=elec;
  
} //END of function

void coulomb_switch(ATOM *atom,FORCEFIELD *ff,ENERGY *ener,SIMULPARAMS *simulCond,PBC *box)
{
  
/*****************************************************************************
 * Switched electrostatic potential :
 * elecSwitch=elecPot(r)*switchFunc(r)
 * elecPot=cte*qi*qj/r
 * switchFunc=(rc**2+2r**2-3ro**2)*(rc**2-r**2)**2/(rc**2-ro**2)**3
 * delecSwitch=delecPot(r)*switchFunc(r)+elecPot(r)*dswitchFunc(r)
 * delecPot(r)=-elecPot(r)/r
 * dswitchFunc(r)=-12*r*(rc**2-r**2)*(ro**2-r**2)/(rc**2-ro**2)**3
 ****************************************************************************/
  
  int i,j,k;
  double elec=0.,pelec,delec,switchFunc,dswitchFunc;
  double r,fx,fy,fz,fxi,fyi,fzi;
  double delta[3];
  
  #ifndef _OPENMP
  int ipr=0;
  #else
  #pragma omp parallel default(none) shared(atom,ff,ener,simulCond) private(i,j,k,pelec,delec,switchFunc,dswitchFunc,r,fx,fy,fz,fxi,fyi,fzi,delta) reduction(+:elec)
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
	
	r=distance(i,j,atom,delta,simulCond,box);
	
	if(r<=simulCond->cuton)
	{
	  
	  pelec=simulCond->chargeConst*atom[i].q*atom[j].q/r;
	  elec+=pelec;
	  delec=-pelec/r;
	  
	  fx=delec*delta[0]/r;
	  fy=delec*delta[1]/r;
	  fz=delec*delta[2]/r;
	  
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
	  
	  pelec=simulCond->chargeConst*atom[i].q*atom[j].q/r;
	  elec+=pelec*switchFunc;
	  delec=pelec*(dswitchFunc-switchFunc/r);
	  
	  fx=delec*delta[0]/r;
	  fy=delec*delta[1]/r;
	  fz=delec*delta[2]/r;
	  
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
  ener->elec+=elec;
}

void coulomb14_full(ATOM *atom,FORCEFIELD *ff,ENERGY *ener,SIMULPARAMS *simulCond,PBC *box)
{
  
  int i,j,k;
  double elec=0.,pelec,delec;
  double r,fx,fy,fz;
  double delta[3];
  
  for(k=0;k<ff->npr14;k++)
  {
    i=ff->ver14[k][0];
    j=ff->ver14[k][1];
    
    r=distance(i,j,atom,delta,simulCond,box);
    pelec=ff->scal14*simulCond->chargeConst*atom[i].q*atom[j].q/r;
    elec+=pelec;
    delec=-pelec/r;
    
    fx=delec*delta[0]/r;
    fy=delec*delta[1]/r;
    fz=delec*delta[2]/r;
    
    atom[i].fx+=fx;
    atom[i].fy+=fy;
    atom[i].fz+=fz;
    
    atom[j].fx+=-fx;
    atom[j].fy+=-fy;
    atom[j].fz+=-fz;
    
  }
  ener->elec+=elec;
}

void coulomb14_shift1(ATOM *atom,FORCEFIELD *ff,ENERGY *ener,SIMULPARAMS *simulCond,PBC *box)
{
  
 /*****************************************************************************
 * Shifted electrostatic potential with the shift fuctional form 1:
 * elecShift=elecPot(r)*shiftFunc(r)
 * elecPot=cte*qi*qj/r
 * shiftFunc=1-2r/rc+r**2/rc**2
 * delecShift=delecPot(r)*shiftFunc(r)+elecPot(r)*dshiftFunc(r)
 * delecPot(r)=-elecPot(r)/r
 * dshiftFunc(r)=-2/rc+2r/rc**2
 ****************************************************************************/
  
  int i,j,k;
  double elec=0.,pelec,delec,shiftFunc,dshiftFunc;
  double r,fx,fy,fz;
  double delta[3];
  
  for(k=0;k<ff->npr14;k++)
  {
    i=ff->ver14[k][0];
    j=ff->ver14[k][1];
    
    r=distance(i,j,atom,delta,simulCond,box);
      
    if(r<=simulCond->cutoff)
    {
      shiftFunc=1.-2.*r/simulCond->cutoff+X2(r/simulCond->cutoff);
      dshiftFunc=-2./simulCond->cutoff+2.*r/X2(simulCond->cutoff);
      
      pelec=ff->scal14*simulCond->chargeConst*atom[i].q*atom[j].q/r;
      elec+=pelec*shiftFunc;
      delec=pelec*(dshiftFunc-shiftFunc/r);
      
      fx=delec*delta[0]/r;
      fy=delec*delta[1]/r;
      fz=delec*delta[2]/r;
    
      atom[i].fx+=fx;
      atom[i].fy+=fy;
      atom[i].fz+=fz;
    
      atom[j].fx+=-fx;
      atom[j].fy+=-fy;
      atom[j].fz+=-fz;
      
    }
  }
  ener->elec+=elec;
}

void coulomb14_shift2(ATOM *atom,FORCEFIELD *ff,ENERGY *ener,SIMULPARAMS *simulCond,PBC *box)
{
  
/*****************************************************************************
 * Shifted electrostatic potential with the shift fuctional form 2 :
 * elecShift=elecPot(r)*shiftFunc(r)
 * elecPot=cte*qi*qj/r
 * shiftFunc=1-2r**2/rc**2+r**4/rc**4
 * delecShift=delecPot(r)*shiftFunc(r)+elecPot(r)*dshiftFunc(r)
 * delecPot(r)=-elecPot(r)/r
 * dshiftFunc(r)=-4r/rc**2+4r**3/rc**4
 ****************************************************************************/
  
  int i,j,k;
  double elec=0.,pelec,delec,shiftFunc,dshiftFunc;
  double r,fx,fy,fz;
  double delta[3];
  
  for(k=0;k<ff->npr14;k++)
  {
    i=ff->ver14[k][0];
    j=ff->ver14[k][1];
    
    r=distance(i,j,atom,delta,simulCond,box);
    
    if(r<=simulCond->cutoff)
    {
      shiftFunc=1.-2.*X2(r/simulCond->cutoff)+X4(r/simulCond->cutoff);
      dshiftFunc=-4.*r/X2(simulCond->cutoff)+4.*X3(r)/X4(simulCond->cutoff);
      
      pelec=ff->scal14*simulCond->chargeConst*atom[i].q*atom[j].q/r;
      elec+=pelec*shiftFunc;
      delec=pelec*(dshiftFunc-shiftFunc/r);
      
      fx=delec*delta[0]/r;
      fy=delec*delta[1]/r;
      fz=delec*delta[2]/r;
    
      atom[i].fx+=fx;
      atom[i].fy+=fy;
      atom[i].fz+=fz;
    
      atom[j].fx+=-fx;
      atom[j].fy+=-fy;
      atom[j].fz+=-fz;
      
    }     
  }
    
  ener->elec+=elec;
}

void coulomb14_switch(ATOM *atom,FORCEFIELD *ff,ENERGY *ener,SIMULPARAMS *simulCond,PBC *box)
{
  
/*****************************************************************************
 * Switched electrostatic potential :
 * elecSwitch=elecPot(r)*switchFunc(r)
 * elecPot=cte*qi*qj/r
 * switchFunc=(rc**2+2r**2-3ro**2)*(rc**2-r**2)**2/(rc**2-ro**2)**3
 * delecSwitch=delecPot(r)*switchFunc(r)+elecPot(r)*dswitchFunc(r)
 * delecPot(r)=-elecPot(r)/r
 * dswitchFunc(r)=-12*r*(rc**2-r**2)*(ro**2-r**2)/(rc**2-ro**2)**3
 ****************************************************************************/
  
  int i,j,k;
  double elec=0.,pelec,delec,switchFunc,dswitchFunc;
  double r,fx,fy,fz;
  double delta[3];
  
  for(k=0;k<ff->npr14;k++)
  {
    i=ff->ver14[k][0];
    j=ff->ver14[k][1];
    
    r=distance(i,j,atom,delta,simulCond,box);
    
    if(r<=simulCond->cuton)
    {
	pelec=ff->scal14*simulCond->chargeConst*atom[i].q*atom[j].q/r;
	elec+=pelec;
	delec=-pelec/r;
	
	fx=delec*delta[0]/r;
	fy=delec*delta[1]/r;
	fz=delec*delta[2]/r;
    
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
	
	pelec=ff->scal14*simulCond->chargeConst*atom[i].q*atom[j].q/r;
	elec+=pelec*switchFunc;
	delec=pelec*(dswitchFunc-switchFunc/r);
	
	fx=delec*delta[0]/r;
	fy=delec*delta[1]/r;
	fz=delec*delta[2]/r;
    
	atom[i].fx+=fx;
	atom[i].fy+=fy;
	atom[i].fz+=fz;
    
	atom[j].fx+=-fx;
	atom[j].fy+=-fy;
	atom[j].fz+=-fz;
	
    }     
  }
  
  ener->elec+=elec;
}
