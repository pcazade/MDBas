#include "global.h"
#include "utils.h"

void coulomb_full(ATOM *atom,FORCEFIELD *ff,ENERGYFORCE *enerFor,SIMULPARAMS *simulCond)
{
  
  int i,j;
  double elec=0.,pelec,delec;
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
      
      pelec=simulCond->chargeConst*ff->q[i]*ff->q[j]/r;
      elec+=pelec;
      delec=-pelec/r;
      
      fxi+=delec*delta[0]/r;
      fyi+=delec*delta[1]/r;
      fzi+=delec*delta[2]/r;
      
      atom->fx[j]+=-delec*delta[0]/r;
      atom->fy[j]+=-delec*delta[1]/r;
      atom->fz[j]+=-delec*delta[2]/r;
      
    }
    
    atom->fx[i]+=fxi;
    atom->fy[i]+=fyi;
    atom->fz[i]+=fzi;
  }
  enerFor->energyElec+=elec;
  
}

void coulomb_shift1(ATOM *atom,FORCEFIELD *ff,ENERGYFORCE *enerFor,SIMULPARAMS *simulCond)
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
  
  int i,j,k,ipr;
  double elec=0.,pelec,delec,shiftFunc,dshiftFunc;
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
      
      if(r<=simulCond->cutoff)
      {
	shiftFunc=1.-2.*r/simulCond->cutoff+X2(r/simulCond->cutoff);
	dshiftFunc=-2./simulCond->cutoff+2.*r/X2(simulCond->cutoff);
	  
	pelec=simulCond->chargeConst*ff->q[i]*ff->q[j]/r;
	elec+=pelec*shiftFunc;
	delec=pelec*(dshiftFunc-shiftFunc/r);
	
	fxi+=delec*delta[0]/r;
	fyi+=delec*delta[1]/r;
	fzi+=delec*delta[2]/r;
	
	atom->fx[j]+=-delec*delta[0]/r;
	atom->fy[j]+=-delec*delta[1]/r;
	atom->fz[j]+=-delec*delta[2]/r;
      } 
    }
    atom->fx[i]+=fxi;
    atom->fy[i]+=fyi;
    atom->fz[i]+=fzi;
  }
  enerFor->energyElec+=elec;
}

void coulomb_shift2(ATOM *atom,FORCEFIELD *ff,ENERGYFORCE *enerFor,SIMULPARAMS *simulCond)
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
  
  int i,j,k,ipr;
  double elec=0.,pelec,delec,shiftFunc,dshiftFunc;
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
      
      if(r<=simulCond->cutoff)
      {
	shiftFunc=1.-2.*X2(r/simulCond->cutoff)+X4(r/simulCond->cutoff);
	dshiftFunc=-4.*r/X2(simulCond->cutoff)+4.*X3(r)/X4(simulCond->cutoff);
	
	pelec=simulCond->chargeConst*ff->q[i]*ff->q[j]/r;
	elec+=pelec*shiftFunc;
	delec=pelec*(dshiftFunc-shiftFunc/r);
	
	fxi+=delec*delta[0]/r;
	fyi+=delec*delta[1]/r;
	fzi+=delec*delta[2]/r;
	
	atom->fx[j]+=-delec*delta[0]/r;
	atom->fy[j]+=-delec*delta[1]/r;
	atom->fz[j]+=-delec*delta[2]/r;
      }     
    }
    atom->fx[i]+=fxi;
    atom->fy[i]+=fyi;
    atom->fz[i]+=fzi;
  }
  enerFor->energyElec+=elec;
  
}

void coulomb_switch(ATOM *atom,FORCEFIELD *ff,ENERGYFORCE *enerFor,SIMULPARAMS *simulCond)
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
  
  int i,j,k,ipr;
  double elec=0.,pelec,delec,switchFunc,dswitchFunc;
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
	
	pelec=simulCond->chargeConst*ff->q[i]*ff->q[j]/r;
	elec+=pelec;
	delec=-pelec/r;
	
	fxi+=delec*delta[0]/r;
	fyi+=delec*delta[1]/r;
	fzi+=delec*delta[2]/r;
	
	atom->fx[j]+=-delec*delta[0]/r;
	atom->fy[j]+=-delec*delta[1]/r;
	atom->fz[j]+=-delec*delta[2]/r;
      }
      else if(r<=simulCond->cutoff)
      {
	switchFunc=X2(X2(simulCond->cutoff)-X2(r))*
	  (X2(simulCond->cutoff)+2.*X2(r)-3.*X2(simulCond->cuton))/
	  X3(X2(simulCond->cutoff)-X2(simulCond->cuton));
	  
	dswitchFunc=12.*r*(X2(simulCond->cutoff)-X2(r))*
	  (X2(simulCond->cuton)-X2(r))/
	  X3(X2(simulCond->cutoff)-X2(simulCond->cuton));
	
	pelec=simulCond->chargeConst*ff->q[i]*ff->q[j]/r;
	elec+=pelec*switchFunc;
	delec=pelec*(dswitchFunc-switchFunc/r);
	
	fxi+=delec*delta[0]/r;
	fyi+=delec*delta[1]/r;
	fzi+=delec*delta[2]/r;
	
	atom->fx[j]+=-delec*delta[0]/r;
	atom->fy[j]+=-delec*delta[1]/r;
	atom->fz[j]+=-delec*delta[2]/r;
      }     
    }
    atom->fx[i]+=fxi;
    atom->fy[i]+=fyi;
    atom->fz[i]+=fzi;
  }
  enerFor->energyElec+=elec;
}

void coulomb14_full(ATOM *atom,FORCEFIELD *ff,ENERGYFORCE *enerFor,SIMULPARAMS *simulCond)
{
  
  int i,j,k;
  double elec=0.,pelec,delec;
  double r;
  double delta[3];
  
  for(k=0;k<ff->npr14;k++)
  {
    i=ff->ver14[k][0];
    j=ff->ver14[k][1];
    
    r=distance(i,j,atom,delta,simulCond);
    pelec=ff->scal14*simulCond->chargeConst*ff->q[i]*ff->q[j]/r;
    elec+=pelec;
    delec=-pelec/r;
    
    atom->fx[i]+=delec*delta[0]/r;
    atom->fy[i]+=delec*delta[1]/r;
    atom->fz[i]+=delec*delta[2]/r;
    
    atom->fx[j]+=-delec*delta[0]/r;
    atom->fy[j]+=-delec*delta[1]/r;
    atom->fz[j]+=-delec*delta[2]/r;
  }
  enerFor->energyElec+=elec;
}

void coulomb14_shift1(ATOM *atom,FORCEFIELD *ff,ENERGYFORCE *enerFor,SIMULPARAMS *simulCond)
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
  double r;
  double delta[3];
  
  for(k=0;k<ff->npr14;k++)
  {
    i=ff->ver14[k][0];
    j=ff->ver14[k][1];
    
    r=distance(i,j,atom,delta,simulCond);
      
    if(r<=simulCond->cutoff)
    {
      shiftFunc=1.-2.*r/simulCond->cutoff+X2(r/simulCond->cutoff);
      dshiftFunc=-2./simulCond->cutoff+2.*r/X2(simulCond->cutoff);
      
      pelec=ff->scal14*simulCond->chargeConst*ff->q[i]*ff->q[j]/r;
      elec+=pelec*shiftFunc;
      delec=pelec*(dshiftFunc-shiftFunc/r);
      
      atom->fx[i]+=delec*delta[0]/r;
      atom->fy[i]+=delec*delta[1]/r;
      atom->fz[i]+=delec*delta[2]/r;
      
      atom->fx[j]+=-delec*delta[0]/r;
      atom->fy[j]+=-delec*delta[1]/r;
      atom->fz[j]+=-delec*delta[2]/r;
    }
  }
  enerFor->energyElec+=elec;
}

void coulomb14_shift2(ATOM *atom,FORCEFIELD *ff,ENERGYFORCE *enerFor,SIMULPARAMS *simulCond)
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
  double r;
  double delta[3];
  
  for(k=0;k<ff->npr14;k++)
  {
    i=ff->ver14[k][0];
    j=ff->ver14[k][1];
    
    r=distance(i,j,atom,delta,simulCond);
    
    if(r<=simulCond->cutoff)
    {
      shiftFunc=1.-2.*X2(r/simulCond->cutoff)+X4(r/simulCond->cutoff);
      dshiftFunc=-4.*r/X2(simulCond->cutoff)+4.*X3(r)/X4(simulCond->cutoff);
      
      pelec=ff->scal14*simulCond->chargeConst*ff->q[i]*ff->q[j]/r;
      elec+=pelec*shiftFunc;
      delec=pelec*(dshiftFunc-shiftFunc/r);
      
      atom->fx[i]+=delec*delta[0]/r;
      atom->fy[i]+=delec*delta[1]/r;
      atom->fz[i]+=delec*delta[2]/r;
      
      atom->fx[j]+=-delec*delta[0]/r;
      atom->fy[j]+=-delec*delta[1]/r;
      atom->fz[j]+=-delec*delta[2]/r;
      
    }     
  }
    
  enerFor->energyElec+=elec;
}

void coulomb14_switch(ATOM *atom,FORCEFIELD *ff,ENERGYFORCE *enerFor,SIMULPARAMS *simulCond)
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
  double r;
  double delta[3];
  
  for(k=0;k<ff->npr14;k++)
  {
    i=ff->ver14[k][0];
    j=ff->ver14[k][1];
    
    r=distance(i,j,atom,delta,simulCond);
    
    if(r<=simulCond->cuton)
    {
	pelec=ff->scal14*simulCond->chargeConst*ff->q[i]*ff->q[j]/r;
	elec+=pelec;
	delec=-pelec/r;
	
	atom->fx[i]+=delec*delta[0]/r;
	atom->fy[i]+=delec*delta[1]/r;
	atom->fz[i]+=delec*delta[2]/r;
	
	atom->fx[j]+=-delec*delta[0]/r;
	atom->fy[j]+=-delec*delta[1]/r;
	atom->fz[j]+=-delec*delta[2]/r;
    }
    else if(r<=simulCond->cutoff)
    {
	switchFunc=X2(X2(simulCond->cutoff)-X2(r))*
	  (X2(simulCond->cutoff)+2.*X2(r)-3.*X2(simulCond->cuton))/
	  X3(X2(simulCond->cutoff)-X2(simulCond->cuton));
	  
	dswitchFunc=12.*r*(X2(simulCond->cutoff)-X2(r))*
	  (X2(simulCond->cuton)-X2(r))/
	  X3(X2(simulCond->cutoff)-X2(simulCond->cuton));
	
	pelec=ff->scal14*simulCond->chargeConst*ff->q[i]*ff->q[j]/r;
	elec+=pelec*switchFunc;
	delec=pelec*(dswitchFunc-switchFunc/r);
	
	atom->fx[i]+=delec*delta[0]/r;
	atom->fy[i]+=delec*delta[1]/r;
	atom->fz[i]+=delec*delta[2]/r;
	
	atom->fx[j]+=-delec*delta[0]/r;
	atom->fy[j]+=-delec*delta[1]/r;
	atom->fz[j]+=-delec*delta[2]/r;
    }     
  }
  
  enerFor->energyElec+=elec;
}
