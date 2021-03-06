/**
 * \file elec.c
 * \brief Contains functions for evaluating electrostatic energies and forces.
 * \author Pierre-Andre Cazade and Florent Hedin
 * \version alpha-branch
 * \date 2012
 */

#include <math.h>

#include "global.h"
#include "utils.h"

/**
 * \brief Empty function called when electrostatic energy and force are disabled.
 * 
 * \remarks delec is set to 0.0 .
 *
 * \return On return 0.0 as no energy evaluated.
 */
double coulomb_none(const PARAM *param,double *delec,const double qel,
		    const double r2,const double rt)
{
  *delec=0.;
  return 0.;
}

/**
 * \param atom Array of structure ATOM (coordinates, forces, etc...).
 * \param ff Pointer to structure FORCEFIELD containing forcefield parameters.
 * \param ener Pointer to structure ENERGY containing values of the different energies.
 * \param simulCond Pointer to structure SIMULPARAMS containing parameters of the current simulation.
 * \param box Pointer to structure PBC containing Periodic Boundaries Conditions parameters.
 *
 * \brief Function called for a full evaluation of the electrostatic energy and force.
 */
void coulomb_full(ENERGY *ener,PARAM *param,PBC *box,double *x,double *y,
		  double *z,double *fx, double *fy, double *fz,double *q,
		  int **exclList,int *exclPair)
{
  
  int i,j,k,exclude;
  double elec=0.,pelec,delec;
  double r,r2,rt,fxi,fyi,fzi,fxj,fyj,fzj;
  double delta[3];
  
  for(i=0;i<param->nAtom-1;i++)
  {
    fxi=0.;
    fyi=0.;
    fzi=0.;
    
    for(j=i+1;j<param->nAtom;j++)
    {
      
      exclude=0;
      for (k=0;k<exclPair[i];k++)
      {
	if(exclList[i][k]==j)
	{
	  exclude=1;
	  break;
	}
      }
      
      if(!exclude)
      {
	
	delta[0]=x[j]-x[i];
	delta[1]=y[j]-y[i];
	delta[2]=z[j]-z[i];
	
	r2=dist(box,delta);
	r=sqrt(r2);
	rt=1./r;
	
	pelec=param->chargeConst*q[i]*q[j]*rt;
	elec+=pelec;
	delec=-pelec*rt;
	
	fxj=delec*delta[0]*rt;
	fyj=delec*delta[1]*rt;
	fzj=delec*delta[2]*rt;
	
	fxi+=fxj;
	fyi+=fyj;
	fzi+=fzj;
	
	fx[j]+=-fxj;
	fy[j]+=-fyj;
	fz[j]+=-fzj;
	
      }
      
    }
    
    fx[i]+=fxi;
    fy[i]+=fyi;
    fz[i]+=fzi;
  }
  ener->elec+=elec;
  
}

/**
 * \param atom Array of structure ATOM (coordinates, forces, etc...).
 * \param ff Pointer to structure FORCEFIELD containing forcefield parameters.
 * \param simulCond Pointer to structure SIMULPARAMS containing parameters of the current simulation.
 * \param box Pointer to structure PBC containing Periodic Boundaries Conditions parameters.
 * \param i Index of first atom.
 * \param j Index of second atom.
 * \param r Distance between the two atoms.
 * \param delec Pointer to derivative of energy used for force evaluation.
 * 
 * \brief Evaluates the electrostatic energy and force for a pair when using the SHIFT_1 cutoff.
 *
 * \return On return the electrostatic energy.
 *
 *
 * Shifted electrostatic potential with the shift functional form 1:
 *
 * \f$ elecShift=elecPot(r)*shiftFunc(r) \f$
 *
 * \f$ elecPot=cte*qi*qj/r \f$
 * 
 * Where cte is chgcharmm or chgnamd
 *
 * \f$ shiftFunc=1-2r/rc+r^2/rc^2 \f$
 *
 * \f$ delecShift=delecPot(r)*shiftFunc(r)+elecPot(r)*dshiftFunc(r) \f$
 *
 * \f$ delecPot(r)=-elecPot(r)/r \f$
 * 
 * \f$ dshiftFunc(r)=-2/rc+2r/rc^2 \f$
 *
 */
double coulomb_shift1(const PARAM *param,double *delec,const double qel,
		      const double r2,const double rt)
{ 
  double elec=0.,pelec,shift1,shift2,shiftFunc,dshiftFunc;
  
  shift1=rt*r2*param->rcutOff;
  shift2=r2*param->rcutOff2;

  shiftFunc=1.-2.*shift1+shift2;
  dshiftFunc=2.*(shift2-shift1);
    
  pelec=param->chargeConst*qel*rt;
  elec=pelec*shiftFunc;
  *delec=pelec*rt*(dshiftFunc-shiftFunc);
  
  return elec;
}

/**
 * \param atom Array of structure ATOM (coordinates, forces, etc...).
 * \param ff Pointer to structure FORCEFIELD containing forcefield parameters.
 * \param simulCond Pointer to structure SIMULPARAMS containing parameters of the current simulation.
 * \param box Pointer to structure PBC containing Periodic Boundaries Conditions parameters.
 * \param i Index of first atom.
 * \param j Index of second atom.
 * \param r Distance between the two atoms.
 * \param delec Pointer to derivative of energy used for force evaluation.
 * 
 * \brief Evaluates the electrostatic energy and force for a pair when using the SHIFT_2 cutoff.
 *
 * \return On return the electrostatic energy.
 *
 * Shifted electrostatic potential with the shift functional form 2 :
 *
 * \f$ elecShift=elecPot(r)*shiftFunc(r) \f$
 *
 * \f$ elecPot=cte*qi*qj/r \f$
 *
 * Where cte is chgcharmm or chgnam
 *
 * \f$ shiftFunc=1-2r^2/rc^2+r^4/rc^4 \f$
 *
 * \f$ delecShift=delecPot(r)*shiftFunc(r)+elecPot(r)*dshiftFunc(r) \f$
 *
 * \f$ delecPot(r)=-elecPot(r)/r \f$
 *
 * \f$ dshiftFunc(r)=-4r/rc^2+4r^3/rc^4 \f$
 *
 */
double coulomb_shift2(const PARAM *param,double *delec,const double qel,
		      const double r2,const double rt)
{
  double elec=0.,pelec,shift1,shift2,shiftFunc,dshiftFunc;
  
  shift1=r2*param->rcutOff2;
  shift2=X2(shift1);
  
  shiftFunc=1.-2.*shift1+shift2;
  dshiftFunc=4.*(shift2-shift1);
  
  pelec=param->chargeConst*qel*rt;
  elec=pelec*shiftFunc;
  *delec=pelec*rt*(dshiftFunc-shiftFunc);
  
  return elec;
} //END of function

/**
 * \param atom Array of structure ATOM (coordinates, forces, etc...).
 * \param ff Pointer to structure FORCEFIELD containing forcefield parameters.
 * \param simulCond Pointer to structure SIMULPARAMS containing parameters of the current simulation.
 * \param box Pointer to structure PBC containing Periodic Boundaries Conditions parameters.
 * \param i Index of first atom.
 * \param j Index of second atom.
 * \param r Distance between the two atoms.
 * \param delec Pointer to derivative of energy used for force evaluation.
 * 
 * \brief Evaluates the electrostatic energy and force for a pair when using the SWITCH cutoff.
 *
 * \return On return the electrostatic energy.
 *
 * Switched electrostatic potential :
 *
 * \f$ elecSwitch=elecPot(r)*switchFunc(r) \f$
 *
 * \f$ elecPot=cte*qi*qj/r \f$
 *
 * Where cte is chgcharmm or chgnam
 *
 * \f$ switchFunc=(rc^2+2r^2-3ro^2)*(rc^2-r^2)^2/(rc^2-ro^2)^3 \f$
 *
 * \f$ delecSwitch=delecPot(r)*switchFunc(r)+elecPot(r)*dswitchFunc(r) \f$
 *
 * \f$ delecPot(r)=-elecPot(r)/r \f$
 *
 * \f$ dswitchFunc(r)=-12*r*(rc^2-r^2)*(ro^2-r^2)/(rc^2-ro^2)^3 \f$
 *
 */
double coulomb_switch(const PARAM *param,double *delec,const double qel,
		      const double r2,const double rt)
{
  double elec=0.,pelec,switch1,switchFunc,dswitchFunc;
  
  if(r2<=param->cutOn2)
  {
    
    pelec=param->chargeConst*qel*rt;
    elec=pelec;
    *delec=-pelec*rt;
    
  }
  else
  {
    switch1=param->cutOff2-r2;
    
    switchFunc=X2(switch1)*(param->cutOff2+2.*r2-3.*param->cutOn2)*param->switch2;
      
    dswitchFunc=12.*r2*switch1*(param->cutOn2-r2)*param->switch2;
    
    pelec=param->chargeConst*qel*rt;
    elec=pelec*switchFunc;
    *delec=pelec*rt*(dswitchFunc-switchFunc);
    
  }     
  
  return elec;
}

/**
 * \brief Empty function called when 1-4 electrostatic energy and force are disabled.
 *
 * \remarks delec is set to 0.0 .
 *
 * \return On return 0.0 as no energy evaluated.
 */
double coulomb14_none(const PARAM *param,double *delec,const double qel,
		      const double r2,const double rt)
{
  *delec=0.;
  return 0.;
}

/**
 * \param atom Array of structure ATOM (coordinates, forces, etc...).
 * \param ff Pointer to structure FORCEFIELD containing forcefield parameters.
 * \param simulCond Pointer to structure SIMULPARAMS containing parameters of the current simulation.
 * \param box Pointer to structure PBC containing Periodic Boundaries Conditions parameters.
 * \param i Index of first atom.
 * \param j Index of second atom.
 * \param r Distance between the two atoms.
 * \param delec Pointer to derivative of energy used for force evaluation.
 *
 * \brief Function called for a full evaluation of the 1-4 electrostatic energy and force.
 *
 * \return On return the electrostatic energy.
 */
double coulomb14_full(const PARAM *param,double *delec,const double qel,
		      const double r2,const double rt)
{
  
  double elec=0.;
  
  elec=param->scal14*param->chargeConst*qel*rt;
  *delec=-elec*rt;
  
  return elec;
}

/**
 * \param atom Array of structure ATOM (coordinates, forces, etc...).
 * \param ff Pointer to structure FORCEFIELD containing forcefield parameters.
 * \param simulCond Pointer to structure SIMULPARAMS containing parameters of the current simulation.
 * \param box Pointer to structure PBC containing Periodic Boundaries Conditions parameters.
 * \param i Index of first atom.
 * \param j Index of second atom.
 * \param r Distance between the two atoms.
 * \param delec Pointer to derivative of energy used for force evaluation.
 * 
 * \brief Evaluates the 1-4 electrostatic energy and force for a pair when using the SHIFT_1 cutoff.
 *
 * \return On return the electrostatic energy.
 *
 * Shifted 1-4 electrostatic potential with the shift functional form 1:
 *
 * \f$ elecShift=elecPot(r)*shiftFunc(r) \f$
 *
 * \f$ elecPot=cte*qi*qj/r \f$
 *
 * \f$ shiftFunc=1-2r/rc+r^2/rc^2 \f$
 *
 * \f$ delecShift=delecPot(r)*shiftFunc(r)+elecPot(r)*dshiftFunc(r) \f$
 *
 * \f$ delecPot(r)=-elecPot(r)/r \f$
 *
 * \f$ dshiftFunc(r)=-2/rc+2r/rc^2 \f$
 * 
 */
double coulomb14_shift1(const PARAM *param,double *delec,const double qel,
			const double r2,const double rt)
{
  double elec=0.,pelec,shift1,shift2,shiftFunc,dshiftFunc;
    
  if(r2<=param->cutOff2)
  {
    
    shift1=rt*r2*param->rcutOff;
    shift2=r2*param->rcutOff2;
   
    shiftFunc=1.-2.*shift1+shift2;
    dshiftFunc=2.*(shift2-shift1);
    
    pelec=param->scal14*param->chargeConst*qel*rt;
    elec=pelec*shiftFunc;
    *delec=pelec*rt*(dshiftFunc-shiftFunc);
    
  }
  
  return elec;
}

/**
 * \param atom Array of structure ATOM (coordinates, forces, etc...).
 * \param ff Pointer to structure FORCEFIELD containing forcefield parameters.
 * \param simulCond Pointer to structure SIMULPARAMS containing parameters of the current simulation.
 * \param box Pointer to structure PBC containing Periodic Boundaries Conditions parameters.
 * \param i Index of first atom.
 * \param j Index of second atom.
 * \param r Distance between the two atoms.
 * \param delec Pointer to derivative of energy used for force evaluation.
 * 
 * \brief Evaluates the 1-4 electrostatic energy and force for a pair when using the SHIFT_2 cutoff.
 *
 * \return On return the electrostatic energy.
 *
 * Shifted 1-4 electrostatic potential with the shift functional form 2 :
 *
 * \f$ elecShift=elecPot(r)*shiftFunc(r) \f$
 *
 * \f$ elecPot=cte*qi*qj/r \f$
 *
 * \f$ shiftFunc=1-2r^2/rc^2+r^4/rc^4 \f$
 *
 * \f$ delecShift=delecPot(r)*shiftFunc(r)+elecPot(r)*dshiftFunc(r) \f$
 *
 * \f$ delecPot(r)=-elecPot(r)/r \f$
 *
 * \f$ dshiftFunc(r)=-4r/rc^2+4r^3/rc^4 \f$
 *
 */
double coulomb14_shift2(const PARAM *param,double *delec,const double qel,
			const double r2,const double rt)
{
  double elec=0.,pelec,shift1,shift2,shiftFunc,dshiftFunc;
    
  if(r2<=param->cutOff2)
  {
    
    shift1=r2*param->rcutOff2;
    shift2=X2(shift1);
    
    shiftFunc=1.-2.*shift1+shift2;
    dshiftFunc=4.*(shift2-shift1);
    
    pelec=param->scal14*param->chargeConst*qel*rt;
    elec=pelec*shiftFunc;
    *delec=pelec*rt*(dshiftFunc-shiftFunc);
    
  }
  
  return elec;
}

/**
 * \param atom Array of structure ATOM (coordinates, forces, etc...).
 * \param ff Pointer to structure FORCEFIELD containing forcefield parameters.
 * \param simulCond Pointer to structure SIMULPARAMS containing parameters of the current simulation.
 * \param box Pointer to structure PBC containing Periodic Boundaries Conditions parameters.
 * \param i Index of first atom.
 * \param j Index of second atom.
 * \param r Distance between the two atoms.
 * \param delec Pointer to derivative of energy used for force evaluation.
 * 
 * \brief Evaluates the 1-4 electrostatic energy and force for a pair when using the SWITCH cutoff.
 *
 * \return On return the electrostatic energy.
 *
 * Switched 1-4 electrostatic potential :
 *
 * \f$ elecSwitch=elecPot(r)*switchFunc(r) \f$
 *
 * \f$ elecPot=cte*qi*qj/r \f$
 *
 * \f$ switchFunc=(rc^2+2r^2-3ro^2)*(rc^2-r^2)^2/(rc^2-ro^2)^3 \f$
 *
 * \f$ delecSwitch=delecPot(r)*switchFunc(r)+elecPot(r)*dswitchFunc(r) \f$
 *
 * \f$ delecPot(r)=-elecPot(r)/r \f$
 *
 * \f$ dswitchFunc(r)=-12*r*(rc^2-r^2)*(ro^2-r^2)/(rc^2-ro^2)^3 \f$
 * 
 */
double coulomb14_switch(const PARAM *param,double *delec,const double qel,
			const double r2,const double rt)
{
  double elec=0.,pelec,switch1,switchFunc,dswitchFunc;
    
  if(r2<=param->cutOn2)
  {
    pelec=param->scal14*param->chargeConst*qel*rt;
    elec=pelec;
    *delec=-pelec*rt;
     
  }
  else if(r2<=param->cutOff2)
  {
    
    switch1=param->cutOff2-r2;
    
    switchFunc=X2(switch1)*(param->cutOff2+2.*r2-3.*param->cutOn2)*param->switch2;
      
    dswitchFunc=12.*r2*switch1*(param->cutOn2-r2)*param->switch2;
    
    pelec=param->scal14*param->chargeConst*qel*rt;
    elec=pelec*switchFunc;
    *delec=pelec*rt*(dswitchFunc-switchFunc);
      
  }     
  return elec;
}
