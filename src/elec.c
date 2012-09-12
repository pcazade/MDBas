/**
 * \file elec.c
 * \brief Contains functions for evaluating electrostatic energies and forces.
 * \author Pierre-Andre Cazade and Florent Hedin
 * \version alpha-branch
 * \date 2012
 */

#include "global.h"
#include "utils.h"

/**
 * \brief Empty function called when electrostatic energy and force are disabled.
 * 
 * \remarks delec is set to 0.0 .
 *
 * \return On return 0.0 as no energy evaluated.
 */
double coulomb_none(ATOM atom[],FORCEFIELD *ff,SIMULPARAMS *simulCond,PBC *box,
		    int i, int j, double r, double *delec)
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
void coulomb_full(ATOM atom[],FORCEFIELD *ff,ENERGY *ener,SIMULPARAMS *simulCond,PBC *box)
{
  
  int i,j,k,exclude;
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
      
    }
    
    atom[i].fx+=fxi;
    atom[i].fy+=fyi;
    atom[i].fz+=fzi;
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
 */
double coulomb_shift1(ATOM atom[],FORCEFIELD *ff,SIMULPARAMS *simulCond,PBC *box,
		    int i, int j, double r, double *delec)
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
  
  double elec=0.,pelec,shiftFunc,dshiftFunc;

  shiftFunc=1.-2.*r/simulCond->cutoff+X2(r/simulCond->cutoff);
  dshiftFunc=-2./simulCond->cutoff+2.*r/X2(simulCond->cutoff);
    
  pelec=simulCond->chargeConst*atom[i].q*atom[j].q/r;
  elec=pelec*shiftFunc;
  *delec=pelec*(dshiftFunc-shiftFunc/r);
  
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
 */
double coulomb_shift2(ATOM atom[],FORCEFIELD *ff,SIMULPARAMS *simulCond,PBC *box,
		    int i, int j, double r, double *delec)
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
  
  double elec=0.,pelec,shiftFunc,dshiftFunc;
  
  shiftFunc=1.-2.*X2(r/simulCond->cutoff)+X4(r/simulCond->cutoff);
  dshiftFunc=-4.*r/X2(simulCond->cutoff)+4.*X3(r)/X4(simulCond->cutoff);
  
  pelec=simulCond->chargeConst*atom[i].q*atom[j].q/r;
  elec=pelec*shiftFunc;
  *delec=pelec*(dshiftFunc-shiftFunc/r);
  
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
 */
double coulomb_switch(ATOM atom[],FORCEFIELD *ff,SIMULPARAMS *simulCond,PBC *box,
		    int i, int j, double r, double *delec)
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
  
  double elec=0.,pelec,switchFunc,dswitchFunc;
  
  if(r<=simulCond->cuton)
  {
    
    pelec=simulCond->chargeConst*atom[i].q*atom[j].q/r;
    elec=pelec;
    *delec=-pelec/r;
    
  }
  else
  {
    
    switchFunc=X2(X2(simulCond->cutoff)-X2(r))*
      (X2(simulCond->cutoff)+2.*X2(r)-3.*X2(simulCond->cuton))/
      X3(X2(simulCond->cutoff)-X2(simulCond->cuton));
      
    dswitchFunc=12.*r*(X2(simulCond->cutoff)-X2(r))*
      (X2(simulCond->cuton)-X2(r))/
      X3(X2(simulCond->cutoff)-X2(simulCond->cuton));
    
    pelec=simulCond->chargeConst*atom[i].q*atom[j].q/r;
    elec=pelec*switchFunc;
    *delec=pelec*(dswitchFunc-switchFunc/r);
    
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
double coulomb14_none(ATOM atom[],FORCEFIELD *ff,SIMULPARAMS *simulCond,PBC *box,
		      int i, int j, double r, double *delec)
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
double coulomb14_full(ATOM atom[],FORCEFIELD *ff,SIMULPARAMS *simulCond,PBC *box,
		    int i, int j, double r, double *delec)
{
  
  double elec=0.;
  
  elec=ff->scal14*simulCond->chargeConst*atom[i].q*atom[j].q/r;
  *delec=-elec/r;
  
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
 */
double coulomb14_shift1(ATOM atom[],FORCEFIELD *ff,SIMULPARAMS *simulCond,PBC *box,
		      int i, int j, double r, double *delec)
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
  
  double elec=0.,pelec,shiftFunc,dshiftFunc;
    
  if(r<=simulCond->cutoff)
  {
    
    shiftFunc=1.-2.*r/simulCond->cutoff+X2(r/simulCond->cutoff);
    dshiftFunc=-2./simulCond->cutoff+2.*r/X2(simulCond->cutoff);
    
    pelec=ff->scal14*simulCond->chargeConst*atom[i].q*atom[j].q/r;
    elec=pelec*shiftFunc;
    *delec=pelec*(dshiftFunc-shiftFunc/r);
    
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
 */
double coulomb14_shift2(ATOM atom[],FORCEFIELD *ff,SIMULPARAMS *simulCond,PBC *box,
		      int i, int j, double r, double *delec)
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
  
  double elec=0.,pelec,shiftFunc,dshiftFunc;
    
  if(r<=simulCond->cutoff)
  {
    
    shiftFunc=1.-2.*X2(r/simulCond->cutoff)+X4(r/simulCond->cutoff);
    dshiftFunc=-4.*r/X2(simulCond->cutoff)+4.*X3(r)/X4(simulCond->cutoff);
    
    pelec=ff->scal14*simulCond->chargeConst*atom[i].q*atom[j].q/r;
    elec=pelec*shiftFunc;
    *delec=pelec*(dshiftFunc-shiftFunc/r);
    
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
 */
double coulomb14_switch(ATOM atom[],FORCEFIELD *ff,SIMULPARAMS *simulCond,PBC *box,
		      int i, int j, double r, double *delec)
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
  
  double elec=0.,pelec,switchFunc,dswitchFunc;
    
  if(r<=simulCond->cuton)
  {
    
    pelec=ff->scal14*simulCond->chargeConst*atom[i].q*atom[j].q/r;
    elec=pelec;
    *delec=-pelec/r;
      
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
    elec=pelec*switchFunc;
    *delec=pelec*(dswitchFunc-switchFunc/r);
      
  }     
  return elec;
}
