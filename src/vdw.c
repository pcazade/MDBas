/*
 * Copyright (c) 2013 Pierre-Andre Cazade
 * Copyright (c) 2013 Florent hedin
 *
 * This file is part of MDBas.
 *
 * MDBas is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MDBas is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with MDBas.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * \file vdw.c
 * \brief Contains functions for evaluating Van der Waals energies and forces.
 * \author Pierre-Andre Cazade and Florent Hedin
 * \version alpha-branch
 * \date 2012
 */

#include <math.h>

#include "global.h"
#include "utils.h"

/**
 * \param atom Array of structure ATOM (coordinates, forces, etc...).
 * \param ff Pointer to structure FORCEFIELD containing forcefield parameters.
 * \param simulCond Pointer to structure SIMULPARAMS containing parameters of the current simulation.
 * \param box Pointer to structure PBC containing Periodic Boundaries Conditions parameters.
 * \param i Index of first atom.
 * \param j Index of second atom.
 * \param r Distance between the two atoms.
 * \param dvdw Pointer to derivative of energy used for force evaluation.
 *
 * \brief Empty function called when Van der Waals energy and force are disabled.
 *
 * \remarks dvdw is set to 0.0 .
 *
 * \return On return 0.0 as no energy evaluated.
 */
double vdw_none(const PARAM *param,double *dvdw,const double veps,
                const double vsig,const double r2, const double rt)
{
    *dvdw=0.;
    return 0.;
}

/**
 * \param atom Array of structure ATOM (coordinates, forces, etc...).
 * \param ff Pointer to structure FORCEFIELD containing forcefield parameters.
 * \param ener Pointer to structure ENERGY containing values of the different energies.
 * \param simulCond Pointer to structure SIMULPARAMS containing parameters of the current simulation.
 * \param box Pointer to structure PBC containing Periodic Boundaries Conditions parameters.
 *
 * \brief Function called for a full evaluation of the Van der Waals energy and force.
 */
void vdw_full(const PARAM *param,PARALLEL *parallel, ENERGY *ener, const PBC *box,double *x,
              double *y,double *z,double *fx, double *fy, double *fz,double *eps,double *sig,
              int **exclList,int *exclPair)
{

    int i,j,k,exclude;;
    double evdw=0.,pvdw,dvdw;
    double r,r2,rt,fxj,fyj,fzj,fxi,fyi,fzi;
    double delta[3];

    for(i=0; i<param->nAtom-1; i++)
    {
        fxi=0.;
        fyi=0.;
        fzi=0.;
	
	k=0;

        for(j=i+1; j<param->nAtom; j++)
        {
	  
	    exclude=0;
	    if( (exclPair[i]>0) && (exclList[i][k]==j) )
	    {
	      exclude=1;
	      k++;
	      
	      if(k>=exclPair[i])
		k=exclPair[i]-1;
	    }

            if(!exclude)
            {

                delta[0]=x[j]-x[i];
                delta[1]=y[j]-y[i];
                delta[2]=z[j]-z[i];

                r2=dist(box,delta);
                r=sqrt(r2);
                rt=1./r;

                pvdw=4.*eps[i]*eps[j]*(X12((sig[i]+sig[j])*rt)-X6((sig[i]+sig[j])*rt));

                dvdw=24.*eps[i]*eps[j]*rt*(X6((sig[i]+sig[j])*rt)-2.*X12((sig[i]+sig[j])*rt));

                evdw+=pvdw;

                fxj=dvdw*delta[0]*rt;
                fyj=dvdw*delta[1]*rt;
                fzj=dvdw*delta[2]*rt;

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
    ener->vdw+=evdw;
}

/**
 * \param atom Array of structure ATOM (coordinates, forces, etc...).
 * \param ff Pointer to structure FORCEFIELD containing forcefield parameters.
 * \param simulCond Pointer to structure SIMULPARAMS containing parameters of the current simulation.
 * \param box Pointer to structure PBC containing Periodic Boundaries Conditions parameters.
 * \param i Index of first atom.
 * \param j Index of second atom.
 * \param r Distance between the two atoms.
 * \param dvdw Pointer to derivative of energy used for force evaluation.
 *
 * \brief Evaluates the Van der Waals energy and force for a pair when using the SWITCH cutoff.
 *
 * \return On return the Van der Waals energy.
 *
 * Switched Van der Waals potential :
 *
 * \f$ vdwSwitch=elecPot(r)*switchFunc(r) \f$
 *
 * \f$ vdwPot=4*eps*((sig/r)^12-(sig/r)^6) \f$
 *
 * \f$ switchFunc=(rc^2+2r^2-3ro^2)*(rc^2-r^2)^2/(rc^2-ro^2)^3 \f$
 *
 * \f$ dvdwSwitch=delecPot(r)*switchFunc(r)+elecPot(r)*dswitchFunc(r) \f$
 *
 * \f$ dvdwPot(r)=-elecPot(r)/r \f$
 *
 * \f$ dswitchFunc(r)=-12*r*(rc^2-r^2)*(ro^2-r^2)/(rc^2-ro^2)^3 \f$
 *
 */
double vdw_switch(const PARAM *param,double *dvdw,const double veps,
                  const double vsig,const double r2, const double rt)
{
    double evdw=0.,pvdw,dpvdw,switch1,switchFunc,dswitchFunc;
    double vsig6,vsig12;

    vsig6=vsig*rt;
    vsig6=X6(vsig6);
    vsig12=X2(vsig6);

    if(r2<=param->cutOn2)
    {

        pvdw=4.*veps*(vsig12-vsig6);

        *dvdw=24.*veps*rt*(vsig6-2.*vsig12);

        evdw=pvdw;

    }
    else
    {

        switch1=param->cutOff2-r2;

        switchFunc=X2(switch1)*(param->cutOff2+2.*r2-3.*param->cutOn2)*param->switch2;

        dswitchFunc=12.*r2*switch1*(param->cutOn2-r2)*param->switch2;

        pvdw=4.*veps*(vsig12-vsig6);

        dpvdw=24.*veps*(vsig6-2.*vsig12);

        evdw=pvdw*switchFunc;

        *dvdw=pvdw*dswitchFunc+dpvdw*switchFunc;
        *dvdw*=rt;
    }

    return evdw;
}

/**
 * \param atom Array of structure ATOM (coordinates, forces, etc...).
 * \param ff Pointer to structure FORCEFIELD containing forcefield parameters.
 * \param simulCond Pointer to structure SIMULPARAMS containing parameters of the current simulation.
 * \param box Pointer to structure PBC containing Periodic Boundaries Conditions parameters.
 * \param i Index of first atom.
 * \param j Index of second atom.
 * \param r Distance between the two atoms.
 * \param dvdw Pointer to derivative of energy used for force evaluation.
 *
 * \brief Empty function called when 1-4 Van der Waals energy and force are disabled.
 *
 * \remarks dvdw is set to 0.0 .
 *
 * \return On return 0.0 as no energy evaluated.
 */
double vdw14_none(const PARAM *param,double *dvdw,const double veps,
                  const double vsig,const double r2, const double rt)
{
    *dvdw=0;
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
 * \param dvdw Pointer to derivative of energy used for force evaluation.
 *
 * \brief Function called for a full evaluation of the 1-4 Van der Waals energy and force.
 *
 * \return On return the electrostatic energy.
 */
double vdw14_full(const PARAM *param,double *dvdw,const double veps,
                  const double vsig,const double r2, const double rt)
{

    double evdw=0.;
    double vsig6,vsig12;

    vsig6=vsig*rt;
    vsig6=X6(vsig6);
    vsig12=X2(vsig6);

    evdw=4.*param->scal14*veps*(vsig12-vsig6);

    *dvdw=24.*param->scal14*veps*rt*(vsig6-2.*vsig12);

    return evdw;
}

/**
 * \param atom Array of structure ATOM (coordinates, forces, etc...).
 * \param ff Pointer to structure FORCEFIELD containing forcefield parameters.
 * \param simulCond Pointer to structure SIMULPARAMS containing parameters of the current simulation.
 * \param box Pointer to structure PBC containing Periodic Boundaries Conditions parameters.
 * \param i Index of first atom.
 * \param j Index of second atom.
 * \param r Distance between the two atoms.
 * \param dvdw Pointer to derivative of energy used for force evaluation.
 *
 * \brief Evaluates the 1-4 Van der Waals energy and force for a pair when using the SWITCH cutoff.
 *
 * \return On return the Van der Waals energy.
 *
 * Switched 1-4 Van der Waals potential :
 *
 * \f$ vdwSwitch=elecPot(r)*switchFunc(r) \f$
 *
 * \f$ vdwPot=4*eps*((sig/r)^12-(sig/r)^6) \f$
 *
 * \f$ switchFunc=(rc^2+2r^2-3ro^2)*(rc^2-r^2)^2/(rc^2-ro^2)^3 \f$
 *
 * \f$ dvdwSwitch=delecPot(r)*switchFunc(r)+elecPot(r)*dswitchFunc(r) \f$
 *
 * \f$ dvdwPot(r)=-elecPot(r)/r \f$
 *
 * \f$ dswitchFunc(r)=-12*r*(rc^2-r^2)*(ro^2-r**2)/(rc^2-ro^2)^3 \f$
 *
 */
double vdw14_switch(const PARAM *param,double *dvdw,const double veps,
                    const double vsig,const double r2, const double rt)
{
    double evdw=0.,pvdw,dpvdw,switch1,switchFunc,dswitchFunc;
    double vsig6,vsig12;

    vsig6=vsig*rt;
    vsig6=X6(vsig6);
    vsig12=X2(vsig6);

    if(r2<=param->cutOn2)
    {

        pvdw=4.*param->scal14*veps*(vsig12-vsig6);

        *dvdw=24.*param->scal14*veps*rt*(vsig6-2.*vsig12);

        evdw=pvdw;

    }
    else if(r2<=param->cutOff2)
    {

        switch1=param->cutOff2-r2;

        switchFunc=X2(switch1)*(param->cutOff2+2.*r2-3.*param->cutOn2)*param->switch2;

        dswitchFunc=12.*r2*switch1*(param->cutOn2-r2)*param->switch2;

        pvdw=4.*param->scal14*veps*(vsig12-vsig6);

        dpvdw=24.*param->scal14*veps*(vsig6-2.*vsig12);

        evdw=pvdw*switchFunc;

        *dvdw=pvdw*dswitchFunc+dpvdw*switchFunc;
        *dvdw*=rt;

    }

    return evdw;
}
