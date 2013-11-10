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

#include <math.h>

#include "global.h"
#include "utils.h"
#include "memory.h"

#ifdef USING_MPI
#include "parallel.h"
#else
#include "serial.h"
#endif

static double *elFieldX,*elFieldY,*elFieldZ;
static double *elFieldIndX,*elFieldIndY,*elFieldIndZ;
static double *muIndX,*muIndY,*muIndZ;
static double *muOldX,*muOldY,*muOldZ;

void init_polar(CTRL *ctrl,PARAM *param,PARALLEL *parallel,PBC *box)
{
  
}

void static_field(CTRL *ctrl,PARAM *param,PARALLEL *parallel,PBC *box,
		  const double x[],const double y[],const double z[],
		  const double q[],int **neighList,const int neighPair[],
		  double dBuffer[])
{
  
  int i,j,k,l;
  double qri,qrj;
  double efxi,efyi,efzi;
  double r2,rt,rt2,rt3;
  double delta[3];
  
  for(i=0;i<param->nAtom;i++)
  {
    elFieldX[i]=0.;
    elFieldY[i]=0.;
    elFieldZ[i]=0.;
  }
  
  l=0;
  for(i=parallel->idProc; i<param->nAtom; i+=parallel->nProc)
  {
    
    efxi=0.;
    efyi=0.;
    efzi=0.;
    
    for(k=0; k<neighPair[l]; k++)
    {
      j=neighList[l][k];

      delta[0]=x[j]-x[i];
      delta[1]=y[j]-y[i];
      delta[2]=z[j]-z[i];

      r2=dist(box,delta);

      if(r2<=param->cutOff2)
      {
	rt2=1./r2;
	rt=sqrt(rt2);
	rt3=rt*rt2;
	
	qri=rt3*q[i];
	qrj=rt3*q[j];
	
	elFieldX[j]+=qri*delta[0];
	elFieldY[j]+=qri*delta[1];
	elFieldZ[j]+=qri*delta[2];
	
	efxi+=qrj*delta[0];
	efyi+=qrj*delta[1];
	efzi+=qrj*delta[2];
	
      } // if(r2<=param->cutOff2)
      
    } // end for(k=0; k<neighPair[l]; k++)
    
    elFieldX[i]-=efxi;
    elFieldY[i]-=efyi;
    elFieldZ[i]-=efzi;
    
    l++;
  } // end for(i=parallel->idProc; i<param->nAtom; i+=parallel->nProc)
  
  if(parallel->nProc>1)
  {
    sum_double_para(elFieldX,dBuffer,param->nAtom);
    sum_double_para(elFieldY,dBuffer,param->nAtom);
    sum_double_para(elFieldZ,dBuffer,param->nAtom);
  }
  
}

double epol(CTRL *ctrl,PARAM *param,PARALLEL *parallel,PBC *box, POLAR *polar,
	    const double x[],const double y[],const double z[],
	    const double alPol[],const int polList,int **neighList,
	    const int neighPair[],double dBuffer[])
{
  int i,j,k,l;
  int icycle;
  double epol;
  double muVar,converged;
  double efxi,efyi,efzi;
  double r2,rt,rt2,rt3,rt5;
  double tm1,tm2,tm3,tm4,tm5,tm6;
  double delta[3];
  
  epol=0;
  
  for(i=0;i<param->nAtom;i++)
  {
    muIndX[i]=0.;
    muIndY[i]=0.;
    muIndZ[i]=0.;
    
    muOldX[i]=0.;
    muOldY[i]=0.;
    muOldZ[i]=0.;
  }
  
  for(i=parallel->idProc; i<param->nAtom; i+=parallel->nProc)
  {
    if(!polList[i])
      continue;
    
    muIndX[i]=alPol[i]*elFieldX[i];
    muIndY[i]=alPol[i]*elFieldY[i];
    muIndZ[i]=alPol[i]*elFieldZ[i];
    
  }
  
  if(parallel->nProc>1)
  {
    sum_double_para(muIndX,dBuffer,param->nAtom);
    sum_double_para(muIndY,dBuffer,param->nAtom);
    sum_double_para(muIndZ,dBuffer,param->nAtom);
  }
  
  icycle=1;
  converged=0;
  do
  {
  
    for(i=0;i<param->nAtom;i++)
    {
      elFieldIndX[i]=0.;
      elFieldIndY[i]=0.;
      elFieldIndZ[i]=0.;
    }
    
    l=0;
    for(i=parallel->idProc; i<param->nAtom; i+=parallel->nProc)
    {
      if(!polList[i])
	continue;
      
      efxi=0.;
      efyi=0.;
      efzi=0.;
      
      for(k=0; k<neighPair[l]; k++)
      {
	j=neighList[l][k];
	
	if(!polList[j])
	  continue;

	delta[0]=x[j]-x[i];
	delta[1]=y[j]-y[i];
	delta[2]=z[j]-z[i];

	r2=dist(box,delta);

	if(r2<=param->cutOff2)
	{
	  rt2=1./r2;
	  rt=sqrt(rt2);
	  rt3=rt*rt2;
	  rt5=rt3*rt2;
	  
	  tm1=delta[0]*rt5;
	  tm3=delta[1]*rt5;
	  tm6=delta[2]*rt5;
	  
	  tm2=tm1*delta[1];
	  tm4=tm1*delta[2];
	  tm5=tm3*delta[2];
	  
	  tm1=tm1*delta[0]-rt3;
	  tm3=tm3*delta[1]-rt3;
	  tm6=tm6*delta[2]-rt3;
	  
	  elFieldIndX[j]+=muIndX[i]*tm1+muIndY[i]*tm2+muIndZ[i]*tm4;
	  elFieldIndY[j]+=muIndX[i]*tm2+muIndY[i]*tm3+muIndZ[i]*tm5;
	  elFieldIndZ[j]+=muIndX[i]*tm4+muIndY[i]*tm5+muIndZ[i]*tm6;
	  
	  efxi=muIndX[j]*tm1+muIndY[j]*tm2+muIndZ[j]*tm4;
	  efyi=muIndX[j]*tm2+muIndY[j]*tm3+muIndZ[j]*tm5;
	  efzi=muIndX[j]*tm4+muIndY[j]*tm5+muIndZ[j]*tm6;
	  
	} // if(r2<=param->cutOff2)
	
      } // end for(k=0; k<neighPair[l]; k++)
      
      elFieldIndX[i]+=efxi;
      elFieldIndY[i]+=efyi;
      elFieldIndZ[i]+=efzi;
      
      l++;
    } // end for(i=parallel->idProc; i<param->nAtom; i+=parallel->nProc)
    
    if(parallel->nProc>1)
    {
      sum_double_para(elFieldIndX,dBuffer,param->nAtom);
      sum_double_para(elFieldIndY,dBuffer,param->nAtom);
      sum_double_para(elFieldIndZ,dBuffer,param->nAtom);
    }
    
    muVar=0.;
    for(i=parallel->idProc; i<param->nAtom; i+=parallel->nProc)
    {
      if(!polList[i])
	continue;
      
      muIndX[i]=alPol[i]*(elFieldX[i]+elFieldIndX[i]);
      muIndY[i]=alPol[i]*(elFieldY[i]+elFieldIndY[i]);
      muIndZ[i]=alPol[i]*(elFieldZ[i]+elFieldIndZ[i]);
      
      muVar+=X2(muIndX[i]-muOldX[i])+X2(muIndY[i]-muOldY[i])+
	     X2(muIndZ[i]-muOldZ[i]);
      
    }
    
    if(parallel->nProc>1)
    {
      sum_double_para(muIndX,dBuffer,param->nAtom);
      sum_double_para(muIndY,dBuffer,param->nAtom);
      sum_double_para(muIndZ,dBuffer,param->nAtom);
      
      sum_double_para(&muVar,dBuffer,1);
    }
    
    if(muVar<=polar->tol)
      converged=1;
    
    if(!converged)
    {
      for(i=0;i<param->nAtom;i++)
      {
	muOldX[i]=muIndX[i];
	muOldY[i]=muIndY[i];
	muOldZ[i]=muIndZ[i];
      }
    }
    
    icycle++;
    
  }while((!converged)&&(icycle<=polar->maxCycle));
  
  epol=0.;
  
  for(i=parallel->idProc; i<param->nAtom; i+=parallel->nProc)
  {
    if(!polList[i])
      continue;
    
    epol+=muIndX[i]*elFieldX[i]+muIndY[i]*elFieldY[i]+muIndZ[i]*elFieldZ[i];
    
  }
  
  epol*=0.5*param->chargeConst;
  
  return(epol);
  
}

double polar_ener()
{
  
}