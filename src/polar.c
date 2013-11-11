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
static double *polTensor;

void init_polar(PARAM *param)
{
  
  elFieldX=(double*)my_malloc(param->nAtom*sizeof(*elFieldX));
  elFieldY=(double*)my_malloc(param->nAtom*sizeof(*elFieldY));
  elFieldZ=(double*)my_malloc(param->nAtom*sizeof(*elFieldZ));
  
  elFieldIndX=(double*)my_malloc(param->nAtom*sizeof(*elFieldIndX));
  elFieldIndY=(double*)my_malloc(param->nAtom*sizeof(*elFieldIndY));
  elFieldIndZ=(double*)my_malloc(param->nAtom*sizeof(*elFieldIndZ));
  
  muIndX=(double*)my_malloc(param->nAtom*sizeof(*muIndX));
  muIndY=(double*)my_malloc(param->nAtom*sizeof(*muIndY));
  muIndZ=(double*)my_malloc(param->nAtom*sizeof(*muIndZ));
  
  muOldX=(double*)my_malloc(param->nAtom*sizeof(*muOldX));
  muOldY=(double*)my_malloc(param->nAtom*sizeof(*muOldY));
  muOldZ=(double*)my_malloc(param->nAtom*sizeof(*muOldZ));
  
}

void free_polar()
{
  
  free(elFieldX);
  free(elFieldY);
  free(elFieldZ);
  
  free(elFieldIndX);
  free(elFieldIndY);
  free(elFieldIndZ);
  
  free(muIndX);
  free(muIndY);
  free(muIndZ);
  
  free(muOldX);
  free(muOldY);
  free(muOldZ);
  
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

double polar_ener_iter(CTRL *ctrl,PARAM *param,PARALLEL *parallel,PBC *box, POLAR *polar,
		       const double x[],const double y[],const double z[],double fx[],
		       double fy[],double fz[],const double q[],
		       const double alPol[],const int polList,int **neighList,
		       const int neighPair[],double dBuffer[])
{
  int i,j,k,l;
  int icycle;
  double epol,qi,qj;
  double muVar,converged;
  double efxi,efyi,efzi;
  double fxi,fyi,fzi;
  double r2,rt,rt2,rt3,rt5;
  double tm1,tm2,tm3,tm4,tm5,tm6;
  double tmxi,tmyi,tmzi,tmxj,tmyj,tmzj;
  double mrt2,mrt5,mufi,mufj,mufij;
  double muftx,mufty,muftz;
  double delta[3];
  
  static_field(ctrl,param,parallel,box,x,y,z,q,neighList,neighPair,dBuffer);
  
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
    
    fxi=0.;
    fyi=0.;
    fzi=0.;
    
    qi=param->chargeConst*q[i];
    
    for(k=0; k<neighPair[l]; k++)
    {
      j=neighList[l][k];
      
      if(!polList[j])
	continue;
      
      qj=param->chargeConst*q[j];

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
	
	mrt2=5.0*rt2;
	mrt5=3.0*param->chargeConst*rt5;
	
	tm1=delta[0]*rt5;
	tm3=delta[1]*rt5;
	tm6=delta[2]*rt5;
	
	tm2=tm1*delta[1];
	tm4=tm1*delta[2];
	tm5=tm3*delta[2];
	
	tm1=tm1*delta[0]-rt3;
	tm3=tm3*delta[1]-rt3;
	tm6=tm6*delta[2]-rt3;
	
	tmxi=tm1*muIndX[i]+tm2*muIndY[i]+tm4*muIndZ[i];
	tmyi=tm2*muIndX[i]+tm3*muIndY[i]+tm5*muIndZ[i];
	tmzi=tm4*muIndX[i]+tm5*muIndY[i]+tm6*muIndZ[i];
	
	tmxj=tm1*muIndX[j]+tm2*muIndY[j]+tm4*muIndZ[j];
	tmyj=tm2*muIndX[j]+tm3*muIndY[j]+tm5*muIndZ[j];
	tmzj=tm4*muIndX[j]+tm5*muIndY[j]+tm6*muIndZ[j];
	
	mufi=muIndX[i]*delta[0]+muIndY[i]*delta[1]+muIndZ[i]*delta[2];
	mufj=muIndX[j]*delta[0]+muIndY[j]*delta[1]+muIndZ[j]*delta[2];
	mufij=muIndX[i]*muIndX[j]+muIndY[i]*muIndY[j]+muIndZ[i]*muIndZ[j];
	
	muftx=mrt5*(mrt2*mufi*mufj*delta[0]-mufij*delta[0]-mufi*muIndX[j]-mufj*muIndX[i]);
	mufty=mrt5*(mrt2*mufi*mufj*delta[1]-mufij*delta[1]-mufi*muIndY[j]-mufj*muIndY[i]);
	muftz=mrt5*(mrt2*mufi*mufj*delta[2]-mufij*delta[2]-mufi*muIndZ[j]-mufj*muIndZ[i]);
	
	fx[j]+=qj*tmxi-qi*tmxj+muftx;
	fy[j]+=qj*tmyi-qi*tmyj+mufty;
	fz[j]+=qj*tmzi-qi*tmzj+muftz;
	
	fxi+=qi*tmxj-qj*tmxi-muftx;
	fyi+=qi*tmyj-qj*tmyi-mufty;
	fzi+=qi*tmzj-qj*tmzi-muftz;
      }
      
    }
    
    epol+=muIndX[i]*elFieldX[i]+muIndY[i]*elFieldY[i]+muIndZ[i]*elFieldZ[i];
    
    fx[i]+=fxi;
    fy[i]+=fyi;
    fz[i]+=fzi;
    
  }
  
  if(parallel->nProc>1)
  { 
    sum_double_para(&epol,dBuffer,1);
  }
  
  epol*=0.5*param->chargeConst;
  
  return(epol);
  
}

double polar_ener_iter(CTRL *ctrl,PARAM *param,PARALLEL *parallel,PBC *box, POLAR *polar,
		       const double x[],const double y[],const double z[],double fx[],
		       double fy[],double fz[],const double q[],const double alPol[],
		       const int polList,int **neighList,const int neighPair[],
		       double dBuffer[])
{
  int i,j,k,l;
  double epol,qi,qj;
  double fxi,fyi,fzi;
  double r2,rt,rt2,rt3,rt5;
  double tm1,tm2,tm3,tm4,tm5,tm6;
  double delta[3];
  
  static_field(ctrl,param,parallel,box,x,y,z,q,neighList,neighPair,dBuffer);
  
  
}