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
#include <float.h>

#include "global.h"
#include "polar.h"
#include "utils.h"
#include "memory.h"
#include "errors.h"

#ifdef USING_MPI
#include "parallel.h"
#else
#include "serial.h"
#endif

extern void dtptri_(char* UPLO,char* DIAG,int* N,real *AP,int *INFO);

static int *polList,*polMap;

static real *elFieldX,*elFieldY,*elFieldZ;
static real *elFieldIndX,*elFieldIndY,*elFieldIndZ;
static real *muIndX,*muIndY,*muIndZ;
static real *muOldX,*muOldY,*muOldZ;
static real *polTensor;

void init_polar(CTRL *ctrl,PARAM *param,POLAR *polar,const real alPol[])
{
  
  int i,ii;
  
  elFieldX=(real*)my_malloc(param->nAtom*sizeof(*elFieldX));
  elFieldY=(real*)my_malloc(param->nAtom*sizeof(*elFieldY));
  elFieldZ=(real*)my_malloc(param->nAtom*sizeof(*elFieldZ));
  
  elFieldIndX=(real*)my_malloc(param->nAtom*sizeof(*elFieldIndX));
  elFieldIndY=(real*)my_malloc(param->nAtom*sizeof(*elFieldIndY));
  elFieldIndZ=(real*)my_malloc(param->nAtom*sizeof(*elFieldIndZ));
  
  muIndX=(real*)my_malloc(param->nAtom*sizeof(*muIndX));
  muIndY=(real*)my_malloc(param->nAtom*sizeof(*muIndY));
  muIndZ=(real*)my_malloc(param->nAtom*sizeof(*muIndZ));
  
  muOldX=(real*)my_malloc(param->nAtom*sizeof(*muOldX));
  muOldY=(real*)my_malloc(param->nAtom*sizeof(*muOldY));
  muOldZ=(real*)my_malloc(param->nAtom*sizeof(*muOldZ));
  
  polList=(int*)my_malloc(param->nAtom*sizeof(*polList));
  polMap=(int*)my_malloc(param->nAtom*sizeof(*polMap));
  
  for(i=0;i<param->nAtom;i++)
  {
    if(alPol[i]>DBL_EPSILON)
      polList[i]=1;
    else
      polList[i]=0;
  }
  
  if(ctrl->keyPolInv)
  {
    
    ii=0;
    polar->nAtPol=0;
    for(i=0;i<param->nAtom;i++)
    {
      
      if(polList[i])
      {
	
	polMap[i]=ii;
	
	polar->nAtPol++;
	ii++;
	
      }// end if(polList[i])
      
    } // end for(i=0;i<param->nAtom;i++)
    
    polar->nPolTensor=polar->nAtPol*(polar->nAtPol+1)/2;
    
    polTensor=(real*)my_malloc(polar->nPolTensor*sizeof(*polTensor));
    
  } // end if(ctrl->keyPolInv)
  
}

void free_polar(CTRL *ctrl)
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
  
  free(polList);
  free(polMap);
  
  if(ctrl->keyPolInv)
    free(polTensor);
  
}

void static_field(PARAM *param,PARALLEL *parallel,PBC *box,
		  const real x[],const real y[],const real z[],
		  const real q[],int **neighList,const int neighPair[],
		  real dBuffer[])
{
  
  int i,j,k,l;
  real qri,qrj;
  real efxi,efyi,efzi;
  real r2,rt,rt2,rt3;
  real delta[3];
  
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

void polar_ener_iter(PARAM *param,PARALLEL *parallel,ENERGY *ener,PBC *box, POLAR *polar,
		     const real x[],const real y[],const real z[],real fx[],
		     real fy[],real fz[],const real q[],const real alPol[],
		     int **neighList,const int neighPair[],real dBuffer[])
{
  int i,j,k,l;
  int icycle;
  real epol,virpol;
  real muVar,converged;
  real efxi,efyi,efzi;
  real r2,rt,rt2,rt3,rt5;
  real tm1,tm2,tm3,tm4,tm5,tm6;
  real delta[3];
  
  static_field(param,parallel,box,x,y,z,q,neighList,neighPair,dBuffer);
  
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
  
  polar_forces(param,parallel,box,x,y,z,fx,fy,fz,q,virpol,neighList,neighPair);
  
  epol=0.;
  for(i=parallel->idProc; i<param->nAtom; i+=parallel->nProc)
  {
    if(!polList[i])
      continue;
    
    epol+=muIndX[i]*elFieldX[i]+muIndY[i]*elFieldY[i]+muIndZ[i]*elFieldZ[i];
    
  }
  
  epol*=0.5*param->chargeConst;
  
  ener->epol=epol;
  ener->virpol=virpol;
  
}

void polar_ener_inv(PARAM *param,PARALLEL *parallel,ENERGY *ener,PBC *box, POLAR *polar,
		    const real x[],const real y[],const real z[],real fx[],
		    real fy[],real fz[],const real q[],const real alPol[],
		    int **neighList,const int neighPair[],real dBuffer[])
{
  char UPLO='U',DIAG='N';
  
  int i,ii,j,jj,k,l;
  int mix,miy,miz;
  int mjx,mjy,mjz;
  int m1,m2,m3;
  int ierr;
  
  real alPolInv,epol,virpol;
  real r2,rt,rt2,rt3,rt5;
  real tm1,tm2,tm3,tm4,tm5,tm6;
  real delta[3];
  
  static_field(param,parallel,box,x,y,z,q,neighList,neighPair,dBuffer);
  
  for(i=0;i<polar->nPolTensor;i++)
  {
    polTensor[i]=0.;
  }
  
  for(i=parallel->idProc; i<param->nAtom; i+=parallel->nProc)
  {
    if(!polList[i])
      continue;
    
    ii=polMap[i];
    
    mix=3*ii;
    miy=mix+1;
    miz=mix+2;
    
    alPolInv=1./alPol[i];
    
    m1=mix+mix*(mix+1)/2;
    polTensor[m1]=alPolInv;
    
    m2=miy+miy*(miy+1)/2;
    polTensor[m2]=alPolInv;
    
    m3=miz+miz*(miz+1)/2;
    polTensor[m3]=alPolInv;
    
  }
  
  l=0;
  for(i=parallel->idProc; i<param->nAtom; i+=parallel->nProc)
  {
    if(!polList[i])
      continue;
    
    ii=polMap[i];
    
    for(k=0; k<neighPair[l]; k++)
    {
      j=neighList[l][k];
      
      if(!polList[j])
	continue;
      
      jj=polMap[j];

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
	
	if(i<j)
	{
	  mix=3*ii;
	  mjx=3*jj;
	}
	else
	{
	  mix=3*jj;
	  mjx=3*ii;
	}
	
	miy=mix+1;
	miz=mix+2;
	
	mjy=mjx+1;
	mjz=mjx+2;
	
	m1=mix+mjx*(mjx+1);
	polTensor[m1]=-tm1;
	
	m1++; // m1=miy+mjx*(mjx+1)
	polTensor[m1]=-tm2;
	
	m1++; // m1=miz+mjx*(mjx+1)
	polTensor[m1]=-tm4;
	
	m2=mix+mjy*(mjy+1);
	polTensor[m2]=-tm2;
	
	m2++; // m2=miy+mjy*(mjy+1)
	polTensor[m2]=-tm3;
	
	m2++; // m2=miz+mjy*(mjy+1)
	polTensor[m2]=-tm5;
	
	m3=mix+mjz*(mjz+1);
	polTensor[m3]=-tm4;
	
	m3++; // m3=miy+mjz*(mjz+1)
	polTensor[m3]=-tm5;
	
	m3++; // m3=miz+mjz*(mjz+1)
	polTensor[m3]=-tm6;
	
      } // end if(r2<=param->cutOff2)
      
    } // end for(k=0; k<neighPair[l]; k++)
    
    l++;
    
  } // end for(i=parallel->idProc; i<param->nAtom; i+=parallel->nProc)
  
  if(parallel->nProc>1)
  {
    sum_double_para(polTensor,dBuffer,polar->nPolTensor);
  }
  
  dtptri_(&UPLO,&DIAG,&(polar->nAtPol),polTensor,&ierr);
  
  if(ierr>0)
    my_error(POLAR_DIAG_ERROR,__FILE__,__LINE__,0);
  else if(ierr<0)
    my_error(POLAR_VALU_ERROR,__FILE__,__LINE__,0);
  
  for(i=0;i<param->nAtom;i++)
  {
    muIndX[i]=0.;
    muIndY[i]=0.;
    muIndZ[i]=0.;
  }
  
  l=0;
  for(i=parallel->idProc; i<param->nAtom; i+=parallel->nProc)
  {
    if(!polList[i])
      continue;
    
    ii=polMap[i];
    
    mix=3*ii;
    miy=mix+1;
    miz=mix+2;
    
    m1=mix+mix*(mix+1);
    m2=miy+miy*(miy+1);
    m3=miz+miz*(miz+1);
    
    muIndX[i]+=polTensor[m1]*elFieldX[i];
    muIndY[i]+=polTensor[m2]*elFieldY[i];
    muIndZ[i]+=polTensor[m3]*elFieldZ[i];
    
    for(k=0; k<neighPair[l]; k++)
    {
      j=neighList[l][k];
      
      if(!polList[j])
	continue;
      
      jj=polMap[j];
	
      if(i<j)
      {
	mix=3*ii;
	mjx=3*jj;
      }
      else
      {
	mix=3*jj;
	mjx=3*ii;
      }
	
      miy=mix+1;
      miz=mix+2;
      
      mjy=mjx+1;
      mjz=mjx+2;
      
      m1=mix+mjx*(mjx+1);
      m2=mix+mjy*(mjy+1);
      m3=mix+mjz*(mjz+1);
      
      muIndX[j]+=polTensor[m1]*elFieldX[i]+polTensor[m2]*elFieldY[i]+polTensor[m3]*elFieldZ[i];
      muIndX[i]+=polTensor[m1]*elFieldX[j]+polTensor[m2]*elFieldY[j]+polTensor[m3]*elFieldZ[j];
      
      m1++;
      m2++;
      m3++;
      
      muIndY[j]+=polTensor[m1]*elFieldX[i]+polTensor[m2]*elFieldY[i]+polTensor[m3]*elFieldZ[i];
      muIndY[i]+=polTensor[m1]*elFieldX[j]+polTensor[m2]*elFieldY[j]+polTensor[m3]*elFieldZ[j];
      
      m1++;
      m2++;
      m3++;
      
      muIndZ[j]+=polTensor[m1]*elFieldX[i]+polTensor[m2]*elFieldY[i]+polTensor[m3]*elFieldZ[i];
      muIndZ[i]+=polTensor[m1]*elFieldX[j]+polTensor[m2]*elFieldY[j]+polTensor[m3]*elFieldZ[j];
      
    } // end for(k=0; k<neighPair[l]; k++)
    
    l++;
    
  } // end for(i=parallel->idProc; i<param->nAtom; i+=parallel->nProc)
  
  if(parallel->nProc>1)
  {
    sum_double_para(muIndX,dBuffer,param->nAtom);
    sum_double_para(muIndY,dBuffer,param->nAtom);
    sum_double_para(muIndZ,dBuffer,param->nAtom);
  }
  
  polar_forces(param,parallel,box,x,y,z,fx,fy,fz,q,virpol,neighList,neighPair);
  
  for(i=parallel->idProc; i<param->nAtom; i+=parallel->nProc)
  {
    if(!polList[i])
      continue;
    
    epol+=muIndX[i]*elFieldX[i]+muIndY[i]*elFieldY[i]+muIndZ[i]*elFieldZ[i];
    
  }
  
  epol*=0.5*param->chargeConst;
  
  ener->epol=epol;
  ener->virpol=virpol;
  
}

void polar_forces(PARAM *param,PARALLEL *parallel,PBC *box,const real x[],
		  const real y[],const real z[],real fx[],real fy[],
		  real fz[],const real q[],real *virpol,int **neighList,
		  const int neighPair[])
{
  int i,j,k,l;
  
  real qi,qj;
  real fxi,fyi,fzi;
  real r2,rt,rt2,rt3,rt5;
  real tm1,tm2,tm3,tm4,tm5,tm6;
  real tmxi,tmyi,tmzi,tmxj,tmyj,tmzj;
  real mrt2,mrt5,mufi,mufj,mufij;
  real muftx,mufty,muftz;
  real fpolx,fpoly,fpolz;
  real delta[3];
  
  *virpol=0.;

  l=0;
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
	
	fpolx=qj*tmxi-qi*tmxj+muftx;
	fpoly=qj*tmyi-qi*tmyj+mufty;
	fpolz=qj*tmzi-qi*tmzj+muftz;
	
	fx[j]+=fpolx;
	fy[j]+=fpoly;
	fz[j]+=fpolz;
	
	fxi+=-fpolx;
	fyi+=-fpoly;
	fzi+=-fpolz;
	
	*virpol-=( (fpolx*delta[0]) + (fpoly*delta[1]) + (fpolz*delta[2]) );
      }
      
    }
    
    fx[i]+=fxi;
    fy[i]+=fyi;
    fz[i]+=fzi;
    
    l++;
  }
  
}