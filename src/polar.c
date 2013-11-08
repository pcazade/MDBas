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

void init_polar(CTRL *ctrl,PARAM *param,PARALLEL *parallel,PBC *box)
{
  
}

void static_field(CTRL *ctrl,PARAM *param,PARALLEL *parallel,PBC *box,
		  const double x[],const double y[],const double z[],
		  const double q[],double elFieldX[],double elFieldY[],
		  double elFieldZ[],int **neighList,const int neighPair[],
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
      rt2=1./r2;
      rt=sqrt(rt2);
      rt3=rt*rt2;

      if(r2<=param->cutOff2)
      {
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

void induced_filed(CTRL *ctrl,PARAM *param,PARALLEL *parallel,PBC *box,
		  const double x[],const double y[],const double z[],
		  const double q[],double elFieldX[],double elFieldY[],
		  double elFieldZ[],int **neighList,const int neighPair[],
		  double dBuffer[])
{
  
}

double polar_ener()
{
  
}