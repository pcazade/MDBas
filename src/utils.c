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
 * \file utils.c
 * \brief Contains various and general utilitary functions.
 * \author Pierre-Andre Cazade and Florent Hedin
 * \version alpha-branch
 * \date 2012
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <float.h>
#include <complex.h>

#include "global.h"
#include "utils.h"
#include "io.h"
#include "errors.h"
#include "memory.h"
#include "parallel.h"

#ifdef _OPENMP
#undef _OPENMP
#endif

/**
 * \param delta Vector of length 3, representing the vector ij.
 * \param simulCond Pointer to structure SIMULPARAMS containing parameters of the current simulation.
 * \param box Pointer to structure PBC containing Periodic Boundaries Conditions parameters.
 *
 * \brief Evaluates the distance between two atoms i and j.
 * \remarks Periodic Boundaries Conditions are internally applied.
 *
 * \return The distance in Angstroems between two atoms.
 */

double dist(const PBC *box, double delta[3])
{
  double r2,xt,yt,zt;
  
  switch(box->type)
  {
    case CUBIC:
      delta[0]-=box->a1*nint(delta[0]*box->u1);
      delta[1]-=box->a1*nint(delta[1]*box->u1);
      delta[2]-=box->a1*nint(delta[2]*box->u1);
    break;
    
    case ORBIC:
      delta[0]-=box->a1*nint(delta[0]*box->u1);
      delta[1]-=box->b2*nint(delta[1]*box->v2);
      delta[2]-=box->c3*nint(delta[2]*box->w3);
    break;
    
    case TCLIN:
      xt=delta[0]*box->u1+delta[1]*box->u2+delta[2]*box->u3;
      yt=delta[0]*box->v1+delta[1]*box->v2+delta[2]*box->v3;
      zt=delta[0]*box->w1+delta[1]*box->w2+delta[2]*box->w3;
      
      xt-=nint(xt);
      yt-=nint(yt);
      zt-=nint(zt);
      
      delta[0]=xt*box->a1+yt*box->b1+zt*box->c1;
      delta[1]=xt*box->a2+yt*box->b2+zt*box->c2;
      delta[2]=xt*box->a3+yt*box->b3+zt*box->c3;
    break;
    
    default:
    break;
  }
  
  //r=sqrt(X2(delta[0])+X2(delta[1])+X2(delta[2]));
  
  r2=(X2(delta[0])+X2(delta[1])+X2(delta[2]));
  
  return r2;
  
}

double dist2(const PBC *box, double *dx, double *dy,double *dz)
{
  double r2,xt,yt,zt;
  
  switch(box->type)
  {
    case CUBIC:
      *dx-=box->a1*nint(*dx*box->u1);
      *dy-=box->a1*nint(*dy*box->u1);
      *dz-=box->a1*nint(*dz*box->u1);
      break;
      
    case ORBIC:
      *dx-=box->a1*nint(*dx*box->u1);
      *dy-=box->b2*nint(*dy*box->v2);
      *dz-=box->c3*nint(*dz*box->w3);
    break;
      
    case TCLIN:
      xt=*dx*box->u1+*dy*box->u2+*dz*box->u3;
      yt=*dx*box->v1+*dy*box->v2+*dz*box->v3;
      zt=*dx*box->w1+*dy*box->w2+*dz*box->w3;
      
      xt-=nint(xt);
      yt-=nint(yt);
      zt-=nint(zt);
      
      *dx=xt*box->a1+yt*box->b1+zt*box->c1;
      *dy=xt*box->a2+yt*box->b2+zt*box->c2;
      *dz=xt*box->a3+yt*box->b3+zt*box->c3;
    break;
      
    default:
      break;
  }
  
  r2=( X2(*dx)+X2(*dy)+X2(*dz) );
  
  return r2;
  
}

/**
 * \param param Pointer to structure PARAM containing parameters of the current simulation.
 * \param atom Array of structure ATOM (coordinates, forces, etc...).
 *
 * \brief Set velocities and forces of frozen atoms to zero.
 * \remarks Called by energy if frozen atoms are used.
 */
void freeze_atoms(const PARAM *param, double *vx, double *vy,double *vz,
		  double *fx,double *fy,double *fz,const int *frozen)
{
  int i;
  
  for(i=0;i<param->nAtom;i++)
  {
    
    if(frozen[i])
    {
      vx[i]=0.;
      vy[i]=0.;
      vz[i]=0.;
      
      fx[i]=0.;
      fy[i]=0.;
      fz[i]=0.;
      
    }
    
  }
  
}

/**
 * \param atom Array of structure ATOM (coordinates, forces, etc...).
 * \param simulCond Pointer to structure SIMULPARAMS containing parameters of the current simulation.
 * \param box Pointer to structure PBC containing Periodic Boundaries Conditions parameters.
 *
 * \brief Applies Periodic Boundaries Conditions to atom[], according to the PBC type.
 */
void image_update(const PARALLEL *parallel,const PBC *box,double x[],double y[], double z[])
{

  int i;
  double xt,yt,zt;
  
  switch(box->type)
  {
    case CUBIC:
      #ifdef _OPENMP
      #pragma omp parallel for default(none) shared(param,atom,box) private(i)
      #endif
      for(i=parallel->fAtProc;i<parallel->lAtProc;i++)
      {
	x[i]=x[i]-box->a1*nint(x[i]*box->u1);
	y[i]=y[i]-box->a1*nint(y[i]*box->u1);
	z[i]=z[i]-box->a1*nint(z[i]*box->u1);
      }
    break;
      
    case ORBIC:
      #ifdef _OPENMP
      #pragma omp parallel for default(none) shared(param,atom,box) private(i)
      #endif
      for(i=parallel->fAtProc;i<parallel->lAtProc;i++)
      {
	x[i]=x[i]-box->a1*nint(x[i]*box->u1);
	y[i]=y[i]-box->b2*nint(y[i]*box->v2);
	z[i]=z[i]-box->c3*nint(z[i]*box->w3);
      }
    break;
      
    case TCLIN:
      #ifdef _OPENMP
      #pragma omp parallel for default(none) shared(param,atom,box) private(i,xt,yt,zt)
      #endif
      for(i=parallel->fAtProc;i<parallel->lAtProc;i++)
      {
	xt=x[i]*box->u1+y[i]*box->u2+z[i]*box->u3;
	yt=x[i]*box->v1+y[i]*box->v2+z[i]*box->v3;
	zt=x[i]*box->w1+y[i]*box->w2+z[i]*box->w3;
	
	xt-=nint(xt);
	yt-=nint(yt);
	zt-=nint(zt);
	
	x[i]=xt*box->a1+yt*box->b1+zt*box->c1;
	y[i]=xt*box->a2+yt*box->b2+zt*box->c2;
	z[i]=xt*box->a3+yt*box->b3+zt*box->c3;
      }
    break;
      
    default:
    break;
  }
  
}

/**
 * \param size_array Size of the array d[].
 * \param d Array of structure DELTA.
 * \param simulCond Pointer to structure SIMULPARAMS containing parameters of the current simulation.
 * \param box Pointer to structure PBC containing Periodic Boundaries Conditions parameters.
 *
 * \brief Applies Periodic Boundaries Conditions to d[], according to the PBC type.
 */
void image_array(const PBC *box, double dx[],double dy[],double dz[], const int size_array)
{

  int i;
  double xt,yt,zt;
  
  switch(box->type)
  {
    case CUBIC:
      #ifdef _OPENMP
      #pragma omp parallel for default(none) shared(size_array,box,d) private(i)
      #endif
      for(i=0;i<size_array;i++)
      {
	dx[i]-=box->a1*nint(dx[i]*box->u1);
	dy[i]-=box->a1*nint(dy[i]*box->u1);
	dz[i]-=box->a1*nint(dz[i]*box->u1);
      }
    break;
      
    case ORBIC:
      #ifdef _OPENMP
      #pragma omp parallel for default(none) shared(size_array,box,d) private(i)
      #endif
      for(i=0;i<size_array;i++)
      {
	dx[i]-=box->a1*nint(dx[i]*box->u1);
	dy[i]-=box->b2*nint(dy[i]*box->v2);
	dz[i]-=box->c3*nint(dz[i]*box->w3);
      }
    break;
      
    case TCLIN:
      #ifdef _OPENMP
      #pragma omp parallel for default(none) shared(size_array,box,d) private(i,xt,yt,zt)
      #endif
      for(i=0;i<size_array;i++)
      {
	xt=dx[i]*box->u1+dy[i]*box->u2+dz[i]*box->u3;
	yt=dx[i]*box->v1+dy[i]*box->v2+dz[i]*box->v3;
	zt=dx[i]*box->w1+dy[i]*box->w2+dz[i]*box->w3;
	
	xt-=nint(xt);
	yt-=nint(yt);
	zt-=nint(zt);
	
	dx[i]=xt*box->a1+yt*box->b1+zt*box->c1;
	dy[i]=xt*box->a2+yt*box->b2+zt*box->c2;
	dz[i]=xt*box->a3+yt*box->b3+zt*box->c3;
      }
    break;
      
    default:
    break;
  }
  
}

/**
 * \param box Pointer to structure PBC containing Periodic Boundaries Conditions parameters.
 *
 * \brief Computes norms, Wigner-Seitz cell, determinant, volume for a given box.
 */
void scale_box(PBC *box, const double cell0[9], const double scale)
{
  
  box->a1=scale*cell0[0];
  box->a2=scale*cell0[1];
  box->a3=scale*cell0[2];
  box->b1=scale*cell0[3];
  box->b2=scale*cell0[4];
  box->b3=scale*cell0[5];
  box->c1=scale*cell0[6];
  box->c2=scale*cell0[7];
  box->c3=scale*cell0[8];
  
  // get norm
  box->a=sqrt(X2(box->a1)+X2(box->a2)+X2(box->a3));
  box->b=sqrt(X2(box->b1)+X2(box->b2)+X2(box->b3));
  box->c=sqrt(X2(box->c1)+X2(box->c2)+X2(box->c3));
  
  // computes Wigner-Seitz cell
  box->u1=box->b2*box->c3-box->c2*box->b3;
  box->u2=box->c1*box->b3-box->b1*box->c3;
  box->u3=box->b1*box->c2-box->c1*box->b2;
  
  box->v1=box->c2*box->a3-box->a2*box->c3;
  box->v2=box->a1*box->c3-box->c1*box->a3;
  box->v3=box->c1*box->a2-box->a1*box->c2;
  
  box->w1=box->a2*box->b3-box->b2*box->a3;
  box->w2=box->b1*box->a3-box->a1*box->b3;
  box->w3=box->a1*box->b2-box->b1*box->a2;
  
  // determinant
  box->det=box->a1*box->u1+box->a2*box->u2+box->a3*box->u3;
  // volume
  box->vol=fabs(box->det);
  
  /* obtain parameters of a virutal orthorombic cell containing
   * the triclinic cell
   */
  box->pa=box->vol/sqrt(X2(box->u1)+X2(box->u2)+X2(box->u3));
  box->pb=box->vol/sqrt(X2(box->v1)+X2(box->v2)+X2(box->v3));
  box->pc=box->vol/sqrt(X2(box->w1)+X2(box->w2)+X2(box->w3));
  
  box->u=1.0/box->pa;
  box->v=1.0/box->pb;
  box->w=1.0/box->pc;
  
  // normalise Wigner-Seitz vectors
  if(fabs(box->det)>0.)
  {
    box->u1/=box->det;
    box->u2/=box->det;
    box->u3/=box->det;
    
    box->v1/=box->det;
    box->v2/=box->det;
    box->v3/=box->det;
    
    box->w1/=box->det;
    box->w2/=box->det;
    box->w3/=box->det;
  }
  else
  {
    // avoid div by 0
    box->u1=0.0;
    box->v1=0.0;
    box->w1=0.0;
    
    box->u2=0.0;
    box->v2=0.0;
    box->w2=0.0;
    
    box->u3=0.0;
    box->v3=0.0;
    box->w3=0.0;
  }
  
}

/**
 * \param box Pointer to structure PBC containing Periodic Boundaries Conditions parameters.
 *
 * \brief Computes norms, Wigner-Seitz cell, determinant, volume for a given box.
 */
void vv_scale_box(PBC *box, const double scale)
{
  
  box->a1*=scale;
  box->a2*=scale;
  box->a3*=scale;
  box->b1*=scale;
  box->b2*=scale;
  box->b3*=scale;
  box->c1*=scale;
  box->c2*=scale;
  box->c3*=scale;
  
  // get norm
  box->a=sqrt(X2(box->a1)+X2(box->a2)+X2(box->a3));
  box->b=sqrt(X2(box->b1)+X2(box->b2)+X2(box->b3));
  box->c=sqrt(X2(box->c1)+X2(box->c2)+X2(box->c3));
  
  // computes Wigner-Seitz cell
  box->u1=box->b2*box->c3-box->c2*box->b3;
  box->u2=box->c1*box->b3-box->b1*box->c3;
  box->u3=box->b1*box->c2-box->c1*box->b2;
  
  box->v1=box->c2*box->a3-box->a2*box->c3;
  box->v2=box->a1*box->c3-box->c1*box->a3;
  box->v3=box->c1*box->a2-box->a1*box->c2;
  
  box->w1=box->a2*box->b3-box->b2*box->a3;
  box->w2=box->b1*box->a3-box->a1*box->b3;
  box->w3=box->a1*box->b2-box->b1*box->a2;
  
  // determinant
  box->det=box->a1*box->u1+box->a2*box->u2+box->a3*box->u3;
  // volume
  box->vol=fabs(box->det);
  
  /* obtain parameters of a virutal orthorombic cell containing
   * the triclinic cell
   */
  box->pa=box->vol/sqrt(X2(box->u1)+X2(box->u2)+X2(box->u3));
  box->pb=box->vol/sqrt(X2(box->v1)+X2(box->v2)+X2(box->v3));
  box->pc=box->vol/sqrt(X2(box->w1)+X2(box->w2)+X2(box->w3));
  
  box->u=1.0/box->pa;
  box->v=1.0/box->pb;
  box->w=1.0/box->pc;
  
  // normalise Wigner-Seitz vectors
  if(fabs(box->det)>0.)
  {
    box->u1/=box->det;
    box->u2/=box->det;
    box->u3/=box->det;
    
    box->v1/=box->det;
    box->v2/=box->det;
    box->v3/=box->det;
    
    box->w1/=box->det;
    box->w2/=box->det;
    box->w3/=box->det;
  }
  else
  {
    // avoid div by 0
    box->u1=0.0;
    box->v1=0.0;
    box->w1=0.0;
    
    box->u2=0.0;
    box->v2=0.0;
    box->w2=0.0;
    
    box->u3=0.0;
    box->v3=0.0;
    box->w3=0.0;
  }
  
}

void box_to_lattice(const PBC *box, double lattice[6])
{
  double cost;
  
  lattice[0]=box->a;
  lattice[1]=box->b;
  lattice[2]=box->c;
  
  cost = ( (box->b1*box->c1) + (box->b2*box->c2) + (box->b3*box->c3) ) / ( box->b * box->c );
  lattice[3]= acos(cost)*180./PI;
  
  cost = ( (box->a1*box->c1) + (box->a2*box->c2) + (box->a3*box->c3) ) / ( box->a * box->c );
  lattice[4]= acos(cost)*180./PI;
  
  cost = ( (box->a1*box->b1) + (box->a2*box->b2) + (box->a3*box->b3) ) / ( box->a * box->b );
  lattice[5]= acos(cost)*180./PI;
  
}

void box_to_crystal(const PBC *box,double crystal[6])
{
  
  /** crystal is the lower triangle of the symetric matrix of the box vectors
   * a1
   * a2  b2
   * a3  b3  c3
   **/
  
  crystal[0]=box->a1;
  crystal[2]=box->b2;
  crystal[5]=box->c3;
  
  crystal[1]=box->a2;
  crystal[3]=box->a3;
  crystal[4]=box->b3;
  
}

/**
 * \param atom Array of structure ATOM (coordinates, forces, etc...).
 * \param simulCond Pointer to structure SIMULPARAMS containing parameters of the current simulation.
 *
 * \brief Estimates the kinetic energy of a system.
 *
 * \return On return, the kinetic energy : \f$ 1/2*m*v^2 \f$
 */
double kinetic(const PARALLEL *parallel,const double vx[],const double vy[],
	       const double vz[],const double mass[],double *dBuffer)
{
  int i;
  double ekin;
  
  ekin=0.;
  for(i=parallel->fAtProc;i<parallel->lAtProc;i++)
  {
    ekin+=mass[i]*(X2(vx[i])+X2(vy[i])+X2(vz[i]));
  }
  
  if(parallel->nProc>1)
    sum_double_para(&ekin,dBuffer,1);
  
  return ( ekin*0.5 );
}

/**
 * \param atom Array of structure ATOM (coordinates, forces, etc...).
 * \param simulCond Pointer to structure SIMULPARAMS containing parameters of the current simulation.
 * \param stress Array representing (parts of) a Stress Tensor.
 *
 * \brief Adds to the stress tensor the kinetic energy.
 */
void stress_kinetic(const PARALLEL *parallel,const double vx[],const double vy[],
		    const double vz[],const double mass[],double stress[6],
		    double dBuffer[])
{
  
  for(int i=0;i<6;i++)
    stress[i]=0.;
  
  for(int i=parallel->fAtProc;i<parallel->lAtProc;i++)
  {
    stress[0] += mass[i]*X2(vx[i]);
    stress[1] += mass[i]*vx[i]*vy[i];
    stress[2] += mass[i]*vx[i]*vz[i];
    stress[3] += mass[i]*X2(vy[i]);
    stress[4] += mass[i]*vy[i]*vz[i];
    stress[5] += mass[i]*X2(vz[i]);
  }
  
  if(parallel->nProc>1)
    sum_double_para(stress,dBuffer,6);
}

/**
 * \param atom Array of structure ATOM (coordinates, forces, etc...).
 * \param simulCond Pointer to structure SIMULPARAMS containing parameters of the current simulation.
 * \param box Pointer to structure PBC containing Periodic Boundaries Conditions parameters.
 *
 * \brief Get initial kinetic energy from the temperature.
 */
void get_kinfromtemp(PARAM *param, const PBC *box)
{
  double degf;
  
  get_degfree(param,box);
  
  //   Energy in internal units 10 J/mol. rboltzui=R/10.
  
  degf=(double)param->nDegFree;
  param->kinTemp0=0.5*param->temp0*degf*rboltzui;
}

/**
 * \param simulCond Pointer to structure SIMULPARAMS containing parameters of the current simulation.
 * \param box Pointer to structure PBC containing Periodic Boundaries Conditions parameters.
 *
 * \brief Obtains the number of degrees of freedom of the system.
 */
void get_degfree(PARAM *param, const PBC *box)
{
  
//   Atoms degrees of freedom - 3 degrees of freedom of the CoM.
  param->nDegFree=3*(param->nAtom-param->nFrozen)-3;
  
//   -3 degrees of freedpm for the box rotation when no PBC.
  if(box->type==NOBOX)
    param->nDegFree-=3;
  
  if(param->nConst>0)
    param->nDegFree-=param->nConst;
    
}
 
/**
 * \param str A character string of any length.
 *
 * \brief Sets to lower case all the characters of the string str.
 */
void nocase(char *str)
{
  int i,n;
  n=strlen(str);
  for(i=0;i<n;i++)
  {
    str[i]=tolower(str[i]);
  }
}

/**
 * \param x Any double precision number.
 *
 * \brief Rounds its argument to the nearest whole number.
 *
 * \return x with the fractional portion of its magnitude eliminated by rounding to the nearest whole number and with its sign preserved, casted to int.
 */
int nint(const double x)
{
  return ( (x)>=0.?(int)((x)+0.5):(int)((x)-0.5) );
}


/**
 * \param buf Any complex numbers vector of size n a power of 2.
 *
 * \brief Piece of code fro Fast Fourier Transform (FFT) from Rosetta code project: http://rosettacode.org.
 *
 */
 
// void fft2(cplx buf[], cplx out[], int n, int step)
// {
// 	if (step < n) {
// 		fft2(out, buf, n, step * 2);
// 		fft2(out + step, buf + step, n, step * 2);
//  
// 		for (int i = 0; i < n; i += 2 * step) {
// 			cplx t = cexp(-I * PI * i / n) * out[i + step];
// 			buf[i / 2]     = out[i] + t;
// 			buf[(i + n)/2] = out[i] - t;
// 		}
// 	}
// }
//  
// void fft1(cplx buf[], int n)
// {
// 	cplx out[n];
// 	for (int i = 0; i < n; i++) out[i] = buf[i];
//  
// 	fft2(buf, out, n, 1);
// }