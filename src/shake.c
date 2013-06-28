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
 * \file shake.c
 * \brief Contains functions for applying the SHAKE constraints.
 * \author Pierre-Andre Cazade and Florent Hedin
 * \version alpha-branch
 * \date 2012
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "global.h"
#include "utils.h"
#include "io.h"
#include "errors.h"
#include "memory.h"
#include "parallel.h"

#if (defined TIMING && defined __unix__ && !defined __STRICT_ANSI__)
#define TIMER
#include "timing.h"
#endif

static double *xt,*yt,*zt;
static double *rt2;
static double *dtx,*dty,*dtz;

void shake_allocate_arrays(const PARAM *param,const PARALLEL *parallel)
{
  xt=yt=zt=rt2=dtx=dty=dtz=NULL;

  xt=(double*)my_malloc(param->nAtom*sizeof(*xt));
  yt=(double*)my_malloc(param->nAtom*sizeof(*yt));
  zt=(double*)my_malloc(param->nAtom*sizeof(*zt));

  rt2=(double*)my_malloc(param->nConst*sizeof(*rt2));

  dtx=(double*)my_malloc(param->nConst*sizeof(*dtx));
  dty=(double*)my_malloc(param->nConst*sizeof(*dty));
  dtz=(double*)my_malloc(param->nConst*sizeof(*dtz));
}

void shake_free_arrays()
{
  free(xt);     free(yt);       free(zt);
  free(rt2);
  free(dtx);    free(dty);      free(dtz);
}

void lf_shake(PARAM *param,PBC *box,CONSTRAINT constList[],PARALLEL *parallel,
	      double x[],double y[],double z[],
	      double ddx[],double ddy[],double ddz[],double rmass[],
	      int *nAtConst,double stress[6],double *virshake)
{
  int i,l,ia,ib,icycle,converged;
  double ts2,maxdist,dist;
  double lambda,lambdai,lambdaj,t2rmi,t2rmj,nia,nib;
  
#ifdef TIMER
  update_timer_begin(TIMER_SHAKE,__func__);
#endif
  
  icycle=0;
  converged=0;
  
  *virshake=0.;
  for(i=0;i<6;i++)
  {
    stress[i]=0.;
  }
  
  ts2=X2(param->timeStep);

  do
  {
    maxdist=0.;
    
    l=0;
    for(i=parallel->fConst;i<parallel->lConst;i++)
    {
      ia=constList[i].a;
      ib=constList[i].b;
      
      dtx[l]=x[ib]-x[ia];
      dty[l]=y[ib]-y[ia];
      dtz[l]=z[ib]-z[ia];
      
      rt2[l]=dist2(box,&(dtx[l]),&(dty[l]),&(dtz[l]));
      
      dist=fabs(rt2[l]-constList[i].rc2)/sqrt(constList[i].rc2);
      maxdist=MAX(maxdist,dist);
      
      l++;
    }
    
    maxdist=0.5*maxdist;
    
    if(maxdist<param->tolShake)
      converged=1;
    
    if(parallel->nProc>1)
    {
      sum_double_para(converged);
    }
    
    if(!converged)
    {
      for(i=0;i<param->nAtom;i++)
      {
	xt[i]=0.;
	yt[i]=0.;
	zt[i]=0.;
      }
      
      l=0;
      for(i=parallel->fConst;i<parallel->lConst;i++)
      {
	ia=constList[i].a;
	ib=constList[i].b;
	
	t2rmi=ts2*rmass[ia];
	t2rmj=ts2*rmass[ib];
	
	lambda=-(constList[i].rc2-rt2[l])/(2.*(t2rmi+t2rmj)*
	  ((ddx[l]*dtx[l])+(ddy[l]*dty[l])+(ddz[l]*dtz[l])));
	
	*virshake+=lambda*(X2(ddx[l])+X2(ddy[l])+X2(ddz[l]));
	
	stress[0]-=lambda*X2(ddx[l]);
	stress[1]-=lambda*ddx[l]*ddy[l];
	stress[2]-=lambda*ddx[l]*ddz[l];
	stress[3]-=lambda*X2(ddy[l]);
	stress[4]-=lambda*ddy[l]*ddz[l];
	stress[5]-=lambda*X2(ddz[l]);
	
	lambdai=lambda*t2rmi;
	xt[ia]+=ddx[l]*lambdai;
	yt[ia]+=ddy[l]*lambdai;
	zt[ia]+=ddz[l]*lambdai;
            
	lambdaj=-lambda*t2rmj;
	xt[ib]+=ddx[l]*lambdaj;
	yt[ib]+=ddy[l]*lambdaj;
	zt[ib]+=ddz[l]*lambdaj;
	
	l++;
	
      }
      
      if(parallel->nProc>1)
      {
	sum_double_para(xt,buffer,param->nAtom);
	sum_double_para(yt,buffer,param->nAtom);
	sum_double_para(zt,buffer,param->nAtom);
      }
      
//       for(i=parallel->fConst;i<parallel->lConst;i++)
//       {
// 	ia=constList[i].a;
// 	ib=constList[i].b;
// 	
// 	nia=1.0/(double)nAtConst[ia];
// 	nib=1.0/(double)nAtConst[ib];
// 	
// 	x[ia]+=xt[ia]*nia;
// 	y[ia]+=yt[ia]*nia;
// 	z[ia]+=zt[ia]*nia;
// 	
// 	x[ib]+=xt[ib]*nib;
// 	y[ib]+=yt[ib]*nib;
// 	z[ib]+=zt[ib]*nib;
// 	
//       }
      
      for(i=parallel->fAtom;i<parallel->lAtom;i++)
      {
	x[i]+=xt[i];
	y[i]+=yt[i];
	z[i]+=zt[i];
      }
      
    }
    
    icycle++;
    
  }while( (!converged) && (icycle<param->maxCycle) );
  
  if(parallel->nProc>1)
  {
    buffer[0]=virshake;
    sum_double_para(buffer,&(buffer[1]),1);
    virshake=buffer[0];
    
    sum_double_para(stress,buffer,6);
  }
  
  if(!converged)
    my_error(CONVERG_SHAKE_ERROR,__FILE__,__LINE__,0);
  
  box->stress1+=stress[0];
  box->stress2+=stress[1];
  box->stress3+=stress[2];
  box->stress4+=stress[1];
  box->stress5+=stress[3];
  box->stress6+=stress[4];
  box->stress7+=stress[2];
  box->stress8+=stress[4];
  box->stress9+=stress[5];
  
  #ifdef TIMER
  update_timer_end(TIMER_SHAKE,__func__);
  #endif
  
}


void vv_shake_r(PARAM *param,PBC *box,CONSTRAINT constList[],PARALLEL *parallel,
		double x[],double y[],double z[],
		double vx[],double vy[],double vz[],
		double ddx[],double ddy[],double ddz[],double rmass[],
		int *nAtConst,double stress[6],double *virshake)
{
  int i,ia,ib,icycle,converged;
  double maxdist,dist;
  double lambda,lambdai,lambdaj,trmi,trmj,nia,nib;
  
#ifdef TIMER
  update_timer_begin(TIMER_SHAKE,__func__);
#endif
  
  icycle=0;
  converged=0;
  
  *virshake=0.;
  for(i=0;i<6;i++)
  {
    stress[i]=0.;
  }

  do
  {
    maxdist=0.;
    
    for(i=0;i<param->nConst;i++)
    {
      ia=constList[i].a;
      ib=constList[i].b;
      
      dtx[i]=x[ib]-x[ia];
      dty[i]=y[ib]-y[ia];
      dtz[i]=z[ib]-z[ia];
      
      rt2[i]=dist2(box,&(dtx[i]),&(dty[i]),&(dtz[i]));
      
      dist=fabs(rt2[i]-constList[i].rc2)/sqrt(constList[i].rc2);
      maxdist=MAX(maxdist,dist);
    }
    
    maxdist=0.5*maxdist;
    
    if(maxdist<param->tolShake)
      converged=1;
    
    if(!converged)
    {
      for(i=0;i<param->nAtom;i++)
      {
	xt[i]=0.;
	yt[i]=0.;
	zt[i]=0.;
      }
      
      for(i=0;i<param->nConst;i++)
      {
	ia=constList[i].a;
	ib=constList[i].b;
	
	trmi=param->timeStep*rmass[ia];
	trmj=param->timeStep*rmass[ib];
	
	lambda=-(constList[i].rc2-rt2[i])/(param->timeStep*(trmi+trmj)*
	  ((ddx[i]*dtx[i])+(ddy[i]*dty[i])+(ddz[i]*dtz[i])));
	
	*virshake+=lambda*(X2(ddx[i])+X2(ddy[i])+X2(ddz[i]));
	
	stress[0]-=lambda*X2(ddx[i]);
	stress[1]-=lambda*ddx[i]*ddy[i];
	stress[2]-=lambda*ddx[i]*ddz[i];
	stress[3]-=lambda*X2(ddy[i]);
	stress[4]-=lambda*ddy[i]*ddz[i];
	stress[5]-=lambda*X2(ddz[i]);
	
	lambdai=0.5*lambda*trmi;
	xt[ia]+=ddx[i]*lambdai;
	yt[ia]+=ddy[i]*lambdai;
	zt[ia]+=ddz[i]*lambdai;
            
	lambdaj=-0.5*lambda*trmj;
	xt[ib]+=ddx[i]*lambdaj;
	yt[ib]+=ddy[i]*lambdaj;
	zt[ib]+=ddz[i]*lambdaj;
	
      }
      
      for(i=0;i<param->nConst;i++)
      {
	ia=constList[i].a;
	ib=constList[i].b;
	
	nia=1.0/(double)nAtConst[ia];
	nib=1.0/(double)nAtConst[ib];
	
	x[ia]+=param->timeStep*xt[ia]*nia;
	y[ia]+=param->timeStep*yt[ia]*nia;
	z[ia]+=param->timeStep*zt[ia]*nia;
	
	x[ib]+=param->timeStep*xt[ib]*nib;
	y[ib]+=param->timeStep*yt[ib]*nib;
	z[ib]+=param->timeStep*zt[ib]*nib;
	
	vx[ia]+=xt[ia]*nia;
	vy[ia]+=yt[ia]*nia;
	vz[ia]+=zt[ia]*nia;
	
	vx[ib]+=xt[ib]*nib;
	vy[ib]+=yt[ib]*nib;
	vz[ib]+=zt[ib]*nib;
	
      }
      
    }
    
    icycle++;
    
  }while( (!converged) && (icycle<param->maxCycle) );
  
  if(!converged)
    my_error(CONVERG_SHAKE_ERROR,__FILE__,__LINE__,0);
    
#ifdef TIMER
  update_timer_end(TIMER_SHAKE,__func__);
#endif

}

void vv_shake_v(PARAM *param,CONSTRAINT constList[],PARALLEL *parallel,
		double vx[],double vy[],double vz[],double ddx[],
		double ddy[],double ddz[],double rmass[],int *nAtConst)
{
  int i,ia,ib,icycle,converged;
  double maxdist,tolvel;
  double lambda,lambdai,lambdaj,trmi,trmj,nia,nib;
  
#ifdef TIMER
  update_timer_begin(TIMER_SHAKE,__func__);
#endif
  
  icycle=0;
  converged=0;
  
  tolvel=param->tolShake*param->rTimeStep;

  do
  {
    
    for(i=0;i<param->nAtom;i++)
    {
      xt[i]=0.;
      yt[i]=0.;
      zt[i]=0.;
    }
    
    maxdist=0.;
    
    for(i=0;i<param->nConst;i++)
    {
      ia=constList[i].a;
      ib=constList[i].b;
      
      trmi=0.5*param->timeStep*rmass[ia];
      trmj=0.5*param->timeStep*rmass[ib];
      
      lambda=(ddx[i]*(vx[ib]-vx[ia])+ddy[i]*(vy[ib]-vy[ia])+
	ddz[i]*(vz[ib]-vz[ia]))/((trmi+trmj)*
	(X2(ddx[i])+X2(ddy[i])+X2(ddz[i])));
	
      maxdist=MAX(maxdist,fabs(lambda));
      
      lambdai=lambda*trmi;
      xt[ia]+=ddx[i]*lambdai;
      yt[ia]+=ddy[i]*lambdai;
      zt[ia]+=ddz[i]*lambdai;
	  
      lambdaj=-lambda*trmj;
      xt[ib]+=ddx[i]*lambdaj;
      yt[ib]+=ddy[i]*lambdaj;
      zt[ib]+=ddz[i]*lambdaj;
      
    }
    
    if(maxdist<tolvel)
      converged=1;
    
    if(!converged)
    {
      for(i=0;i<param->nConst;i++)
      {
	ia=constList[i].a;
	ib=constList[i].b;
	
	nia=1.0/(double)nAtConst[ia];
	nib=1.0/(double)nAtConst[ib];
	
	vx[ia]+=xt[ia]*nia;
	vy[ia]+=yt[ia]*nia;
	vz[ia]+=zt[ia]*nia;
	
	vx[ib]+=xt[ib]*nib;
	vy[ib]+=yt[ib]*nib;
	vz[ib]+=zt[ib]*nib;
	
      }
    }
    
    icycle++;
    
  }while( (!converged) && (icycle<param->maxCycle) );
  
  if(!converged)
    my_error(CONVERG_SHAKE_ERROR,__FILE__,__LINE__,0);
    
#ifdef TIMER
  update_timer_end(TIMER_SHAKE,__func__);
#endif
  
}