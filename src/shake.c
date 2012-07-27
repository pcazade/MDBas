#include <stdlib.h>
#include <math.h>
#include "global.h"
#include "utils.h"
#include "io.h"

void lf_shake(ATOM *atom,SIMULPARAMS *simulCond,CONSTRAINT *constList,DELTA *dd)
{
  int i,ia,ib,icycle,converged;
  double *xt,*yt,*zt,*rt2,ts2,maxdist,dist;
  double lambda,lambdai,lambdaj,t2rmi,t2rmj,nia,nib;
  DELTA *dt;
  
  xt=(double*)malloc(atom->natom*sizeof(*xt));
  yt=(double*)malloc(atom->natom*sizeof(*yt));
  zt=(double*)malloc(atom->natom*sizeof(*zt));
  
  rt2=(double*)malloc(simulCond->nconst*sizeof(*rt2));
  
  dt=(DELTA*)malloc(simulCond->nconst*sizeof(*dt));
  
  icycle=0;
  converged=0;
  
  ts2=X2(simulCond->timeStep);

  do
  {
    maxdist=0.;
    
    for(i=0;i<simulCond->nconst;i++)
    {
      ia=constList[i].a;
      ib=constList[i].b;
      
      rt2[i]=distance2(ia,ib,atom,&(dt[i]),simulCond);
      
      dist=fabs(rt2[i]-constList[i].rc2)/sqrt(constList[i].rc2);
      maxdist=MAX(maxdist,dist);
    }
    
    maxdist=0.5*maxdist;
    
    if(maxdist<simulCond->tolshake)
      converged=1;
    
    if(!converged)
    {
      for(i=0;i<atom->natom;i++)
      {
	xt[i]=0.;
	yt[i]=0.;
	zt[i]=0.;
      }
      
      for(i=0;i<simulCond->nconst;i++)
      {
	ia=constList[i].a;
	ib=constList[i].b;
	
	t2rmi=ts2/atom->m[ia];
	t2rmj=ts2/atom->m[ib];
	
	lambda=-(constList[i].rc2-rt2[i])/(2.*(t2rmi+t2rmj)*
	  ((dd[i].x*dt[i].x)+(dd[i].y*dt[i].y)+(dd[i].z*dt[i].z)));
	
	lambdai=lambda*t2rmi;
	xt[ia]+=dd[i].x*lambdai;
	yt[ia]+=dd[i].y*lambdai;
	zt[ia]+=dd[i].z*lambdai;
            
	lambdaj=-lambda*t2rmj;
	xt[ib]+=dd[i].x*lambdaj;
	yt[ib]+=dd[i].y*lambdaj;
	zt[ib]+=dd[i].z*lambdaj;
	
      }
      
      for(i=0;i<simulCond->nconst;i++)
      {
	ia=constList[i].a;
	ib=constList[i].b;
	
	nia=(double)atom->inconst[ia];
	nib=(double)atom->inconst[ib];
	
	atom->x[ia]+=xt[ia]/nia;
	atom->y[ia]+=yt[ia]/nia;
	atom->z[ia]+=zt[ia]/nia;
	
	atom->x[ib]+=xt[ib]/nib;
	atom->y[ib]+=yt[ib]/nib;
	atom->z[ib]+=zt[ib]/nib;
	
      }
      
    }
    
    icycle++;
    
  }while( (!converged) && (icycle<simulCond->maxcycle) );
  
  if(!converged)
    error(311);
  
  free(xt);
  free(yt);
  free(zt);
  free(rt2);
  free(dt);
  
}

void vv_shake_r(ATOM *atom,SIMULPARAMS *simulCond,CONSTRAINT *constList,DELTA *dd)
{
  int i,ia,ib,icycle,converged;
  double *xt,*yt,*zt,*rt2,maxdist,dist;
  double lambda,lambdai,lambdaj,trmi,trmj,nia,nib;
  DELTA *dt;
  
  xt=(double*)malloc(atom->natom*sizeof(*xt));
  yt=(double*)malloc(atom->natom*sizeof(*yt));
  zt=(double*)malloc(atom->natom*sizeof(*zt));
  
  rt2=(double*)malloc(simulCond->nconst*sizeof(*rt2));
  
  dt=(DELTA*)malloc(simulCond->nconst*sizeof(*dt));
  
  icycle=0;
  converged=0;

  do
  {
    maxdist=0.;
    
    for(i=0;i<simulCond->nconst;i++)
    {
      ia=constList[i].a;
      ib=constList[i].b;
      
      rt2[i]=distance2(ia,ib,atom,&(dt[i]),simulCond);
      
      dist=fabs(rt2[i]-constList[i].rc2)/sqrt(constList[i].rc2);
      maxdist=MAX(maxdist,dist);
    }
    
    maxdist=0.5*maxdist;
    
    if(maxdist<simulCond->tolshake)
      converged=1;
    
    if(!converged)
    {
      for(i=0;i<atom->natom;i++)
      {
	xt[i]=0.;
	yt[i]=0.;
	zt[i]=0.;
      }
      
      for(i=0;i<simulCond->nconst;i++)
      {
	ia=constList[i].a;
	ib=constList[i].b;
	
	trmi=simulCond->timeStep/atom->m[ia];
	trmj=simulCond->timeStep/atom->m[ib];
	
	lambda=-(constList[i].rc2-rt2[i])/(simulCond->timeStep*(trmi+trmj)*
	  ((dd[i].x*dt[i].x)+(dd[i].y*dt[i].y)+(dd[i].z*dt[i].z)));
	
	lambdai=0.5*lambda*trmi;
	xt[ia]+=dd[i].x*lambdai;
	yt[ia]+=dd[i].y*lambdai;
	zt[ia]+=dd[i].z*lambdai;
            
	lambdaj=-0.5*lambda*trmj;
	xt[ib]+=dd[i].x*lambdaj;
	yt[ib]+=dd[i].y*lambdaj;
	zt[ib]+=dd[i].z*lambdaj;
	
      }
      
      for(i=0;i<simulCond->nconst;i++)
      {
	ia=constList[i].a;
	ib=constList[i].b;
	
	nia=(double)atom->inconst[ia];
	nib=(double)atom->inconst[ib];
	
	atom->x[ia]+=simulCond->timeStep*xt[ia]/nia;
	atom->y[ia]+=simulCond->timeStep*yt[ia]/nia;
	atom->z[ia]+=simulCond->timeStep*zt[ia]/nia;
	
	atom->x[ib]+=simulCond->timeStep*xt[ib]/nib;
	atom->y[ib]+=simulCond->timeStep*yt[ib]/nib;
	atom->z[ib]+=simulCond->timeStep*zt[ib]/nib;
	
	atom->vx[ia]+=xt[ia]/nia;
	atom->vy[ia]+=yt[ia]/nia;
	atom->vz[ia]+=zt[ia]/nia;
	
	atom->vx[ib]+=xt[ib]/nib;
	atom->vy[ib]+=yt[ib]/nib;
	atom->vz[ib]+=zt[ib]/nib;
	
      }
      
    }
    
    icycle++;
    
  }while( (!converged) && (icycle<simulCond->maxcycle) );
  
  if(!converged)
    error(311);
  
  free(xt);
  free(yt);
  free(zt);
  free(rt2);
  free(dt);
  
}

void vv_shake_v(ATOM *atom,SIMULPARAMS *simulCond,CONSTRAINT *constList,DELTA *dd)
{
  int i,ia,ib,icycle,converged;
  double *xt,*yt,*zt,maxdist,tolvel;
  double lambda,lambdai,lambdaj,trmi,trmj,nia,nib;
  
  xt=(double*)malloc(atom->natom*sizeof(*xt));
  yt=(double*)malloc(atom->natom*sizeof(*yt));
  zt=(double*)malloc(atom->natom*sizeof(*zt));
  
  icycle=0;
  converged=0;
  
  tolvel=simulCond->tolshake/simulCond->timeStep;

  do
  {
    
    for(i=0;i<atom->natom;i++)
    {
      xt[i]=0.;
      yt[i]=0.;
      zt[i]=0.;
    }
    
    maxdist=0.;
    
    for(i=0;i<simulCond->nconst;i++)
    {
      ia=constList[i].a;
      ib=constList[i].b;
      
      trmi=0.5*simulCond->timeStep/atom->m[ia];
      trmj=0.5*simulCond->timeStep/atom->m[ib];
      
      lambda=(dd[i].x*(atom->vx[ib]-atom->vx[ia])+dd[i].y*(atom->vy[ib]-atom->vy[ia])+
	dd[i].z*(atom->vz[ib]-atom->vz[ia]))/((trmi+trmj)*
	(X2(dd[i].x)+X2(dd[i].y)+X2(dd[i].z)));
	
      maxdist=MAX(maxdist,fabs(lambda));
      
      lambdai=lambda*trmi;
      xt[ia]+=dd[i].x*lambdai;
      yt[ia]+=dd[i].y*lambdai;
      zt[ia]+=dd[i].z*lambdai;
	  
      lambdaj=-lambda*trmj;
      xt[ib]+=dd[i].x*lambdaj;
      yt[ib]+=dd[i].y*lambdaj;
      zt[ib]+=dd[i].z*lambdaj;
      
    }
    
    if(maxdist<tolvel)
      converged=1;
    
    if(!converged)
    {
      for(i=0;i<simulCond->nconst;i++)
      {
	ia=constList[i].a;
	ib=constList[i].b;
	
	nia=(double)atom->inconst[ia];
	nib=(double)atom->inconst[ib];
	
	atom->vx[ia]+=xt[ia]/nia;
	atom->vy[ia]+=yt[ia]/nia;
	atom->vz[ia]+=zt[ia]/nia;
	
	atom->vx[ib]+=xt[ib]/nib;
	atom->vy[ib]+=yt[ib]/nib;
	atom->vz[ib]+=zt[ib]/nib;
	
      }
    }
    
    icycle++;
    
  }while( (!converged) && (icycle<simulCond->maxcycle) );
  
  if(!converged)
    error(311);
  
  free(xt);
  free(yt);
  free(zt);
  
}