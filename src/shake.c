#include <stdlib.h>
#include <math.h>
#include "global.h"
#include "utils.h"
#include "io.h"

void lf_shake(ATOM *atom,SIMULPARAMS *simulCond,CONSTRAINT *constList,DELTA *dd,PBC *box)
{
  int i,ia,ib,icycle,converged;
  double *xt,*yt,*zt,*rt2,ts2,maxdist,dist;
  double lambda,lambdai,lambdaj,t2rmi,t2rmj,nia,nib;
  double virshake=0.,stress[6]={0.};
  DELTA *dt;
  
  xt=(double*)malloc(simulCond->natom*sizeof(*xt));
  yt=(double*)malloc(simulCond->natom*sizeof(*yt));
  zt=(double*)malloc(simulCond->natom*sizeof(*zt));
  
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
      
      rt2[i]=distance2(ia,ib,atom,&(dt[i]),simulCond,box);
      
      dist=fabs(rt2[i]-constList[i].rc2)/sqrt(constList[i].rc2);
      maxdist=MAX(maxdist,dist);
    }
    
    maxdist=0.5*maxdist;
    
    if(maxdist<simulCond->tolshake)
      converged=1;
    
    if(!converged)
    {
      for(i=0;i<simulCond->natom;i++)
      {
	xt[i]=0.;
	yt[i]=0.;
	zt[i]=0.;
      }
      
      for(i=0;i<simulCond->nconst;i++)
      {
	ia=constList[i].a;
	ib=constList[i].b;
	
	t2rmi=ts2/atom[ia].m;
	t2rmj=ts2/atom[ib].m;
	
	lambda=-(constList[i].rc2-rt2[i])/(2.*(t2rmi+t2rmj)*
	  ((dd[i].x*dt[i].x)+(dd[i].y*dt[i].y)+(dd[i].z*dt[i].z)));
	
	virshake+=lambda*(X2(dd[i].x)+X2(dd[i].y)+X2(dd[i].z));
	
	stress[0]-=lambda*X2(dd[i].x);
	stress[1]-=lambda*dd[i].x*dd[i].y;
	stress[2]-=lambda*dd[i].x*dd[i].z;
	stress[3]-=lambda*X2(dd[i].y);
	stress[4]-=lambda*dd[i].y*dd[i].z;
	stress[5]-=lambda*X2(dd[i].z);
	
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
	
	nia=(double)atom[ia].inconst;
	nib=(double)atom[ib].inconst;
	
	atom[ia].x+=xt[ia]/nia;
	atom[ia].y+=yt[ia]/nia;
	atom[ia].z+=zt[ia]/nia;
	
	atom[ib].x+=xt[ib]/nib;
	atom[ib].y+=yt[ib]/nib;
	atom[ib].z+=zt[ib]/nib;
	
      }
      
    }
    
    icycle++;
    
  }while( (!converged) && (icycle<simulCond->maxcycle) );
  
  if(!converged)
    error(311);
  /*
  ener->virshake+=virshake;
  
  box->stress1+=stress[0];
  box->stress2+=stress[1];
  box->stress3+=stress[2];
  box->stress4+=stress[1];
  box->stress5+=stress[3];
  box->stress6+=stress[4];
  box->stress7+=stress[2];
  box->stress8+=stress[4];
  box->stress9+=stress[5];
  */
  
  free(xt);
  free(yt);
  free(zt);
  free(rt2);
  free(dt);
  
}

void vv_shake_r(ATOM *atom,SIMULPARAMS *simulCond,CONSTRAINT *constList,DELTA *dd,PBC *box)
{
  int i,ia,ib,icycle,converged;
  double *xt,*yt,*zt,*rt2,maxdist,dist;
  double lambda,lambdai,lambdaj,trmi,trmj,nia,nib;
  DELTA *dt;
  
  xt=(double*)malloc(simulCond->natom*sizeof(*xt));
  yt=(double*)malloc(simulCond->natom*sizeof(*yt));
  zt=(double*)malloc(simulCond->natom*sizeof(*zt));
  
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
      
      rt2[i]=distance2(ia,ib,atom,&(dt[i]),simulCond,box);
      
      dist=fabs(rt2[i]-constList[i].rc2)/sqrt(constList[i].rc2);
      maxdist=MAX(maxdist,dist);
    }
    
    maxdist=0.5*maxdist;
    
    if(maxdist<simulCond->tolshake)
      converged=1;
    
    if(!converged)
    {
      for(i=0;i<simulCond->natom;i++)
      {
	xt[i]=0.;
	yt[i]=0.;
	zt[i]=0.;
      }
      
      for(i=0;i<simulCond->nconst;i++)
      {
	ia=constList[i].a;
	ib=constList[i].b;
	
	trmi=simulCond->timeStep/atom[ia].m;
	trmj=simulCond->timeStep/atom[ib].m;
	
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
	
	nia=(double)atom[ia].inconst;
	nib=(double)atom[ib].inconst;
	
	atom[ia].x+=simulCond->timeStep*xt[ia]/nia;
	atom[ia].y+=simulCond->timeStep*yt[ia]/nia;
	atom[ia].z+=simulCond->timeStep*zt[ia]/nia;
	
	atom[ib].x+=simulCond->timeStep*xt[ib]/nib;
	atom[ib].y+=simulCond->timeStep*yt[ib]/nib;
	atom[ib].z+=simulCond->timeStep*zt[ib]/nib;
	
	atom[ia].vx+=xt[ia]/nia;
	atom[ia].vy+=yt[ia]/nia;
	atom[ia].vz+=zt[ia]/nia;
	
	atom[ib].vx+=xt[ib]/nib;
	atom[ib].vy+=yt[ib]/nib;
	atom[ib].vz+=zt[ib]/nib;
	
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
  
  xt=(double*)malloc(simulCond->natom*sizeof(*xt));
  yt=(double*)malloc(simulCond->natom*sizeof(*yt));
  zt=(double*)malloc(simulCond->natom*sizeof(*zt));
  
  icycle=0;
  converged=0;
  
  tolvel=simulCond->tolshake/simulCond->timeStep;

  do
  {
    
    for(i=0;i<simulCond->natom;i++)
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
      
      trmi=0.5*simulCond->timeStep/atom[ia].m;
      trmj=0.5*simulCond->timeStep/atom[ib].m;
      
      lambda=(dd[i].x*(atom[ib].vx-atom[ia].vx)+dd[i].y*(atom[ib].vy-atom[ia].vy)+
	dd[i].z*(atom[ib].vz-atom[ia].vz))/((trmi+trmj)*
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
	
	nia=(double)atom[ia].inconst;
	nib=(double)atom[ib].inconst;
	
	atom[ia].vx+=xt[ia]/nia;
	atom[ia].vy+=yt[ia]/nia;
	atom[ia].vz+=zt[ia]/nia;
	
	atom[ib].vx+=xt[ib]/nib;
	atom[ib].vy+=yt[ib]/nib;
	atom[ib].vz+=zt[ib]/nib;
	
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