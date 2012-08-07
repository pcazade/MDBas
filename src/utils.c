#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include "global.h"
#include "utils.h"
#include "rand.h"
#include "io.h"

double distance(int i,int j, ATOM *atom,double *delta,SIMULPARAMS *simulCond)
{
  double r;
  
  delta[0]=atom[j].x-atom[i].x;
  delta[1]=atom[j].y-atom[i].y;
  delta[2]=atom[j].z-atom[i].z;
  
  if(simulCond->periodicType==1)
  {
    delta[0]=delta[0]-simulCond->periodicBox[0][0]*nint(delta[0]/simulCond->periodicBox[0][0]);
    delta[1]=delta[1]-simulCond->periodicBox[0][0]*nint(delta[1]/simulCond->periodicBox[0][0]);
    delta[2]=delta[2]-simulCond->periodicBox[0][0]*nint(delta[2]/simulCond->periodicBox[0][0]);
  }
  else if(simulCond->periodicType==2)
  {
    
    delta[0]=delta[0]-simulCond->periodicBox[0][0]*nint(delta[0]/simulCond->periodicBox[0][0]);
    delta[1]=delta[1]-simulCond->periodicBox[1][1]*nint(delta[1]/simulCond->periodicBox[1][1]);
    delta[2]=delta[2]-simulCond->periodicBox[2][2]*nint(delta[2]/simulCond->periodicBox[2][2]);
  }
  
  r=sqrt(X2(delta[0])+X2(delta[1])+X2(delta[2]));
  
  return r;
  
}

double distance2(int i,int j, ATOM *atom,DELTA *d,SIMULPARAMS *simulCond)
{
  double r2;
  
  d->x=atom[j].x-atom[i].x;
  d->y=atom[j].y-atom[i].y;
  d->z=atom[j].z-atom[i].z;
  
  if(simulCond->periodicType==1)
  {
    d->x-=simulCond->periodicBox[0][0]*nint(d->x/simulCond->periodicBox[0][0]);
    d->y-=simulCond->periodicBox[0][0]*nint(d->y/simulCond->periodicBox[0][0]);
    d->z-=simulCond->periodicBox[0][0]*nint(d->z/simulCond->periodicBox[0][0]);
  }
  else if(simulCond->periodicType==2)
  {
    d->x-=simulCond->periodicBox[0][0]*nint(d->x/simulCond->periodicBox[0][0]);
    d->y-=simulCond->periodicBox[1][1]*nint(d->y/simulCond->periodicBox[1][1]);
    d->z-=simulCond->periodicBox[2][2]*nint(d->z/simulCond->periodicBox[2][2]);
  }
  
  r2=( X2(d->x)+X2(d->y)+X2(d->z) );
  
  return r2;
  
}

void init_vel(ATOM *atom,SIMULPARAMS *simulCond,CONSTRAINT *constList)
{
  int i,natoms;
  double cmvx=0.,cmvy=0.,cmvz=0.,cmm=0.,initKin=0.;
  double factor;  
  
  natoms=2*(simulCond->natom/2);
  
//   Obtain gaussian distribution of velocities via Box-Muller method.
  
  for(i=0;i<natoms;i+=2)
  {
    get_BoxMuller(&(atom[i].vx),&(atom[i+1].vx));
    get_BoxMuller(&(atom[i].vy),&(atom[i+1].vy));
    get_BoxMuller(&(atom[i].vz),&(atom[i+1].vz));
  }
  
  if(natoms!=simulCond->natom)
  {
    double v;
    get_BoxMuller(&(atom[natoms].vx),&v);
    get_BoxMuller(&(atom[natoms].vy),&v);
    get_BoxMuller(&(atom[natoms].vz),&v);
  }
  
//   Random numers provided by the generator have to be weighted
//   by the square root of the atoms mass to obtain velocities:
//   P(vx_i) ~ exp(-m_i*vx_i**2/(2*kb*T))
  
  for(i=0;i<simulCond->natom;i++)
  {
    atom[i].vx/=sqrt(atom[i].m);
    atom[i].vy/=sqrt(atom[i].m);
    atom[i].vz/=sqrt(atom[i].m);
  }
  
  if(simulCond->nconst>0)
    init_constvel(atom,simulCond,constList);
  
//   Set the system total momentum to zero to acheive momentum conservation.
  
  for(i=0;i<simulCond->natom;i++)
  {
    cmvx+=atom[i].vx*atom[i].m;
    cmvy+=atom[i].vy*atom[i].m;
    cmvz+=atom[i].vz*atom[i].m;
    
    cmm+=atom[i].m;
  }
  
  cmvx/=cmm;
  cmvy/=cmm;
  cmvz/=cmm;
  
  for(i=0;i<simulCond->natom;i++)
  {
    atom[i].vx-=cmvx;
    atom[i].vy-=cmvy;
    atom[i].vz-=cmvz;
  }
  
  initKin=kinetic(atom,simulCond);
  
  factor=sqrt(simulCond->kintemp0/initKin);
  
  for(i=0;i<simulCond->natom;i++)
  {
    atom[i].vx*=factor;
    atom[i].vy*=factor;
    atom[i].vz*=factor;
  }
  
}

void init_constvel(ATOM *atom,SIMULPARAMS *simulCond,CONSTRAINT *constList)
{
  int i,ia,ib,icycle,converged;
  double *vxu=NULL,*vyu=NULL,*vzu=NULL;
  double rt,maxdv,dv,w1,w2,nia,nib;
  DELTA *dt=NULL;
  
  vxu=(double*)malloc(simulCond->natom*sizeof(*vxu));
  vyu=(double*)malloc(simulCond->natom*sizeof(*vyu));
  vzu=(double*)malloc(simulCond->natom*sizeof(*vzu));
  
  dt=(DELTA*)malloc(simulCond->nconst*sizeof(*dt));
  
  for(i=0;i<simulCond->nconst;i++)
  {
    ia=constList[i].a;
    ib=constList[i].b;
      
    rt=sqrt(distance2(ia,ib,atom,&(dt[i]),simulCond));
    
    dt[i].x/=rt;
    dt[i].y/=rt;
    dt[i].z/=rt;
    
  }
  
  icycle=0;
  converged=0;
  
  do
  {
    for(i=0;i<simulCond->natom;i++)
    {
      vxu[i]=0.;
      vyu[i]=0.;
      vzu[i]=0.;
    }
    
    maxdv=0.;
    
    for(i=0;i<simulCond->nconst;i++)
    {
      
      ia=constList[i].a;
      ib=constList[i].b;
      
      dv=dt[i].x*(atom[ib].vx-atom[ia].vx)+dt[i].y*(atom[ib].vy-atom[ia].vy)+
	 dt[i].z*(atom[ib].vz-atom[ia].vz);
      
      maxdv=MAX(maxdv,fabs(dv));
      
      w1=atom[ib].m*dv/(atom[ia].m+atom[ib].m);
      w2=atom[ia].m*dv/(atom[ia].m+atom[ib].m);
      
      vxu[ia]+=w1*dt[i].x;
      vyu[ia]+=w1*dt[i].y;
      vzu[ia]+=w1*dt[i].z;
      
      vxu[ib]-=w2*dt[i].x;
      vyu[ib]-=w2*dt[i].y;
      vzu[ib]-=w2*dt[i].z;
      
    }
    
    if(maxdv<simulCond->tolshake)
      converged=1;
    
    if(!converged)
    {
      for(i=0;i<simulCond->nconst;i++)
      {
	ia=constList[i].a;
	ib=constList[i].b;
	
	nia=(double)atom[ia].inconst;
	nib=(double)atom[ib].inconst;
	
	atom[ia].vx+=vxu[ia]/nia;
	atom[ia].vy+=vyu[ia]/nia;
	atom[ia].vz+=vzu[ia]/nia;
	
	atom[ib].vx+=vxu[ib]/nib;
	atom[ib].vy+=vyu[ib]/nib;
	atom[ib].vz+=vzu[ib]/nib;
      }
    }
    
    icycle++;
    
  }while( (!converged) && (icycle<simulCond->maxcycle) );
  
  if(!converged)
    error(310);
  
  free(vxu);
  free(vyu);
  free(vzu);
  free(dt);
}

void image_update(ATOM *atom,SIMULPARAMS *simulCond)
{

  int i;
  
  if(simulCond->periodicType==1)
  {
    
    for(i=0;i<simulCond->natom;i++)
    {
      atom[i].x=atom[i].x-simulCond->periodicBox[0][0]*nint(atom[i].x/simulCond->periodicBox[0][0]);
      atom[i].y=atom[i].y-simulCond->periodicBox[0][0]*nint(atom[i].y/simulCond->periodicBox[0][0]);
      atom[i].z=atom[i].z-simulCond->periodicBox[0][0]*nint(atom[i].z/simulCond->periodicBox[0][0]);
    }
    
  }
  else if(simulCond->periodicType==2)
  {
    
    for(i=0;i<simulCond->natom;i++)
    {
      atom[i].x=atom[i].x-simulCond->periodicBox[0][0]*nint(atom[i].x/simulCond->periodicBox[0][0]);
      atom[i].y=atom[i].y-simulCond->periodicBox[1][1]*nint(atom[i].y/simulCond->periodicBox[1][1]);
      atom[i].z=atom[i].z-simulCond->periodicBox[2][2]*nint(atom[i].z/simulCond->periodicBox[2][2]);
    }
    
  }
  
}

void image_array(int size_array,DELTA *d,SIMULPARAMS *simulCond)
{

  int i;
  
  if(simulCond->periodicType==1)
  {
    
    for(i=0;i<size_array;i++)
    {
      d[i].x-=simulCond->periodicBox[0][0]*nint(d[i].x/simulCond->periodicBox[0][0]);
      d[i].y-=simulCond->periodicBox[0][0]*nint(d[i].y/simulCond->periodicBox[0][0]);
      d[i].z-=simulCond->periodicBox[0][0]*nint(d[i].z/simulCond->periodicBox[0][0]);
    }
    
  }
  else if(simulCond->periodicType==2)
  {
    
    for(i=0;i<size_array;i++)
    {
      d[i].x-=simulCond->periodicBox[0][0]*nint(d[i].x/simulCond->periodicBox[0][0]);
      d[i].y-=simulCond->periodicBox[1][1]*nint(d[i].y/simulCond->periodicBox[1][1]);
      d[i].z-=simulCond->periodicBox[2][2]*nint(d[i].z/simulCond->periodicBox[2][2]);
    }
    
  }
  
}

double kinetic(ATOM *atom,SIMULPARAMS *simulCond)
{
  int i;
  double ekin;
  
  ekin=0.;
  for(i=0;i<simulCond->natom;i++)
  {
    ekin+=atom[i].m*(X2(atom[i].vx)+X2(atom[i].vy)+X2(atom[i].vz));
  }
  
  return ( ekin*0.5 );
}

void get_kinfromtemp(ATOM *atom,SIMULPARAMS *simulCond)
{
  double degf;
  
  get_degfree(atom,simulCond);
  
  //   Energy in internal units 10 J/mol. rboltzui=R/10.
  
  degf=(double)simulCond->degfree;
  simulCond->kintemp0=0.5*simulCond->temp*degf*rboltzui;
}

void get_degfree(ATOM *atom,SIMULPARAMS *simulCond)
{
  
//   Atoms degrees of freedom - 3 degrees of freedom of the CoM.
  simulCond->degfree=3*simulCond->natom-3;
  
//   -3 degrees of freedpm for the box rotation when no PBC.
  if(simulCond->periodicType==0)
    simulCond->degfree-=3;
  
  if(simulCond->nconst>0)
    simulCond->degfree-=simulCond->nconst;
    
}

void nocase(char *str)
{
  int i,n;
  n=strlen(str);
  for(i=0;i<n;i++)
  {
    str[i]=tolower(str[i]);
  }
}

int nint(double x)
{
  return ( (x)>=0.?(int)((x)+0.5):(int)((x)-0.5) );
}
