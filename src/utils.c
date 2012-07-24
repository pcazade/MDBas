#include <math.h>
#include <string.h>
#include <ctype.h>
#include "global.h"
#include "utils.h"
#include "rand.h"

double distance(int i,int j, ATOM *atom,double *delta,SIMULPARAMS *simulCond)
{
  double r;
  
  delta[0]=atom->x[j]-atom->x[i];
  delta[1]=atom->y[j]-atom->y[i];
  delta[2]=atom->z[j]-atom->z[i];
  
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
  
  d->x=atom->x[j]-atom->x[i];
  d->y=atom->y[j]-atom->y[i];
  d->z=atom->z[j]-atom->z[i];
  
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

void init_vel(ATOM *atom,SIMULPARAMS *simulCond)
{
  int i,natoms;
  double cmvx=0.,cmvy=0.,cmvz=0.,cmm=0.,initKin=0.;
  double factor;  
  
  natoms=2*(atom->natom/2);
  
//   Obtain gaussian distribution of velocities via Box-Muller method.
  
  for(i=0;i<natoms;i+=2)
  {
    get_BoxMuller(&(atom->vx[i]),&(atom->vx[i+1]));
    get_BoxMuller(&(atom->vy[i]),&(atom->vy[i+1]));
    get_BoxMuller(&(atom->vz[i]),&(atom->vz[i+1]));
  }
  
  if(natoms!=atom->natom)
  {
    double v;
    get_BoxMuller(&(atom->vx[natoms]),&v);
    get_BoxMuller(&(atom->vy[natoms]),&v);
    get_BoxMuller(&(atom->vz[natoms]),&v);
  }
  
//   Random numers provided by the generator have to be weighted
//   by the square root of the atoms mass to obtain velocities:
//   P(vx_i) ~ exp(-m_i*vx_i**2/(2*kb*T))
  
  for(i=0;i<atom->natom;i++)
  {
    atom->vx[i]/=sqrt(atom->m[i]);
    atom->vy[i]/=sqrt(atom->m[i]);
    atom->vz[i]/=sqrt(atom->m[i]);
  }
  
//   Set the system total momentum to zero to acheive momentum conservation.
  
  for(i=0;i<atom->natom;i++)
  {
    cmvx+=atom->vx[i]*atom->m[i];
    cmvy+=atom->vy[i]*atom->m[i];
    cmvz+=atom->vz[i]*atom->m[i];
    
    cmm+=atom->m[i];
  }
  
  cmvx/=cmm;
  cmvy/=cmm;
  cmvz/=cmm;
  
  for(i=0;i<atom->natom;i++)
  {
    atom->vx[i]-=cmvx;
    atom->vy[i]-=cmvy;
    atom->vz[i]-=cmvz;
  }
  
  initKin=kinetic(atom);
  
  factor=sqrt(simulCond->kintemp0/initKin);
  
  for(i=0;i<atom->natom;i++)
  {
    atom->vx[i]*=factor;
    atom->vy[i]*=factor;
    atom->vz[i]*=factor;
  }
  
}

void image_update(ATOM *atom,SIMULPARAMS *simulCond)
{

  int i;
  
  if(simulCond->periodicType==1)
  {
    
    for(i=0;i<atom->natom;i++)
    {
      atom->x[i]=atom->x[i]-simulCond->periodicBox[0][0]*nint(atom->x[i]/simulCond->periodicBox[0][0]);
      atom->y[i]=atom->y[i]-simulCond->periodicBox[0][0]*nint(atom->y[i]/simulCond->periodicBox[0][0]);
      atom->z[i]=atom->z[i]-simulCond->periodicBox[0][0]*nint(atom->z[i]/simulCond->periodicBox[0][0]);
    }
    
  }
  else if(simulCond->periodicType==2)
  {
    
    for(i=0;i<atom->natom;i++)
    {
      atom->x[i]=atom->x[i]-simulCond->periodicBox[0][0]*nint(atom->x[i]/simulCond->periodicBox[0][0]);
      atom->y[i]=atom->y[i]-simulCond->periodicBox[1][1]*nint(atom->y[i]/simulCond->periodicBox[1][1]);
      atom->z[i]=atom->z[i]-simulCond->periodicBox[2][2]*nint(atom->z[i]/simulCond->periodicBox[2][2]);
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

double kinetic(ATOM *atom)
{
  int i;
  double ekin;
  
  ekin=0.;
  for(i=0;i<atom->natom;i++)
  {
    ekin+=atom->m[i]*(X2(atom->vx[i])+X2(atom->vy[i])+X2(atom->vz[i]));
  }
  
  return ( ekin*0.5 );
}

void get_kinfromtemp(ATOM *atom,SIMULPARAMS *simulCond)
{
  get_degfree(atom,simulCond);
  
  //   Energy in internal units 10 J/mol. rboltzui=R/10.
  
  simulCond->kintemp0=0.5*simulCond->temp*simulCond->degfree*rboltzui;
}

void get_degfree(ATOM *atom,SIMULPARAMS *simulCond)
{
  
//   Atoms degrees of freedom - 3 degrees of freedom of the CoM.
  simulCond->degfree=3*atom->natom-3;
  
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