/**
 * \file utils.c
 * \brief Contains various and general utilitary functions.
 * \author Pierre-Andre Cazade and Florent Hedin
 * \version alpha-branch
 * \date 2012
 */

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>

#include "global.h"
#include "utils.h"
#include "rand.h"
#include "io.h"

/**
 * \param i Index of first atom.
 * \param j Index of second atom.
 * \param atom Array of structure ATOM (coordinates, forces, etc...).
 * \param delta Vector of length 3, representing the vector ij.
 * \param simulCond Pointer to structure SIMULPARAMS containing parameters of the current simulation.
 * \param box Pointer to structure PBC containing Periodic Boundaries Conditions parameters.
 *
 * \brief Evaluates the distance between two atoms i and j.
 * \remarks Periodic Boundaries Conditions are internally applied.
 *
 * \return The distance in Angstroems between two atoms.
 */
double distance(int i,int j, ATOM atom[],double delta[3],SIMULPARAMS *simulCond,PBC *box)
{
  double r,xt,yt,zt;
  
  delta[0]=atom[j].x-atom[i].x;
  delta[1]=atom[j].y-atom[i].y;
  delta[2]=atom[j].z-atom[i].z;
  
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
  
  r=sqrt(X2(delta[0])+X2(delta[1])+X2(delta[2]));
  
  return r;
  
}

/**
 * \param i Index of first atom.
 * \param j Index of second atom.
 * \param atom Array of structure ATOM (coordinates, forces, etc...).
 * \param d Pointer to structure DELTA representing the vector ij.
 * \param simulCond Pointer to structure SIMULPARAMS containing parameters of the current simulation.
 * \param box Pointer to structure PBC containing Periodic Boundaries Conditions parameters.
 *
 * \brief Evaluates the square of the distance between two atoms i and j.
 * \remarks Periodic Boundaries Conditions are internally applied.
 *
 * \return The squared distance in Angstroems^2 between two atoms.
 */
double distance2(int i,int j, ATOM atom[],DELTA *d,SIMULPARAMS *simulCond,PBC *box)
{
  double r2,xt,yt,zt;
  
  d->x=atom[j].x-atom[i].x;
  d->y=atom[j].y-atom[i].y;
  d->z=atom[j].z-atom[i].z;
  
  switch(box->type)
  {
    case CUBIC:
      d->x-=box->a1*nint(d->x*box->u1);
      d->y-=box->a1*nint(d->y*box->u1);
      d->z-=box->a1*nint(d->z*box->u1);
      break;
      
    case ORBIC:
      d->x-=box->a1*nint(d->x*box->u1);
      d->y-=box->b2*nint(d->y*box->v2);
      d->z-=box->c3*nint(d->z*box->w3);
    break;
      
    case TCLIN:
      xt=d->x*box->u1+d->y*box->u2+d->z*box->u3;
      yt=d->x*box->v1+d->y*box->v2+d->z*box->v3;
      zt=d->x*box->w1+d->y*box->w2+d->z*box->w3;
      
      xt-=nint(xt);
      yt-=nint(yt);
      zt-=nint(zt);
      
      d->x=xt*box->a1+yt*box->b1+zt*box->c1;
      d->y=xt*box->a2+yt*box->b2+zt*box->c2;
      d->z=xt*box->a3+yt*box->b3+zt*box->c3;
    break;
      
    default:
      break;
  }
  
  r2=( X2(d->x)+X2(d->y)+X2(d->z) );
  
  return r2;
  
}

/**
 * \param atom Array of structure ATOM (coordinates, forces, etc...).
 * \param simulCond Pointer to structure SIMULPARAMS containing parameters of the current simulation.
 * \param constList Pointer to structure CONSTRAINT containing parameters of constraints.
 * \param box Pointer to structure PBC containing Periodic Boundaries Conditions parameters.
 *
 * \brief Initialise velocity of atoms according to a random normal distribution.
 * \remarks If constraints are used, init_constvel is internally called.
 */
void init_vel(ATOM atom[],SIMULPARAMS *simulCond,CONSTRAINT *constList,PBC *box)
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
    init_constvel(atom,simulCond,constList,box);
  
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

/**
 * \param atom Array of structure ATOM (coordinates, forces, etc...).
 * \param simulCond Pointer to structure SIMULPARAMS containing parameters of the current simulation.
 * \param constList Pointer to structure CONSTRAINT containing parameters of constraints.
 * \param box Pointer to structure PBC containing Periodic Boundaries Conditions parameters.
 *
 * \brief Initialise velocity of atoms according to a random normal distribution.
 * \remarks Called by init_vel if constraints are used.
 */
void init_constvel(ATOM atom[],SIMULPARAMS *simulCond,CONSTRAINT *constList,PBC *box)
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
      
    rt=sqrt(distance2(ia,ib,atom,&(dt[i]),simulCond,box));
    
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

/**
 * \param atom Array of structure ATOM (coordinates, forces, etc...).
 * \param simulCond Pointer to structure SIMULPARAMS containing parameters of the current simulation.
 * \param box Pointer to structure PBC containing Periodic Boundaries Conditions parameters.
 *
 * \brief Applies Periodic Boundaries Conditions to atom[], according to the PBC type.
 */
void image_update(ATOM atom[],SIMULPARAMS *simulCond,PBC *box)
{

  int i;
  double xt,yt,zt;
  
  switch(box->type)
  {
    case CUBIC:
      #ifdef _OPENMP
      #pragma omp parallel for default(none) shared(simulCond,atom,box) private(i)
      #endif
      for(i=0;i<simulCond->natom;i++)
      {
	atom[i].x=atom[i].x-box->a1*nint(atom[i].x*box->u1);
	atom[i].y=atom[i].y-box->a1*nint(atom[i].y*box->u1);
	atom[i].z=atom[i].z-box->a1*nint(atom[i].z*box->u1);
      }
    break;
      
    case ORBIC:
      #ifdef _OPENMP
      #pragma omp parallel for default(none) shared(simulCond,atom,box) private(i)
      #endif
      for(i=0;i<simulCond->natom;i++)
      {
	atom[i].x=atom[i].x-box->a1*nint(atom[i].x*box->u1);
	atom[i].y=atom[i].y-box->b2*nint(atom[i].y*box->v2);
	atom[i].z=atom[i].z-box->c3*nint(atom[i].z*box->w3);
      }
    break;
      
    case TCLIN:
      #ifdef _OPENMP
      #pragma omp parallel for default(none) shared(simulCond,atom,box) private(i,xt,yt,zt)
      #endif
      for(i=0;i<simulCond->natom;i++)
      {
	xt=atom[i].x*box->u1+atom[i].y*box->u2+atom[i].z*box->u3;
	yt=atom[i].x*box->v1+atom[i].y*box->v2+atom[i].z*box->v3;
	zt=atom[i].x*box->w1+atom[i].y*box->w2+atom[i].z*box->w3;
	
	xt-=nint(xt);
	yt-=nint(yt);
	zt-=nint(zt);
	
	atom[i].x=xt*box->a1+yt*box->b1+zt*box->c1;
	atom[i].y=xt*box->a2+yt*box->b2+zt*box->c2;
	atom[i].z=xt*box->a3+yt*box->b3+zt*box->c3;
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
void image_array(int size_array,DELTA d[],SIMULPARAMS *simulCond,PBC *box)
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
	d[i].x-=box->a1*nint(d[i].x*box->u1);
	d[i].y-=box->a1*nint(d[i].y*box->u1);
	d[i].z-=box->a1*nint(d[i].z*box->u1);
      }
    break;
      
    case ORBIC:
      #ifdef _OPENMP
      #pragma omp parallel for default(none) shared(size_array,box,d) private(i)
      #endif
      for(i=0;i<size_array;i++)
      {
	d[i].x-=box->a1*nint(d[i].x*box->u1);
	d[i].y-=box->b2*nint(d[i].y*box->v2);
	d[i].z-=box->c3*nint(d[i].z*box->w3);
      }
    break;
      
    case TCLIN:
      #ifdef _OPENMP
      #pragma omp parallel for default(none) shared(size_array,box,d) private(i,xt,yt,zt)
      #endif
      for(i=0;i<size_array;i++)
      {
	xt=d[i].x*box->u1+d[i].y*box->u2+d[i].z*box->u3;
	yt=d[i].x*box->v1+d[i].y*box->v2+d[i].z*box->v3;
	zt=d[i].x*box->w1+d[i].y*box->w2+d[i].z*box->w3;
	
	xt-=nint(xt);
	yt-=nint(yt);
	zt-=nint(zt);
	
	d[i].x=xt*box->a1+yt*box->b1+zt*box->c1;
	d[i].y=xt*box->a2+yt*box->b2+zt*box->c2;
	d[i].z=xt*box->a3+yt*box->b3+zt*box->c3;
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
void init_box(PBC *box)
{
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
  box->vol0=box->vol; //for NPT and NST purpose
  
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
void scale_box(PBC *box,double scale,double cell0[9])
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
  
  printf("%lf %lf %lf",box->a1,box->b2,box->c3);
  
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
void vv_scale_box(PBC *box,double scale)
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
  
  printf("%lf %lf %lf",box->a1,box->b2,box->c3);
  
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
 * \param atom Array of structure ATOM (coordinates, forces, etc...).
 * \param simulCond Pointer to structure SIMULPARAMS containing parameters of the current simulation.
 *
 * \brief Estimates the kinetic energy of a system.
 *
 * \return On return, the kinetic energy : \f$ 1/2*m*v^2 \f$
 */
double kinetic(ATOM atom[],SIMULPARAMS *simulCond)
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

/**
 * \param atom Array of structure ATOM (coordinates, forces, etc...).
 * \param simulCond Pointer to structure SIMULPARAMS containing parameters of the current simulation.
 * \param stress Array representing (parts of) a Stress Tensor.
 *
 * \brief Adds to the stress tensor the kinetic energy.
 */
void stress_kinetic(ATOM atom[],SIMULPARAMS *simulCond,double stress[6])
{
  int i;
  
  for(i=0;i<6;i++)
    stress[i]=0.;
  
  for(i=0;i<simulCond->natom;i++)
  {
    stress[0] += atom[i].m*X2(atom[i].vx);
    stress[1] += atom[i].m*atom[i].vx*atom[i].vy;
    stress[2] += atom[i].m*atom[i].vx*atom[i].vz;
    stress[3] += atom[i].m*X2(atom[i].vy);
    stress[4] += atom[i].m*atom[i].vy*atom[i].vz;
    stress[5] += atom[i].m*X2(atom[i].vz);
  }
  
}

/**
 * \param atom Array of structure ATOM (coordinates, forces, etc...).
 * \param simulCond Pointer to structure SIMULPARAMS containing parameters of the current simulation.
 * \param box Pointer to structure PBC containing Periodic Boundaries Conditions parameters.
 *
 * \brief Get initial kinetic energy from the temperature.
 */
void get_kinfromtemp(ATOM atom[],SIMULPARAMS *simulCond,PBC *box)
{
  double degf;
  
  get_degfree(simulCond,box);
  
  //   Energy in internal units 10 J/mol. rboltzui=R/10.
  
  degf=(double)simulCond->degfree;
  simulCond->kintemp0=0.5*simulCond->temp*degf*rboltzui;
}

/**
 * \param simulCond Pointer to structure SIMULPARAMS containing parameters of the current simulation.
 * \param box Pointer to structure PBC containing Periodic Boundaries Conditions parameters.
 *
 * \brief Obtains the number of degrees of freedom of the system.
 */
void get_degfree(SIMULPARAMS *simulCond,PBC *box)
{
  
//   Atoms degrees of freedom - 3 degrees of freedom of the CoM.
  simulCond->degfree=3*simulCond->natom-3;
  
//   -3 degrees of freedpm for the box rotation when no PBC.
  if(box->type==NOBOX)
    simulCond->degfree-=3;
  
  if(simulCond->nconst>0)
    simulCond->degfree-=simulCond->nconst;
    
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
int nint(double x)
{
  return ( (x)>=0.?(int)((x)+0.5):(int)((x)-0.5) );
}
