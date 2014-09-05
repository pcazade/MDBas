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

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#include "global.h"
#include "errors.h"
#include "memory.h"

#include "parallel.h"

#define STRBUFSIZ 8192

extern FILE *outFile;

void init_para(int *argc, char ***argv)
{
  int err;
  
  err=MPI_Init(argc,argv);
  
  MPI_Errhandler_set(MPI_COMM_WORLD, MPI_ERRORS_RETURN);
  
  if(err!=MPI_SUCCESS)
    mpi_error(err,__FILE__,__LINE__);
  
}

int my_proc()
{
  
  int idProc,err;
  
  err=MPI_Comm_rank(MPI_COMM_WORLD,&idProc);
  
  if(err!=MPI_SUCCESS)
    mpi_error(err,__FILE__,__LINE__);
  
  return (idProc);
  
}

int num_proc()
{
  int numProc,err;
  
  err=MPI_Comm_size(MPI_COMM_WORLD, &numProc);
  
  if(err!=MPI_SUCCESS)
    mpi_error(err,__FILE__,__LINE__);
  
  return(numProc);
  
}

void barrier_para()
{
  int err;
  
  err=MPI_Barrier(MPI_COMM_WORLD);
  
  if(err!=MPI_SUCCESS)
    mpi_error(err,__FILE__,__LINE__);
  
}

void sum_double_para(real *buf1,real *buf2,int size)
{
  int i,err;
  
  err=MPI_Allreduce(buf1,buf2,size,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  
  if(err!=MPI_SUCCESS)
    mpi_error(err,__FILE__,__LINE__);
  
  for(i=0;i<size;i++)
  {
    buf1[i]=buf2[i];
  }
  
}

void sum_int_para(int *buf1,int *buf2,int size)
{
  int i,err;
  
  err=MPI_Allreduce(buf1,buf2,size,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  
  if(err!=MPI_SUCCESS)
    mpi_error(err,__FILE__,__LINE__);
  
  for(i=0;i<size;i++)
  {
    buf1[i]=buf2[i];
  }
  
}

void update_double_para(PARAM *param,PARALLEL *parallel,real *buf1,real *buf2)
{
  int i,err;
  
  for(i=0;i<parallel->fAtProc;i++)
  {
    buf2[i]=0.;
  }
  
  for(i=parallel->fAtProc;i<parallel->lAtProc;i++)
  {
    buf2[i]=buf1[i];
  }
  
  for(i=parallel->lAtProc;i<param->nAtom;i++)
  {
    buf2[i]=0.;
  }
  
  err=MPI_Allreduce(buf2,buf1,param->nAtom,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  
  if(err!=MPI_SUCCESS)
    mpi_error(err,__FILE__,__LINE__);
  
}

void test_para(int *buf1)
{
  int test,err;
  
  err=MPI_Allreduce(buf1,&test,1,MPI_INT,MPI_LAND,MPI_COMM_WORLD);
  
  if(err!=MPI_SUCCESS)
    mpi_error(err,__FILE__,__LINE__);
  
    *buf1=test;
  
}

void parallel_allocate_buffers(PARAM *param,PARALLEL *parallel,real **dBuffer,
			       int **iBuffer)
{
  parallel->iBufferSize=128;
  parallel->dBufferSize=128;
  
  free(*iBuffer);
  free(*dBuffer);
  
  parallel->iBufferSize=MAX(param->nAtom,parallel->iBufferSize);
  parallel->dBufferSize=MAX(param->nAtom,parallel->dBufferSize);
  
  if(param->nConst>0)
  {
    parallel->iBufferSize=MAX(2*param->nConst,parallel->iBufferSize);
    parallel->dBufferSize=MAX(param->nConst,parallel->dBufferSize);
  }
  
  if(param->nBond>0)
  {
    parallel->iBufferSize=MAX(3*param->nBond,parallel->iBufferSize);
    parallel->dBufferSize=MAX(3*param->nBond,parallel->dBufferSize);
  }
  
  if(param->nAngle>0)
  {
    parallel->iBufferSize=MAX(3*param->nAngle,parallel->iBufferSize);
    parallel->dBufferSize=MAX(2*param->nAngle,parallel->dBufferSize);
  }
  
  if(param->nUb>0)
  {
    parallel->iBufferSize=MAX(3*param->nUb,parallel->iBufferSize);
    parallel->dBufferSize=MAX(3*param->nUb,parallel->dBufferSize);
  }
  
  if(param->nDihedral>0)
  {
    parallel->iBufferSize=MAX(6*param->nDihedral,parallel->iBufferSize);
    parallel->dBufferSize=MAX(3*param->nDihedral,parallel->dBufferSize);
  }
  
  if(param->nImproper>0)
  {
    parallel->iBufferSize=MAX(6*param->nImproper,parallel->iBufferSize);
    parallel->dBufferSize=MAX(3*param->nImproper,parallel->dBufferSize);
  }
  
  (*iBuffer)=(int*)my_malloc(parallel->iBufferSize*sizeof(**iBuffer));
  (*dBuffer)=(real*)my_malloc(parallel->dBufferSize*sizeof(**dBuffer));
  
}

void parallel_reallocate_buffers(PARALLEL *parallel,EWALD *ewald,real **dBuffer)
{
  
  if(ewald->mmax>parallel->dBufferSize)
  { 
    parallel->dBufferSize=ewald->mmax;
    
    (*dBuffer)=(real*)my_realloc(*dBuffer,parallel->dBufferSize*sizeof(**dBuffer));
  }
  
}

void bcast_int_para(int *buf1,int size,int iNode)
{
  int err;
  
  err=MPI_Bcast(buf1,size,MPI_INT,iNode,MPI_COMM_WORLD);
  
  if(err!=MPI_SUCCESS)
    mpi_error(err,__FILE__,__LINE__);
  
}

void bcast_double_para(real *buf1,int size,int iNode)
{
  int err;
  
  err=MPI_Bcast(buf1,size,MPI_DOUBLE,iNode,MPI_COMM_WORLD);
  
  if(err!=MPI_SUCCESS)
    mpi_error(err,__FILE__,__LINE__);
  
}

void bcast_param_para(PARAM *param,PARALLEL *parallel,real *dBuffer,int *iBuffer)
{
    
  if(parallel->idProc==0)
  {
    iBuffer[0]=param->step;
    iBuffer[1]=param->nSteps;
    iBuffer[2]=param->nAtom;
    iBuffer[3]=param->nDegFree;
    iBuffer[4]=param->nFrozen;
    iBuffer[5]=param->nBond;
    iBuffer[6]=param->nAngle;
    iBuffer[7]=param->nUb;
    iBuffer[8]=param->nDihedral;
    iBuffer[9]=param->nImproper;
    iBuffer[10]=param->nConst;
    iBuffer[11]=param->maxCycle;
    
    dBuffer[0]=param->timeStep;
    dBuffer[1]=param->rTimeStep;
    dBuffer[2]=param->chargeConst;
    dBuffer[3]=param->tolShake;
    dBuffer[4]=param->scal14;
    dBuffer[5]=param->cutOff;
    dBuffer[6]=param->rcutOff;
    dBuffer[7]=param->delr;
    dBuffer[8]=param->cutOff2;
    dBuffer[9]=param->rcutOff2;
    dBuffer[10]=param->cutOn;
    dBuffer[11]=param->cutOn2;
    dBuffer[12]=param->switch2;
    dBuffer[13]=param->temp0;
    dBuffer[14]=param->press0;
    dBuffer[15]=param->kinTemp0;
    dBuffer[16]=param->alpha;
    dBuffer[17]=param->prec;
    dBuffer[18]=param->tol;
    dBuffer[19]=param->damp1;
    dBuffer[20]=param->damp2;
  }
    
  bcast_int_para(iBuffer,12,0);
    
  bcast_double_para(dBuffer,21,0);
    
  if(parallel->idProc>0)
  {
    param->step=iBuffer[0];
    param->nSteps=iBuffer[1];
    param->nAtom=iBuffer[2];
    param->nDegFree=iBuffer[3];
    param->nFrozen=iBuffer[4];
    param->nBond=iBuffer[5];
    param->nAngle=iBuffer[6];
    param->nUb=iBuffer[7];
    param->nDihedral=iBuffer[8];
    param->nImproper=iBuffer[9];
    param->nConst=iBuffer[10];
    param->maxCycle=iBuffer[11];
    
    param->timeStep=dBuffer[0];
    param->rTimeStep=dBuffer[1];
    param->chargeConst=dBuffer[2];
    param->tolShake=dBuffer[3];
    param->scal14=dBuffer[4];
    param->cutOff=dBuffer[5];
    param->rcutOff=dBuffer[6];
    param->delr=dBuffer[7];
    param->cutOff2=dBuffer[8];
    param->rcutOff2=dBuffer[9];
    param->cutOn=dBuffer[10];
    param->cutOn2=dBuffer[11];
    param->switch2=dBuffer[12];
    param->temp0=dBuffer[13];
    param->press0=dBuffer[14];
    param->kinTemp0=dBuffer[15];
    param->alpha=dBuffer[16];
    param->prec=dBuffer[17];
    param->tol=dBuffer[18];
    param->damp1=dBuffer[19];
    param->damp2=dBuffer[20];
  }
  
}

void bcast_ctrl_para(CTRL *ctrl,PARALLEL *parallel,int *iBuffer)
{
  
  if(parallel->idProc==0)
  {
    iBuffer[0]=ctrl->newjob;
    iBuffer[1]=ctrl->keyRand;
    iBuffer[2]=ctrl->keyTraj;
    iBuffer[3]=ctrl->keyProp;
    iBuffer[4]=ctrl->keyForF;
    iBuffer[5]=ctrl->keyRest;
    iBuffer[6]=ctrl->printOut;
    iBuffer[7]=ctrl->printProp;
    iBuffer[8]=ctrl->printTraj;
    iBuffer[9]=ctrl->printRest;
    iBuffer[10]=ctrl->keyMd;
    iBuffer[11]=ctrl->mdType;
    iBuffer[12]=ctrl->seed;
    iBuffer[13]=ctrl->keyMinim;
    iBuffer[14]=ctrl->keyLink;
    iBuffer[15]=ctrl->noLink;
    iBuffer[16]=ctrl->keyConstH;
    iBuffer[17]=ctrl->keyNb14;
    iBuffer[18]=ctrl->keyNumForce;
    iBuffer[19]=ctrl->keyEwald;
    iBuffer[20]=ctrl->keyAlpha;
    iBuffer[21]=ctrl->keyMmax;
    iBuffer[22]=ctrl->elecType;
    iBuffer[23]=ctrl->vdwType;
    iBuffer[24]=ctrl->integrator;
    iBuffer[25]=ctrl->ens;
    iBuffer[26]=ctrl->keyHeuristic;
  }
  
  bcast_int_para(iBuffer,27,0);
  
  if(parallel->idProc>0)
  {
    ctrl->newjob=iBuffer[0];
    ctrl->keyRand=iBuffer[1];
    ctrl->keyTraj=iBuffer[2];
    ctrl->keyProp=iBuffer[3];
    ctrl->keyForF=iBuffer[4];
    ctrl->keyRest=iBuffer[5];
    ctrl->printOut=iBuffer[6];
    ctrl->printProp=iBuffer[7];
    ctrl->printTraj=iBuffer[8];
    ctrl->printRest=iBuffer[9];
    ctrl->keyMd=iBuffer[10];
    ctrl->mdType=iBuffer[11];
    ctrl->seed=iBuffer[12];
    ctrl->keyMinim=iBuffer[13];
    ctrl->keyLink=iBuffer[14];
    ctrl->noLink=iBuffer[15];
    ctrl->keyConstH=iBuffer[16];
    ctrl->keyNb14=iBuffer[17];
    ctrl->keyNumForce=iBuffer[18];
    ctrl->keyEwald=iBuffer[19];
    ctrl->keyAlpha=iBuffer[20];
    ctrl->keyMmax=iBuffer[21];
    ctrl->elecType=iBuffer[22];
    ctrl->vdwType=iBuffer[23];
    ctrl->integrator=iBuffer[24];
    ctrl->ens=iBuffer[25];
    ctrl->keyHeuristic=iBuffer[26];
  }
   
}

void bcast_bath_para(BATH *bath,PARALLEL *parallel,real *dBuffer)
{
  
  if(parallel->idProc==0)
  {
    dBuffer[0]=bath->tauT;
    dBuffer[1]=bath->tauP;
    dBuffer[2]=bath->chiT;
    dBuffer[3]=bath->chiP;
    dBuffer[4]=bath->compress;
  }
  
  bcast_double_para(dBuffer,5,0);
  
  if(parallel->idProc>0)
  {
    bath->tauT=dBuffer[0];
    bath->tauP=dBuffer[1];
    bath->chiT=dBuffer[2];
    bath->chiP=dBuffer[3];
    bath->compress=dBuffer[4];
  }
   
}

void bcast_pbc_para(PBC *box,PARALLEL *parallel,real *dBuffer)
{
  int boxType;
  
  if(parallel->idProc==0)
  {
    boxType=box->type;
    
    dBuffer[0]=box->a1;
    dBuffer[1]=box->a2;
    dBuffer[2]=box->a3;
    
    dBuffer[3]=box->b1;
    dBuffer[4]=box->b2;
    dBuffer[5]=box->b3;
    
    dBuffer[6]=box->c1;
    dBuffer[7]=box->c2;
    dBuffer[8]=box->c3;
    
  }
  
  bcast_int_para(&boxType,1,0);
  bcast_double_para(dBuffer,9,0);
  
  if(parallel->idProc>0)
  {
    box->type=boxType;
    
    box->a1=dBuffer[0];
    box->a2=dBuffer[1];
    box->a3=dBuffer[2];
    
    box->b1=dBuffer[3];
    box->b2=dBuffer[4];
    box->b3=dBuffer[5];
    
    box->c1=dBuffer[6];
    box->c2=dBuffer[7];
    box->c3=dBuffer[8];
  }
   
}

void bcast_neigh_para(NEIGH *neigh,PARALLEL *parallel,int *iBuffer)
{
  
  if(parallel->idProc==0)
  {
    iBuffer[0]=neigh->update;
    iBuffer[1]=neigh->linkRatio;
  }
  
  bcast_int_para(iBuffer,2,0);
  
  if(parallel->idProc>0)
  {
    neigh->update=iBuffer[0];
    neigh->linkRatio=iBuffer[1];
  }
   
}

void bcast_ewald_para(EWALD *ewald,PARALLEL *parallel,int *iBuffer,real *dBuffer)
{
  
  if(parallel->idProc==0)
  {
    iBuffer[0]=ewald->nbsp;
    iBuffer[1]=ewald->mmax;
    iBuffer[2]=ewald->m1max;
    iBuffer[3]=ewald->m2max;
    iBuffer[4]=ewald->m3max;
    
    dBuffer[0]=ewald->prec;
    dBuffer[1]=ewald->alpha;
  }
  
  bcast_int_para(iBuffer,5,0);
  bcast_double_para(dBuffer,2,0);
  
  if(parallel->idProc>0)
  {
    ewald->nbsp=iBuffer[0];
    ewald->mmax=iBuffer[1];
    ewald->m1max=iBuffer[2];
    ewald->m2max=iBuffer[3];
    ewald->m3max=iBuffer[4];
    
    ewald->prec=dBuffer[0];
    ewald->alpha=dBuffer[1];
  }
  
}

void bcast_const_para(CONSTRAINT *constList,PARALLEL *parallel,int *iBuffer,real *dBuffer,int size)
{
  int idx1,idx2;
  
  if(parallel->idProc==0)
  {
    idx1=0;
    idx2=0;
    for(int i=0;i<size;i++)
    {
      iBuffer[idx1++]=constList[i].a;
      iBuffer[idx1++]=constList[i].b;
      
      dBuffer[idx2++]=constList[i].rc2;
    }
  }
  
  bcast_int_para(iBuffer,2*size,0);
  bcast_double_para(dBuffer,size,0);
  
  if(parallel->idProc>0)
  {
    idx1=0;
    idx2=0;
    for(int i=0;i<size;i++)
    {
      constList[i].a=iBuffer[idx1++];
      constList[i].b=iBuffer[idx1++];
      
      constList[i].rc2=dBuffer[idx2++];
    }
  }
  
}

void bcast_bond_para(BOND *bond,PARALLEL *parallel,int *iBuffer,real *dBuffer,int size)
{
  int idx;
  
  if(parallel->idProc==0)
  {
    idx=0;
    for(int i=0;i<size;i++)
    {
      iBuffer[idx]=bond[i].type;
      dBuffer[idx++]=bond[i].k;
      
      iBuffer[idx]=bond[i].a;
      dBuffer[idx++]=bond[i].r0;
      
      iBuffer[idx]=bond[i].b;
      dBuffer[idx++]=bond[i].beta;
    }
  }
  
  bcast_int_para(iBuffer,3*size,0);
  bcast_double_para(dBuffer,3*size,0);
  
  if(parallel->idProc>0)
  {
    idx=0;
    for(int i=0;i<size;i++)
    {
      bond[i].type=iBuffer[idx];
      bond[i].k=dBuffer[idx++];
      
      bond[i].a=iBuffer[idx];
      bond[i].r0=dBuffer[idx++];
      
      bond[i].b=iBuffer[idx];
      bond[i].beta=dBuffer[idx++];
    }
  }
  
}

void bcast_angle_para(ANGLE *angle,PARALLEL *parallel,int *iBuffer,real *dBuffer,int size)
{
  int idx1,idx2;
  
  if(parallel->idProc==0)
  {
    idx1=0;
    idx2=0;
    for(int i=0;i<size;i++)
    {
      iBuffer[idx1++]=angle[i].a;
      iBuffer[idx1++]=angle[i].b;
      iBuffer[idx1++]=angle[i].c;
      
      dBuffer[idx2++]=angle[i].k;
      dBuffer[idx2++]=angle[i].theta0;
    }
  }
  
  bcast_int_para(iBuffer,3*size,0);
  bcast_double_para(dBuffer,2*size,0);
  
  if(parallel->idProc>0)
  {
    idx1=0;
    idx2=0;
    for(int i=0;i<size;i++)
    {
      angle[i].a=iBuffer[idx1++];
      angle[i].b=iBuffer[idx1++];
      angle[i].c=iBuffer[idx1++];
      
      angle[i].k=dBuffer[idx2++];
      angle[i].theta0=dBuffer[idx2++];
    }
  }
  
}

void bcast_dihe_para(DIHE *dihe,PARALLEL *parallel,int *iBuffer,real *dBuffer,int size)
{
  int idx1,idx2;
  
  if(parallel->idProc==0)
  {
    idx1=0;
    idx2=0;
    for(int i=0;i<size;i++)
    {
      iBuffer[idx1++]=dihe[i].type;
      iBuffer[idx1++]=dihe[i].order;
      iBuffer[idx1++]=dihe[i].a;
      iBuffer[idx1++]=dihe[i].b;
      iBuffer[idx1++]=dihe[i].c;
      iBuffer[idx1++]=dihe[i].d;
      
      dBuffer[idx2++]=dihe[i].k;
      dBuffer[idx2++]=dihe[i].phi0;
      dBuffer[idx2++]=dihe[i].mult;
    }
  }
  
  bcast_int_para(iBuffer,6*size,0);
  bcast_double_para(dBuffer,3*size,0);
  
  if(parallel->idProc>0)
  {
    idx1=0;
    idx2=0;
    for(int i=0;i<size;i++)
    {
      dihe[i].type=iBuffer[idx1++];
      dihe[i].order=iBuffer[idx1++];
      dihe[i].a=iBuffer[idx1++];
      dihe[i].b=iBuffer[idx1++];
      dihe[i].c=iBuffer[idx1++];
      dihe[i].d=iBuffer[idx1++];
      
      dihe[i].k=dBuffer[idx2++];
      dihe[i].phi0=dBuffer[idx2++];
      dihe[i].mult=dBuffer[idx2++];
    }
  }
  
}

void setup_para(CTRL *ctrl,PARAM *param,PARALLEL *parallel,BATH *bath,NEIGH *neigh,
		EWALD *ewald,PBC *box,CONSTRAINT **constList,
		BOND **bond,ANGLE **angle,DIHE **dihe,DIHE **impr,BOND **ub,
		real **x,real **y, real **z,
		real **vx,real **vy,real **vz,
		real **fx,real **fy, real **fz,
		real **mass,real **rmass,real **q,
		real **eps,real **sig,real **eps14,real **sig14,
		int **frozen,int **nAtConst,real **dBuffer,int **iBuffer)
{
  
  if(parallel->nProc>1)
  {
  
    (*iBuffer)=(int*)my_malloc(parallel->iBufferSize*sizeof(**iBuffer));
    (*dBuffer)=(real*)my_malloc(parallel->dBufferSize*sizeof(**dBuffer));
      
    bcast_param_para(param,parallel,*dBuffer,*iBuffer);
      
    parallel_allocate_buffers(param,parallel,dBuffer,iBuffer);
      
    bcast_ctrl_para(ctrl,parallel,*iBuffer);
      
    bcast_bath_para(bath,parallel,*dBuffer);
      
    bcast_pbc_para(box,parallel,*dBuffer);
      
    bcast_neigh_para(neigh,parallel,*iBuffer);
      
    bcast_ewald_para(ewald,parallel,*iBuffer,*dBuffer);
  
  }
    
  parallel->maxAtProc=(param->nAtom     + parallel->nProc-1)/parallel->nProc;
  parallel->maxCtProc=(param->nConst    + parallel->nProc-1)/parallel->nProc;
  parallel->maxBdProc=(param->nBond     + parallel->nProc-1)/parallel->nProc;
  parallel->maxAgProc=(param->nAngle    + parallel->nProc-1)/parallel->nProc;
  parallel->maxUbProc=(param->nUb       + parallel->nProc-1)/parallel->nProc;
  parallel->maxDiProc=(param->nDihedral + parallel->nProc-1)/parallel->nProc;
  parallel->maxImProc=(param->nImproper + parallel->nProc-1)/parallel->nProc;
  
  parallel->fAtProc=(parallel->idProc*param->nAtom)/parallel->nProc;
  parallel->lAtProc=((parallel->idProc+1)*param->nAtom)/parallel->nProc;
  parallel->nAtProc=parallel->lAtProc-parallel->fAtProc;
  
  parallel->fCtProc=(parallel->idProc*param->nConst)/parallel->nProc;
  parallel->lCtProc=((parallel->idProc+1)*param->nConst)/parallel->nProc;
  parallel->nCtProc=parallel->lCtProc-parallel->fCtProc;
  
  parallel->fBdProc=(parallel->idProc*param->nBond)/parallel->nProc;
  parallel->lBdProc=((parallel->idProc+1)*param->nBond)/parallel->nProc;
  parallel->nBdProc=parallel->lBdProc-parallel->fBdProc;
  
  parallel->fAgProc=(parallel->idProc*param->nAngle)/parallel->nProc;
  parallel->lAgProc=((parallel->idProc+1)*param->nAngle)/parallel->nProc;
  parallel->nAgProc=parallel->lAgProc-parallel->fAgProc;
  
  parallel->fUbProc=(parallel->idProc*param->nUb)/parallel->nProc;
  parallel->lUbProc=((parallel->idProc+1)*param->nUb)/parallel->nProc;
  parallel->nUbProc=parallel->lUbProc-parallel->fUbProc;
  
  parallel->fDhProc=(parallel->idProc*param->nDihedral)/parallel->nProc;
  parallel->lDhProc=((parallel->idProc+1)*param->nDihedral)/parallel->nProc;
  parallel->nDhProc=parallel->lDhProc-parallel->fDhProc;
  
  parallel->fIpProc=(parallel->idProc*param->nImproper)/parallel->nProc;
  parallel->lIpProc=((parallel->idProc+1)*param->nImproper)/parallel->nProc;
  parallel->nIpProc=parallel->lIpProc-parallel->fIpProc;
  
  if(parallel->idProc>0)
  {
       
    *x=(real*)my_malloc(param->nAtom*sizeof(real));
    *y=(real*)my_malloc(param->nAtom*sizeof(real));
    *z=(real*)my_malloc(param->nAtom*sizeof(real));
    
    *vx=(real*)my_malloc(param->nAtom*sizeof(real));
    *vy=(real*)my_malloc(param->nAtom*sizeof(real));
    *vz=(real*)my_malloc(param->nAtom*sizeof(real));
    
    *fx=(real*)my_malloc(param->nAtom*sizeof(real));
    *fy=(real*)my_malloc(param->nAtom*sizeof(real));
    *fz=(real*)my_malloc(param->nAtom*sizeof(real));
    
    *q=(real*)my_malloc(param->nAtom*sizeof(real));
    
    *mass=(real*)my_malloc(param->nAtom*sizeof(real));
    *rmass=(real*)my_malloc(param->nAtom*sizeof(real));
    
    *frozen=(int*)my_malloc(param->nAtom*sizeof(int));
    
    *nAtConst=(int*)my_malloc(param->nAtom*sizeof(int));
    
    *eps=(real*)my_malloc(param->nAtom*sizeof(real));
    *sig=(real*)my_malloc(param->nAtom*sizeof(real));
    *eps14=(real*)my_malloc(param->nAtom*sizeof(real));
    *sig14=(real*)my_malloc(param->nAtom*sizeof(real));
    
    if(param->nConst>0)
      *constList=(CONSTRAINT*)my_malloc(param->nConst*sizeof(CONSTRAINT));
    
    if(param->nBond>0)
      *bond=(BOND*)my_malloc(param->nBond*sizeof(BOND));
    
    if(param->nUb>0)
      *ub=(BOND*)my_malloc(param->nUb*sizeof(BOND));
    
    if(param->nAngle>0)
      *angle=(ANGLE*)my_malloc(param->nAngle*sizeof(ANGLE));
    
    if(param->nDihedral>0)
      *dihe=(DIHE*)my_malloc(param->nDihedral*sizeof(DIHE));
    
    if(param->nImproper>0)
      *impr=(DIHE*)my_malloc(param->nImproper*sizeof(DIHE));
  }
  
  if(parallel->nProc>1)
  {
  
    bcast_double_para(*x,param->nAtom,0);
    bcast_double_para(*y,param->nAtom,0);
    bcast_double_para(*z,param->nAtom,0);
    
    bcast_double_para(*vx,param->nAtom,0);
    bcast_double_para(*vy,param->nAtom,0);
    bcast_double_para(*vz,param->nAtom,0);
    
    bcast_double_para(*q,param->nAtom,0);
    
    bcast_double_para(*mass,param->nAtom,0);
    bcast_double_para(*rmass,param->nAtom,0);
    
    bcast_double_para(*eps,param->nAtom,0);
    bcast_double_para(*sig,param->nAtom,0);
    bcast_double_para(*eps14,param->nAtom,0);
    bcast_double_para(*sig14,param->nAtom,0);
    
    bcast_int_para(*frozen,param->nAtom,0);
    bcast_int_para(*nAtConst,param->nAtom,0);
    
    bcast_const_para(*constList,parallel,*iBuffer,*dBuffer,param->nConst);
    
    bcast_bond_para(*bond,parallel,*iBuffer,*dBuffer,param->nBond);
    
    bcast_angle_para(*angle,parallel,*iBuffer,*dBuffer,param->nAngle);
    
    bcast_bond_para(*ub,parallel,*iBuffer,*dBuffer,param->nUb);
    
    bcast_dihe_para(*dihe,parallel,*iBuffer,*dBuffer,param->nDihedral);
    
    bcast_dihe_para(*impr,parallel,*iBuffer,*dBuffer,param->nImproper);
  
  }
  
}

void close_para()
{
  int err;
  
  err=MPI_Finalize();
  
}

void mpi_error(int err, char file[],int line)
{
  char errString[STRBUFSIZ];
  int idProc,lenErrString, errClass;
  
  idProc=my_proc();
  
  MPI_Error_class(err, &errClass);
  
  MPI_Error_string(errClass, errString, &lenErrString);
  fprintf(outFile, "%3d: %s\n", idProc, errString);
  
  MPI_Error_string(err, errString, &lenErrString);
  fprintf(outFile, "%3d: %s\n", idProc, errString);
  
  //MPI_Abort(MPI_COMM_WORLD, err);
  
  my_error(MPI_ERROR,file,line,0);
    
}

void abort_para(int err)
{
  
  if(err==0)
  {
    MPI_Abort(MPI_COMM_WORLD,MPI_SUCCESS);
  }
  else
  {
    MPI_Abort(MPI_COMM_WORLD,MPI_ERR_UNKNOWN);
  }
  
  exit(err);
}

