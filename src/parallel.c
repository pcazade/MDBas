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

//#ifdef MPI_VERSION

#include <stdio.h>
#include <mpi.h>

#include "global.h"
#include "errors.h"
#include "parallel.h"
#include "memory.h"

#define BUFSIZ 8192

void init_para(int *argc, char ***argv,PARAM *param)
{
  int err;
    
  MPI_Errhandler_set(MPI_COMM_WORLD, MPI_ERRORS_RETURN);
  
  err=MPI_Init(argc,argv);
  
  if(err!=MPI_SUCCESS)
    mpi_error(err,__FILE__,__LINE__);
  
}

int my_proc()
{
  
  int idProc,err;
  
  err=MPI_Comm_rank(MPI_COMM_WORLD,&parallel->idProc);
  
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

void sum_double_para(double *buf1,double *buf2,int size)
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

void update_double_para(PARAM *param,PARALLEL *parallel,double *buf1,double *buf2)
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
  int test1,test2,err;
  
  test1=0;
  
  if(!*buf1)
    test1=1;
  
  err=MPI_Allreduce(&test1,&test2,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  
  if(err!=MPI_SUCCESS)
    mpi_error(err,__FILE__,__LINE__);
  
  if(test2==0)
    *buf1=1;
  
}

void setup_para(CTRL *ctrl,PARAM *param,PARALLEL *parallel,ENERGY *ener,
		BATH *bath,NEIGH *neigh,EWALD *ewald,PBC *box,ATOM **atom,CONSTRAINT **constList,
		BOND **bond,ANGLE **angle,DIHE **dihe,DIHE **impr,BOND **ub,double **x,
		double **y, double **z,double **vx,double **vy,double **vz,double **fx,
		double **fy, double **fz,double **mass,double **rmass,double **q,
		double **eps,double **sig,double **eps14,double **sig14,int **frozen,
		int **nAtConst,int **neighList,int **neighPair,int **neighList14,
		int ***exclList,int **exclPair)
{
  
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
  
  if(parallel->idProc>0)
  {
       
    *x=(double*)my_malloc(param->nAtom*sizeof(double));
    *y=(double*)my_malloc(param->nAtom*sizeof(double));
    *z=(double*)my_malloc(param->nAtom*sizeof(double));
    
    *vx=(double*)my_malloc(param->nAtom*sizeof(double));
    *vy=(double*)my_malloc(param->nAtom*sizeof(double));
    *vz=(double*)my_malloc(param->nAtom*sizeof(double));
    
    *fx=(double*)my_malloc(param->nAtom*sizeof(double));
    *fy=(double*)my_malloc(param->nAtom*sizeof(double));
    *fz=(double*)my_malloc(param->nAtom*sizeof(double));
    
    *q=(double*)my_malloc(param->nAtom*sizeof(double));
    
    *mass=(double*)my_malloc(param->nAtom*sizeof(double));
    *rmass=(double*)my_malloc(param->nAtom*sizeof(double));
    
    *frozen=(int*)my_malloc(param->nAtom*sizeof(int));
    
    *nAtConst=(int*)my_malloc(param->nAtom*sizeof(int));
    
    *eps=(double*)my_malloc(param->nAtom*sizeof(double));
    *sig=(double*)my_malloc(param->nAtom*sizeof(double));
    *eps14=(double*)my_malloc(param->nAtom*sizeof(double));
    *sig14=(double*)my_malloc(param->nAtom*sizeof(double));
    
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
}

void close_para()
{
  int err;
  
  err=MPI_Finalize();
  
}

void mpi_error(int err, char file[],int line)
{
  char errString[BUFSIZ];
  int parallel->idProc,lenErrString, errClass;
  
  parallel->idProc=myproc();
  
  MPI_Error_class(err, &errClass);
  
  MPI_Error_string(errClass, errString, &lenErrString);
  fprintf(outFile, "%3d: %s\n", parallel->idProc, errString);
  
  MPI_Error_string(err, errString, &lenErrString);
  fprintf(outFile, "%3d: %s\n", parallel->idProc, errString);
  
  MPI_Abort(MPI_COMM_WORLD, err);
  
  my_error(MPI_ERROR,file,line,0);
    
}

//#endif