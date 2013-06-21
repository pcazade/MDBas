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

#ifdef MPI_VERSION

#include <stdio.h>
#include <mpi.h>

#include "global.h"
#include "errors.h"
#include "parallel.h"

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
  
  int parallel->idProc,err;
  
  err=MPI_Comm_rank(MPI_COMM_WORLD,&parallel->idProc);
  
  if(err!=MPI_SUCCESS)
    mpi_error(err,__FILE__,__LINE__);
  
  return (parallel->idProc);
  
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
  
  for(i=0;i<parallel->fAtom;i++)
  {
    buf2[i]=0.;
  }
  
  for(i=parallel->fAtom;i<parallel->lAtom;i++)
  {
    buf2[i]=buf1[i];
  }
  
  for(i=parallel->lAtom;i<param->nAtom;i++)
  {
    buf2[i]=0.;
  }
  
  err=MPI_Allreduce(buf2,buf1,param->nAtom,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  
  if(err!=MPI_SUCCESS)
    mpi_error(err,__FILE__,__LINE__);
  
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

#endif