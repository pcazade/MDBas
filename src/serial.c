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
#include "global.h"

void init_para(int *argc, char ***argv)
{
  
}

int my_proc()
{
  
  return (0);
  
}

int num_proc()
{
  
  return (1);
  
}

void barrier_para()
{
  
}

void sum_double_para(double *buf1,double *buf2,int size)
{
  
}

void sum_int_para(int *buf1,int *buf2,int size)
{
  
}

void update_double_para(PARAM *param,PARALLEL *parallel,double *buf1,double *buf2)
{
  
}

void test_para(int *buf1)
{
  
}

void parallel_allocate_buffers(PARAM *param,PARALLEL *parallel,double **dBuffer,
			       int **iBuffer)
{
  
}

void parallel_reallocate_buffers(PARALLEL *parallel,EWALD *ewald,double **dBuffer)
{
  
}

void bcast_int_para(int *buf1,int size,int iNode)
{
  
}

void bcast_double_para(double *buf1,int size,int iNode)
{

}

void setup_para(CTRL *ctrl,PARAM *param,PARALLEL *parallel,ENERGY *ener,
		BATH *bath,NEIGH *neigh,EWALD *ewald,PBC *box,ATOM **atom,CONSTRAINT **constList,
		BOND **bond,ANGLE **angle,DIHE **dihe,DIHE **impr,BOND **ub,double **x,
		double **y, double **z,double **vx,double **vy,double **vz,double **fx,
		double **fy, double **fz,double **mass,double **rmass,double **q,
		double **eps,double **sig,double **eps14,double **sig14,int **frozen,
		int **nAtConst,double **dBuffer,int **iBuffer)
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
  
}

void close_para()
{
  
}

