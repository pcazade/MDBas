/**
 * \file shake.h
 * \brief Prototypes for file shake.c
 * \author Pierre-Andre Cazade and Florent Hedin
 * \version alpha-branch
 * \date 2012
 */

#ifndef SHAKEH_INCLUDED
#define SHAKEH_INCLUDED

void lf_shake  (ATOM atom[],SIMULPARAMS *simulCond,CONSTRAINT *constList,DELTA *dd,PBC *box,double *virshake,double *stress);
void vv_shake_r(ATOM atom[],SIMULPARAMS *simulCond,CONSTRAINT *constList,DELTA *dd,PBC *box,double *virshake,double *stress);
void vv_shake_v(ATOM atom[],SIMULPARAMS *simulCond,CONSTRAINT *constList,DELTA *dd);

#endif