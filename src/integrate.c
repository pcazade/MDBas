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

/**
 * \file integrate.c
 * \brief Contains functions performing the integration
 * \author Pierre-Andre Cazade and Florent Hedin
 * \version alpha-branch
 * \date 2012
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "global.h"
#include "utils.h"
#include "shake.h"
#include "integrate.h"
#include "memory.h"

#ifdef USING_MPI
#include "parallel.h"
#else
#include "serial.h"
#endif

#if (defined TIMING && defined __unix__ && !defined __STRICT_ANSI__)
#define TIMER
#include "timing.h"
#endif

#ifdef _OPENMP
#undef _OPENMP
#endif

/** Pointer to the output file. **/
extern FILE *outFile;

static double *ddx,*ddy,*ddz;
static double *xo,*yo,*zo;
static double *xt,*yt,*zt;
static double *vxo,*vyo,*vzo;
static double *vxu,*vyu,*vzu;

void integrators_allocate_arrays(CTRL *ctrl,PARAM *param,PARALLEL *parallel)
{
    ddx=ddy=ddz=NULL;
    xo=yo=zo=xt=yt=zt=vxo=vyo=vzo=vxu=vyu=vzu=NULL;

    if(ctrl->integrator == LEAPFROG)
    {
        if (ctrl->ens == NVE)
        {
            xo=(double*)my_malloc(parallel->maxAtProc*sizeof(*xo));
            yo=(double*)my_malloc(parallel->maxAtProc*sizeof(*yo));
            zo=(double*)my_malloc(parallel->maxAtProc*sizeof(*zo));
            vxu=(double*)my_malloc(parallel->maxAtProc*sizeof(*vxu));
            vyu=(double*)my_malloc(parallel->maxAtProc*sizeof(*vyu));
            vzu=(double*)my_malloc(parallel->maxAtProc*sizeof(*vzu));
            if(param->nConst>0)
            {
                ddx=(double*)my_malloc(parallel->maxCtProc*sizeof(*ddx));
                ddy=(double*)my_malloc(parallel->maxCtProc*sizeof(*ddy));
                ddz=(double*)my_malloc(parallel->maxCtProc*sizeof(*ddz));
            }
        }
        else
        {
            xo=(double*)my_malloc(parallel->maxAtProc*sizeof(*xo));
            yo=(double*)my_malloc(parallel->maxAtProc*sizeof(*yo));
            zo=(double*)my_malloc(parallel->maxAtProc*sizeof(*zo));
            vxo=(double*)my_malloc(parallel->maxAtProc*sizeof(*vxo));
            vyo=(double*)my_malloc(parallel->maxAtProc*sizeof(*vyo));
            vzo=(double*)my_malloc(parallel->maxAtProc*sizeof(*vzo));
            vxu=(double*)my_malloc(parallel->maxAtProc*sizeof(*vxu));
            vyu=(double*)my_malloc(parallel->maxAtProc*sizeof(*vyu));
            vzu=(double*)my_malloc(parallel->maxAtProc*sizeof(*vzu));
            if(param->nConst>0)
            {
                ddx=(double*)my_malloc(parallel->maxCtProc*sizeof(*ddx));
                ddy=(double*)my_malloc(parallel->maxCtProc*sizeof(*ddy));
                ddz=(double*)my_malloc(parallel->maxCtProc*sizeof(*ddz));

                xt=(double*)my_malloc(parallel->maxAtProc*sizeof(*xt));
                yt=(double*)my_malloc(parallel->maxAtProc*sizeof(*yt));
                zt=(double*)my_malloc(parallel->maxAtProc*sizeof(*zt));
            }
        }
    }
    else if (ctrl->integrator == VELOCITY)
    {
        if (ctrl->ens == NVE || ctrl->ens == NVT_B || ctrl->ens == NVT_H)
        {
            if(param->nConst>0)
            {
                ddx=(double*)my_malloc(parallel->maxCtProc*sizeof(*ddx));
                ddy=(double*)my_malloc(parallel->maxCtProc*sizeof(*ddy));
                ddz=(double*)my_malloc(parallel->maxCtProc*sizeof(*ddz));
            }
        }
        else
        {
            xo=(double*)my_malloc(parallel->maxAtProc*sizeof(*xo));
            yo=(double*)my_malloc(parallel->maxAtProc*sizeof(*yo));
            zo=(double*)my_malloc(parallel->maxAtProc*sizeof(*zo));
            vxo=(double*)my_malloc(parallel->maxAtProc*sizeof(*vxo));
            vyo=(double*)my_malloc(parallel->maxAtProc*sizeof(*vyo));
            vzo=(double*)my_malloc(parallel->maxAtProc*sizeof(*vzo));
            if(param->nConst>0)
            {
                ddx=(double*)my_malloc(parallel->maxCtProc*sizeof(*ddx));
                ddy=(double*)my_malloc(parallel->maxCtProc*sizeof(*ddy));
                ddz=(double*)my_malloc(parallel->maxCtProc*sizeof(*ddz));
            }
        }
    }
}

void integrators_free_arrays(CTRL *ctrl,PARAM *param)
{
    if(ctrl->integrator == LEAPFROG)
    {
        if (ctrl->ens == NVE)
        {
            free(xo);
            free(yo);
            free(zo);
            free(vxu);
            free(vyu);
            free(vzu);
            if(param->nConst>0)
            {
                free(ddx);
                free(ddy);
                free(ddz);
            }
        }
        else
        {
            free(xo);
            free(yo);
            free(zo);
            free(vxo);
            free(vyo);
            free(vzo);
            free(vxu);
            free(vyu);
            free(vzu);
            if(param->nConst>0)
            {
                free(ddx);
                free(ddy);
                free(ddz);

                free(xt);
                free(yt);
                free(zt);
            }
        }
    }
    else if (ctrl->integrator == VELOCITY)
    {
        if (ctrl->ens == NVE || ctrl->ens == NVT_B || ctrl->ens == NVT_H)
        {
            if(param->nConst>0)
            {
                free(ddx);
                free(ddy);
                free(ddz);
            }
        }
        else
        {
            free(xo);
            free(yo);
            free(zo);
            free(vxo);
            free(vyo);
            free(vzo);
            if(param->nConst>0)
            {
                free(ddx);
                free(ddy);
                free(ddz);
            }
        }
    }

}

void lf_integrate(CTRL *ctrl,PARAM *param,ENERGY *ener,PBC *box,
                  BATH *bath,CONSTRAINT constList[],PARALLEL *parallel,
                  double *x,double *y,double *z,
                  double *vx,double *vy,double *vz,
                  double *fx,double *fy,double *fz,
                  double *mass,double *rmass,int *nAtConst,double dBuffer[])
{

#ifdef TIMER
    update_timer_begin(TIMER_INTEGRATE,__func__);
#endif

    switch (ctrl->ens)
    {
    case NVE:
        lf_nve(param,ener,box,constList,parallel,x,y,z,vx,vy,vz,fx,fy,fz,mass,rmass,nAtConst,dBuffer);
        break;
    case NVT_B:
        lf_nvt_b(param,ener,box,bath,constList,parallel,x,y,z,vx,vy,vz,fx,fy,fz,mass,rmass,nAtConst,dBuffer);
        break;
    case NPT_B:
        lf_npt_b(param,ener,box,bath,constList,parallel,x,y,z,vx,vy,vz,fx,fy,fz,mass,rmass,nAtConst,dBuffer);
        break;
    case NVT_H:
        lf_nvt_h(param,ener,box,bath,constList,parallel,x,y,z,vx,vy,vz,fx,fy,fz,mass,rmass,nAtConst,dBuffer);
        break;
    case NPT_H:
        lf_npt_h(param,ener,box,bath,constList,parallel,x,y,z,vx,vy,vz,fx,fy,fz,mass,rmass,nAtConst,dBuffer);
        break;
    default:
        lf_nve(param,ener,box,constList,parallel,x,y,z,vx,vy,vz,fx,fy,fz,mass,rmass,nAtConst,dBuffer);
        break;
    }

#ifdef TIMER
    update_timer_end(TIMER_INTEGRATE,__func__);
#endif

}

void lf_nve(PARAM *param,ENERGY *ener,PBC *box,
            CONSTRAINT constList[],PARALLEL *parallel,
            double *x,double *y,double *z,
            double *vx,double *vy,double *vz,
            double *fx,double *fy,double *fz,
            double *mass,double *rmass,int *nAtConst,double dBuffer[])
{
    int i,ia,ib,l;
    double virshake=0.,stress[6]= {0.},stresk[6]= {0.};

    if(param->nConst>0)
    {
        l=0;
#ifdef _OPENMP
        #pragma omp parallel for default(none) shared(param,atom,constList,dd) private(i,ia,ib)
#endif
        for(i=parallel->fCtProc; i<parallel->lCtProc; i++)
        {
            ia=constList[i].a;
            ib=constList[i].b;

            ddx[l]=x[ib]-x[ia];
            ddy[l]=y[ib]-y[ia];
            ddz[l]=z[ib]-z[ia];

            l++;
        }

        image_array(box,ddx,ddy,ddz,parallel->nCtProc);

        l=0;
#ifdef _OPENMP
        #pragma omp parallel for default(none) shared(param,xo,yo,zo,atom) private(i)
#endif
        for(i=parallel->fAtProc; i<parallel->lAtProc; i++)
        {

// Store old coordinates.

            xo[l]=x[i];
            yo[l]=y[i];
            zo[l]=z[i];

            l++;

        }

    }

// move atoms by leapfrog algorithm
    l=0;
#ifdef _OPENMP
    #pragma omp parallel for default(none) shared(vxu,vyu,vzu,param,atom) private(i)
#endif
    for(i=parallel->fAtProc; i<parallel->lAtProc; i++)
    {

// update velocities

        vxu[l]=vx[i]+param->timeStep*fx[i]*rmass[i];
        vyu[l]=vy[i]+param->timeStep*fy[i]*rmass[i];
        vzu[l]=vz[i]+param->timeStep*fz[i]*rmass[i];

// update positions

        x[i]+=param->timeStep*vxu[l];
        y[i]+=param->timeStep*vyu[l];
        z[i]+=param->timeStep*vzu[l];

        l++;
    }

    if(param->nConst>0)
    {

        if(parallel->nProc>1)
        {
            update_double_para(param,parallel,x,dBuffer);
            update_double_para(param,parallel,y,dBuffer);
            update_double_para(param,parallel,z,dBuffer);
        }

// Apply constraint with Shake algorithm.

        lf_shake(param,box,constList,parallel,x,y,z,ddx,ddy,ddz,rmass,nAtConst,stress,&virshake,dBuffer);

        l=0;
#ifdef _OPENMP
        #pragma omp parallel for default(none) shared(vxu,vyu,vzu,xo,yo,zo,param,atom) private(i)
#endif
        for(i=parallel->fAtProc; i<parallel->lAtProc; i++)
        {

// Corrected velocities

            vxu[l]=(x[i]-xo[l])*param->rTimeStep;
            vyu[l]=(y[i]-yo[l])*param->rTimeStep;
            vzu[l]=(z[i]-zo[l])*param->rTimeStep;

// Corrected Forces

            fx[i]=(vxu[l]-vx[i])*mass[i]*param->rTimeStep;
            fy[i]=(vyu[l]-vy[i])*mass[i]*param->rTimeStep;
            fz[i]=(vzu[l]-vz[i])*mass[i]*param->rTimeStep;

            l++;

        }
    }

// calculate full timestep velocity
    l=0;
#ifdef _OPENMP
    #pragma omp parallel for default(none) shared(vxu,vyu,vzu,param,atom) private(i)
#endif
    for(i=parallel->fAtProc; i<parallel->lAtProc; i++)
    {

        vx[i]=0.5*(vx[i]+vxu[l]);
        vy[i]=0.5*(vy[i]+vyu[l]);
        vz[i]=0.5*(vz[i]+vzu[l]);

        l++;

    }

// calculate kinetic energy

    ener->kin=kinetic(parallel,vx,vy,vz,mass,dBuffer);

    ener->virshake=virshake;

//   stress_kinetic(atom,simulCond,stresk);
    stress_kinetic(parallel,vx,vy,vz,mass,stresk,dBuffer);

    box->stress1+=stress[0]+stresk[0];
    box->stress2+=stress[1]+stresk[1];
    box->stress3+=stress[2]+stresk[2];
    box->stress4+=stress[1]+stresk[1];
    box->stress5+=stress[3]+stresk[3];
    box->stress6+=stress[4]+stresk[4];
    box->stress7+=stress[2]+stresk[2];
    box->stress8+=stress[4]+stresk[4];
    box->stress9+=stress[5]+stresk[5];

// periodic boundary condition

//   image_update(atom,simulCond,box);
    image_update(parallel,box,x,y,z);

// updated velocity
    l=0;
#ifdef _OPENMP
    #pragma omp parallel for default(none) shared(vxu,vyu,vzu,param,atom) private(i)
#endif
    for(i=parallel->fAtProc; i<parallel->lAtProc; i++)
    {

        vx[i]=vxu[l];
        vy[i]=vyu[l];
        vz[i]=vzu[l];

        l++;

    }

    if(parallel->nProc>1)
    {
        update_double_para(param,parallel,x,dBuffer);
        update_double_para(param,parallel,y,dBuffer);
        update_double_para(param,parallel,z,dBuffer);

        update_double_para(param,parallel,vx,dBuffer);
        update_double_para(param,parallel,vy,dBuffer);
        update_double_para(param,parallel,vz,dBuffer);

        if(param->nConst>0)
        {
            update_double_para(param,parallel,fx,dBuffer);
            update_double_para(param,parallel,fy,dBuffer);
            update_double_para(param,parallel,fz,dBuffer);
        }

    }
}

void lf_nvt_b(PARAM *param,ENERGY *ener,PBC *box,BATH *bath,CONSTRAINT constList[],PARALLEL *parallel,
              double *x,double *y,double *z,
              double *vx,double *vy,double *vz,
              double *fx,double *fy,double *fz,
              double *mass,double *rmass,int *nAtConst,double dBuffer[])
{
    int i,k,l,ia,ib,bercycle;
    double lambda,rts2;
    double virshake=0.,virshakt=0.,stress[6]= {0.},strest[6]= {0.},stresk[6]= {0.};

    l=0;
#ifdef _OPENMP
    #pragma omp parallel for default(none) shared(xo,yo,zo,vxo,vyo,vzo,param,atom) private(i)
#endif
    for(i=parallel->fAtProc; i<parallel->lAtProc; i++)
    {

// Store old coordinates and old velocities.

        xo[l]=x[i];
        yo[l]=y[i];
        zo[l]=z[i];

        vxo[l]=vx[i];
        vyo[l]=vy[i];
        vzo[l]=vz[i];

        l++;

    }

    if(param->nConst>0)
    {
        l=0;
#ifdef _OPENMP
        #pragma omp parallel for default(none) shared(param,atom,constList,dd) private(i,ia,ib)
#endif
        for(i=parallel->fCtProc; i<parallel->lCtProc; i++)
        {
            ia=constList[i].a;
            ib=constList[i].b;

            ddx[l]=x[ib]-x[ia];
            ddy[l]=y[ib]-y[ia];
            ddz[l]=z[ib]-z[ia];
	    
	    l++;
        }

        image_array(box,ddx,ddy,ddz,parallel->nCtProc);

    }

    rts2=1./X2(param->timeStep);

#ifdef _OPENMP
    #pragma omp parallel for default(none) shared(param,atom) private(i)
#endif
    for(i=parallel->fAtProc; i<parallel->lAtProc; i++)
    {
        vx[i]+=0.5*param->timeStep*fx[i]*rmass[i];
        vy[i]+=0.5*param->timeStep*fy[i]*rmass[i];
        vz[i]+=0.5*param->timeStep*fz[i]*rmass[i];
    }

    ener->kin=kinetic(parallel,vx,vy,vz,mass,dBuffer);

    if(param->nConst>0)
        bercycle=2;
    else
        bercycle=3;

    for(k=0; k<bercycle; k++)
    {

        lambda=sqrt(1.0+param->timeStep/bath->tauT*(param->kinTemp0/ener->kin-1.0));

// move atoms by leapfrog algorithm
        l=0;
#ifdef _OPENMP
        #pragma omp parallel for default(none) shared(param,atom,xo,yo,zo,xt,yt,zt,vxo,vyo,vzo,vxu,vyu,vzu,lambda) private(i)
#endif
        for(i=parallel->fAtProc; i<parallel->lAtProc; i++)
        {

// update velocities

            vxu[l]=(vxo[l]+param->timeStep*fx[i]*rmass[i])*lambda;
            vyu[l]=(vyo[l]+param->timeStep*fy[i]*rmass[i])*lambda;
            vzu[l]=(vzo[l]+param->timeStep*fz[i]*rmass[i])*lambda;

// update positions

            x[i]=xo[l]+param->timeStep*vxu[l];
            y[i]=yo[l]+param->timeStep*vyu[l];
            z[i]=zo[l]+param->timeStep*vzu[l];

// Temporary storage of the uncorrected positions

            if(param->nConst>0)
            {
                xt[l]=x[i];
                yt[l]=y[i];
                zt[l]=z[i];
            }

            l++;
        }

        if( (param->nConst>0) && (k==0) )
        {

            if(parallel->nProc>1)
            {
                update_double_para(param,parallel,x,dBuffer);
                update_double_para(param,parallel,y,dBuffer);
                update_double_para(param,parallel,z,dBuffer);
            }

// Apply constraint with Shake algorithm.

            lf_shake(param,box,constList,parallel,x,y,z,ddx,ddy,ddz,rmass,nAtConst,strest,&virshakt,dBuffer);

            virshake+=virshakt;
            for(i=0; i<6; i++)
                stress[i]+=strest[i];

            l=0;
#ifdef _OPENMP
            #pragma omp parallel for default(none) shared(param,atom,xt,yt,zt,vxu,vyu,vzu,rts2) private(i)
#endif
            for(i=parallel->fAtProc; i<parallel->lAtProc; i++)
            {

// Corrected velocities

                vxu[l]+=(x[i]-xt[l])*param->rTimeStep;
                vyu[l]+=(y[i]-yt[l])*param->rTimeStep;
                vzu[l]+=(z[i]-zt[l])*param->rTimeStep;

// Corrected Forces

                fx[i]+=(x[i]-xt[l])*mass[i]*rts2;
                fy[i]+=(y[i]-yt[l])*mass[i]*rts2;
                fz[i]+=(z[i]-zt[l])*mass[i]*rts2;

                l++;

            }
        }

// calculate full timestep velocity
        l=0;
#ifdef _OPENMP
        #pragma omp parallel for default(none) shared(vxo,vyo,vzo,vxu,vyu,vzu,param,atom) private(i)
#endif
        for(i=parallel->fAtProc; i<parallel->lAtProc; i++)
        {

            vx[i]=0.5*(vxo[l]+vxu[l]);
            vy[i]=0.5*(vyo[l]+vyu[l]);
            vz[i]=0.5*(vzo[l]+vzu[l]);

            l++;
        }

// calculate kinetic energy

        ener->kin=kinetic(parallel,vx,vy,vz,mass,dBuffer);

    }

    ener->virshake=virshake;

//   stress_kinetic(atom,simulCond,stresk);
    stress_kinetic(parallel,vx,vy,vz,mass,stresk,dBuffer);

    box->stress1+=stress[0]+stresk[0];
    box->stress2+=stress[1]+stresk[1];
    box->stress3+=stress[2]+stresk[2];
    box->stress4+=stress[1]+stresk[1];
    box->stress5+=stress[3]+stresk[3];
    box->stress6+=stress[4]+stresk[4];
    box->stress7+=stress[2]+stresk[2];
    box->stress8+=stress[4]+stresk[4];
    box->stress9+=stress[5]+stresk[5];

// periodic boundary condition

//   image_update(atom,simulCond,box);
    image_update(parallel,box,x,y,z);

// updated velocity
    l=0;
#ifdef _OPENMP
    #pragma omp parallel for default(none) shared(vxu,vyu,vzu,param,atom) private(i)
#endif
    for(i=parallel->fAtProc; i<parallel->lAtProc; i++)
    {

        vx[i]=vxu[l];
        vy[i]=vyu[l];
        vz[i]=vzu[l];

        l++;

    }

    if(parallel->nProc>1)
    {
        update_double_para(param,parallel,x,dBuffer);
        update_double_para(param,parallel,y,dBuffer);
        update_double_para(param,parallel,z,dBuffer);

        update_double_para(param,parallel,vx,dBuffer);
        update_double_para(param,parallel,vy,dBuffer);
        update_double_para(param,parallel,vz,dBuffer);

        if(param->nConst>0)
        {
            update_double_para(param,parallel,fx,dBuffer);
            update_double_para(param,parallel,fy,dBuffer);
            update_double_para(param,parallel,fz,dBuffer);
        }

    }

}

void lf_npt_b(PARAM *param,ENERGY *ener,PBC *box,BATH *bath,CONSTRAINT constList[],PARALLEL *parallel,
              double *x,double *y,double *z,
              double *vx,double *vy,double *vz,
              double *fx,double *fy,double *fz,
              double *mass,double *rmass,int *nAtConst,double dBuffer[])
{
    int i,k,l,ia,ib,bercycle;
    double lambda,gamma,cbrga,pp,rts2;
    double volume,cell0[9];
    double virshake=0.,virshakt=0.,stress[6]= {0.},strest[6]= {0.},stresk[6]= {0.};

    l=0;
#ifdef _OPENMP
    #pragma omp parallel for default(none) shared(xo,yo,zo,vxo,vyo,vzo,param,atom) private(i)
#endif
    for(i=parallel->fAtProc; i<parallel->lAtProc; i++)
    {

// Store old coordinates and old velocities.

        xo[l]=x[i];
        yo[l]=y[i];
        zo[l]=z[i];

        vxo[l]=vx[i];
        vyo[l]=vy[i];
        vzo[l]=vz[i];

        l++;

    }

    //Store initial box parameters

    cell0[0]=box->a1;
    cell0[1]=box->a2;
    cell0[2]=box->a3;
    cell0[3]=box->b1;
    cell0[4]=box->b2;
    cell0[5]=box->b3;
    cell0[6]=box->c1;
    cell0[7]=box->c2;
    cell0[8]=box->c3;

    volume=box->vol;

    if(param->nConst>0)
    {
        l=0;
#ifdef _OPENMP
        #pragma omp parallel for default(none) shared(param,atom,constList,dd) private(i,ia,ib)
#endif
        for(i=parallel->fCtProc; i<parallel->lCtProc; i++)
        {
            ia=constList[i].a;
            ib=constList[i].b;

            ddx[l]=x[ib]-x[ia];
            ddy[l]=y[ib]-y[ia];
            ddz[l]=z[ib]-z[ia];

            l++;
        }

        image_array(box,ddx,ddy,ddz,parallel->nCtProc);

    }

    rts2=1./X2(param->timeStep);

#ifdef _OPENMP
    #pragma omp parallel for default(none) shared(param,atom) private(i)
#endif
    for(i=parallel->fAtProc; i<parallel->lAtProc; i++)
    {
        vx[i]+=0.5*param->timeStep*fx[i]*rmass[i];
        vy[i]+=0.5*param->timeStep*fy[i]*rmass[i];
        vz[i]+=0.5*param->timeStep*fz[i]*rmass[i];
    }

    ener->kin=kinetic(parallel,vx,vy,vz,mass,dBuffer);

    pp=(2.*ener->kin-ener->virpot-virshake)/(3.*volume);
    gamma=1.+bath->compress*param->timeStep*(pp-param->press0)/bath->tauP;
    cbrga=cbrt(gamma);

    lambda=sqrt(1.0+param->timeStep/bath->tauT*(param->kinTemp0/ener->kin-1.0));

    if(param->nConst>0)
        bercycle=4;
    else
        bercycle=5;

    for(k=0; k<bercycle; k++)
    {

// move atoms by leapfrog algorithm
        l=0;
#ifdef _OPENMP
        #pragma omp parallel for default(none) shared(param,atom,xo,yo,zo,xt,yt,zt,vxo,vyo,vzo,vxu,vyu,vzu,lambda,cbrga) private(i)
#endif
        for(i=parallel->fAtProc; i<parallel->lAtProc; i++)
        {

// update velocities

            vxu[l]=(vxo[l]+param->timeStep*fx[i]*rmass[i])*lambda;
            vyu[l]=(vyo[l]+param->timeStep*fy[i]*rmass[i])*lambda;
            vzu[l]=(vzo[l]+param->timeStep*fz[i]*rmass[i])*lambda;

// update positions

            x[i]=cbrga*xo[l]+param->timeStep*vxu[l];
            y[i]=cbrga*yo[l]+param->timeStep*vyu[l];
            z[i]=cbrga*zo[l]+param->timeStep*vzu[l];

// Temporary storage of the uncorrected positions

            if(param->nConst>0)
            {
                xt[l]=x[i];
                yt[l]=y[i];
                zt[l]=z[i];
            }

            l++;

        }

        if( (param->nConst>0) && (k==0) )
        {

            if(parallel->nProc>1)
            {
                update_double_para(param,parallel,x,dBuffer);
                update_double_para(param,parallel,y,dBuffer);
                update_double_para(param,parallel,z,dBuffer);
            }

            scale_box(box,cell0,cbrga);

// Apply constraint with Shake algorithm.

            lf_shake(param,box,constList,parallel,x,y,z,ddx,ddy,ddz,rmass,nAtConst,strest,&virshakt,dBuffer);

            virshake+=virshakt;
            for(i=0; i<6; i++)
                stress[i]+=strest[i];

            l=0;
#ifdef _OPENMP
            #pragma omp parallel for default(none) shared(param,atom,xt,yt,zt,vxu,vyu,vzu,rts2) private(i)
#endif
            for(i=parallel->fAtProc; i<parallel->lAtProc; i++)
            {

// Corrected velocities

                vxu[l]+=(x[i]-xt[l])*param->rTimeStep;
                vyu[l]+=(y[i]-yt[l])*param->rTimeStep;
                vzu[l]+=(z[i]-zt[l])*param->rTimeStep;

// Corrected Forces

                fx[i]+=(x[i]-xt[l])*mass[i]*rts2;
                fy[i]+=(y[i]-yt[l])*mass[i]*rts2;
                fz[i]+=(z[i]-zt[l])*mass[i]*rts2;

                l++;
            }
        }

// calculate full timestep velocity

        l=0;
#ifdef _OPENMP
        #pragma omp parallel for default(none) shared(vxo,vyo,vzo,vxu,vyu,vzu,param,atom) private(i)
#endif
        for(i=parallel->fAtProc; i<parallel->lAtProc; i++)
        {

            vx[i]=0.5*(vxo[l]+vxu[l]);
            vy[i]=0.5*(vyo[l]+vyu[l]);
            vz[i]=0.5*(vzo[l]+vzu[l]);

            l++;

        }

// calculate kinetic energy

        ener->kin=kinetic(parallel,vx,vy,vz,mass,dBuffer);

        pp=(2.*ener->kin-ener->virpot-virshake)/(3.*volume);
        gamma=1.+watercomp*param->timeStep*(pp-param->press0)/bath->tauP;
        cbrga=cbrt(gamma);

        lambda=sqrt(1.0+param->timeStep/bath->tauT*(param->kinTemp0/ener->kin-1.0));

    }

    ener->virshake=virshake;

//   stress_kinetic(atom,simulCond,stresk);
    stress_kinetic(parallel,vx,vy,vz,mass,stresk,dBuffer);

    box->stress1+=stress[0]+stresk[0];
    box->stress2+=stress[1]+stresk[1];
    box->stress3+=stress[2]+stresk[2];
    box->stress4+=stress[1]+stresk[1];
    box->stress5+=stress[3]+stresk[3];
    box->stress6+=stress[4]+stresk[4];
    box->stress7+=stress[2]+stresk[2];
    box->stress8+=stress[4]+stresk[4];
    box->stress9+=stress[5]+stresk[5];

//   scale_box(box,cbrga,cell0);
    scale_box(box,cell0,cbrga);

// periodic boundary condition

//   image_update(atom,simulCond,box);
    image_update(parallel,box,x,y,z);

// updated velocity
    l=0;
#ifdef _OPENMP
    #pragma omp parallel for default(none) shared(vxu,vyu,vzu,param,atom) private(i)
#endif
    for(i=parallel->fAtProc; i<parallel->lAtProc; i++)
    {

        vx[i]=vxu[l];
        vy[i]=vyu[l];
        vz[i]=vzu[l];

        l++;
    }

    if(parallel->nProc>1)
    {
        update_double_para(param,parallel,x,dBuffer);
        update_double_para(param,parallel,y,dBuffer);
        update_double_para(param,parallel,z,dBuffer);

        update_double_para(param,parallel,vx,dBuffer);
        update_double_para(param,parallel,vy,dBuffer);
        update_double_para(param,parallel,vz,dBuffer);

        if(param->nConst>0)
        {
            update_double_para(param,parallel,fx,dBuffer);
            update_double_para(param,parallel,fy,dBuffer);
            update_double_para(param,parallel,fz,dBuffer);
        }

    }

}

void lf_nvt_h(PARAM *param,ENERGY *ener,PBC *box,BATH *bath,CONSTRAINT constList[],PARALLEL *parallel,
              double *x,double *y,double *z,
              double *vx,double *vy,double *vz,
              double *fx,double *fy,double *fz,
              double *mass,double *rmass,int *nAtConst,double dBuffer[])
{
    int i,k,l,ia,ib,nosecycle;
    double lambda,lambdb,lambdc,rts2,qmass;
    double virshake=0.,virshakt=0.,stress[6]= {0.},strest[6]= {0.},stresk[6]= {0.};

    l=0;
#ifdef _OPENMP
    #pragma omp parallel for default(none) shared(xo,yo,zo,vxo,vyo,vzo,param,atom) private(i)
#endif
    for(i=parallel->fAtProc; i<parallel->lAtProc; i++)
    {

// Store old coordinates and old velocities.

        xo[l]=x[i];
        yo[l]=y[i];
        zo[l]=z[i];

        vxo[l]=vx[i];
        vyo[l]=vy[i];
        vzo[l]=vz[i];

        l++;
    }

    if(param->nConst>0)
    {
        l=0;
#ifdef _OPENMP
        #pragma omp parallel for default(none) shared(param,atom,dd,constList) private(i,ia,ib)
#endif
        for(i=parallel->fCtProc; i<parallel->lCtProc; i++)
        {
            ia=constList[i].a;
            ib=constList[i].b;

            ddx[l]=x[ib]-x[ia];
            ddy[l]=y[ib]-y[ia];
            ddz[l]=z[ib]-z[ia];

            l++;
        }

        image_array(box,ddx,ddy,ddz,parallel->nCtProc);

    }

    //   Mass parameter for Nose-Hoover thermostat

    qmass=2.0*param->kinTemp0*X2(bath->tauT);

    rts2=1./X2(param->timeStep);

#ifdef _OPENMP
    #pragma omp parallel for default(none) shared(param,atom) private(i)
#endif
    for(i=parallel->fAtProc; i<parallel->lAtProc; i++)
    {
        vx[i]+=0.5*param->timeStep*fx[i]*rmass[i];
        vy[i]+=0.5*param->timeStep*fy[i]*rmass[i];
        vz[i]+=0.5*param->timeStep*fz[i]*rmass[i];
    }

    ener->kin=kinetic(parallel,vx,vy,vz,mass,dBuffer);

    lambdb=2.0*(ener->kin-param->kinTemp0)/qmass;
    lambdc=bath->chiT+param->timeStep*lambdb;
    lambda=0.5*(bath->chiT+lambdc);

    if(param->nConst>0)
        nosecycle=3;
    else
        nosecycle=4;

    for(k=0; k<nosecycle; k++)
    {

// move atoms by leapfrog algorithm
        l=0;
#ifdef _OPENMP
        #pragma omp parallel for default(none) shared(xt,yt,zt,xo,yo,zo,vxo,vyo,vzo,vxu,vyu,vzu,param,atom,lambda) private(i)
#endif
        for(i=parallel->fAtProc; i<parallel->lAtProc; i++)
        {

// update velocities

            vxu[l]=vxo[l]+param->timeStep*(fx[i]*rmass[i]-vx[i]*lambda);
            vyu[l]=vyo[l]+param->timeStep*(fy[i]*rmass[i]-vy[i]*lambda);
            vzu[l]=vzo[l]+param->timeStep*(fz[i]*rmass[i]-vz[i]*lambda);

// update positions

            x[i]=xo[l]+param->timeStep*vxu[l];
            y[i]=yo[l]+param->timeStep*vyu[l];
            z[i]=zo[l]+param->timeStep*vzu[l];

// Temporary storage of the uncorrected positions

            if(param->nConst>0)
            {
                xt[l]=x[i];
                yt[l]=y[i];
                zt[l]=z[i];
            }

            l++;
        }

        if( (param->nConst>0) && (k==0) )
        {

            if(parallel->nProc>1)
            {
                update_double_para(param,parallel,x,dBuffer);
                update_double_para(param,parallel,y,dBuffer);
                update_double_para(param,parallel,z,dBuffer);
            }

// Apply constraint with Shake algorithm.

            lf_shake(param,box,constList,parallel,x,y,z,ddx,ddy,ddz,rmass,nAtConst,strest,&virshakt,dBuffer);

            virshake+=virshakt;
            for(i=0; i<6; i++)
                stress[i]+=strest[i];

            l=0;
#ifdef _OPENMP
            #pragma omp parallel for default(none) shared(xt,yt,zt,vxu,vyu,vzu,param,atom,rts2) private(i)
#endif
            for(i=parallel->fAtProc; i<parallel->lAtProc; i++)
            {

// Corrected velocities

                vxu[l]+=(x[i]-xt[l])*param->rTimeStep;
                vyu[l]+=(y[i]-yt[l])*param->rTimeStep;
                vzu[l]+=(z[i]-zt[l])*param->rTimeStep;

// Corrected Forces

                fx[i]+=(x[i]-xt[l])*mass[i]*rts2;
                fy[i]+=(y[i]-yt[l])*mass[i]*rts2;
                fz[i]+=(z[i]-zt[l])*mass[i]*rts2;

                l++;
            }
        }

// calculate full timestep velocity

        l=0;
#ifdef _OPENMP
        #pragma omp parallel for default(none) shared(vxo,vyo,vzo,vxu,vyu,vzu,param,atom) private(i)
#endif
        for(i=parallel->fAtProc; i<parallel->lAtProc; i++)
        {

            vx[i]=0.5*(vxo[l]+vxu[l]);
            vy[i]=0.5*(vyo[l]+vyu[l]);
            vz[i]=0.5*(vzo[l]+vzu[l]);

            l++;
        }

// calculate kinetic energy

        ener->kin=kinetic(parallel,vx,vy,vz,mass,dBuffer);

        lambdb=2.0*(ener->kin-param->kinTemp0)/qmass;
        lambdc=bath->chiT+param->timeStep*lambdb;
        lambda=0.5*(bath->chiT+lambdc);

    }

    ener->virshake=virshake;

    stress_kinetic(parallel,vx,vy,vz,mass,stresk,dBuffer);

    box->stress1+=stress[0]+stresk[0];
    box->stress2+=stress[1]+stresk[1];
    box->stress3+=stress[2]+stresk[2];
    box->stress4+=stress[1]+stresk[1];
    box->stress5+=stress[3]+stresk[3];
    box->stress6+=stress[4]+stresk[4];
    box->stress7+=stress[2]+stresk[2];
    box->stress8+=stress[4]+stresk[4];
    box->stress9+=stress[5]+stresk[5];

    bath->chiT=lambdc;

    ener->conint+=param->timeStep*lambda*qmass/X2(bath->tauT);
    ener->consv=ener->conint+0.5*qmass*X2(lambda);

// periodic boundary condition

    image_update(parallel,box,x,y,z);

// updated velocity

    l=0;
#ifdef _OPENMP
    #pragma omp parallel for default(none) shared(vxu,vyu,vzu,param,atom) private(i)
#endif
    for(i=parallel->fAtProc; i<parallel->lAtProc; i++)
    {

        vx[i]=vxu[l];
        vy[i]=vyu[l];
        vz[i]=vzu[l];

        l++;

    }

    if(parallel->nProc>1)
    {
        update_double_para(param,parallel,x,dBuffer);
        update_double_para(param,parallel,y,dBuffer);
        update_double_para(param,parallel,z,dBuffer);

        update_double_para(param,parallel,vx,dBuffer);
        update_double_para(param,parallel,vy,dBuffer);
        update_double_para(param,parallel,vz,dBuffer);

        if(param->nConst>0)
        {
            update_double_para(param,parallel,fx,dBuffer);
            update_double_para(param,parallel,fy,dBuffer);
            update_double_para(param,parallel,fz,dBuffer);
        }

    }

}

void lf_npt_h(PARAM *param,ENERGY *ener,PBC *box,BATH *bath,CONSTRAINT constList[],PARALLEL *parallel,
              double *x,double *y,double *z,
              double *vx,double *vy,double *vz,
              double *fx,double *fy,double *fz,
              double *mass,double *rmass,int *nAtConst,double dBuffer[])
{
    int i,k,l,ia,ib,nosecycle;
    double lambda,lambdb,lambdc,rts2,qmass;
    double gamma,gammb,gammc,cbrga,pmass;
    double volume,masst=0.,cell0[9],com[3]= {0.},vom[3]= {0.};
    double virshake=0.,virshakt=0.,stress[6]= {0.},strest[6]= {0.},stresk[6]= {0.};

    l=0;
#ifdef _OPENMP
    #pragma omp parallel for default(none) shared(xo,yo,zo,vxo,vyo,vzo,param,atom) private(i)
#endif
    for(i=parallel->fAtProc; i<parallel->lAtProc; i++)
    {

// Store old coordinates and old velocities.

        xo[l]=x[i];
        yo[l]=y[i];
        zo[l]=z[i];

        vxo[l]=vx[i];
        vyo[l]=vy[i];
        vzo[l]=vz[i];

        l++;

    }

    //Store initial box parameters

    cell0[0]=box->a1;
    cell0[1]=box->a2;
    cell0[2]=box->a3;
    cell0[3]=box->b1;
    cell0[4]=box->b2;
    cell0[5]=box->b3;
    cell0[6]=box->c1;
    cell0[7]=box->c2;
    cell0[8]=box->c3;

    volume=box->vol;

    if(param->nConst>0)
    {

        l=0;
#ifdef _OPENMP
        #pragma omp parallel for default(none) shared(param,atom,dd,constList) private(i,ia,ib)
#endif
        for(i=parallel->fCtProc; i<parallel->lCtProc; i++)
        {
            ia=constList[i].a;
            ib=constList[i].b;

            ddx[l]=x[ib]-x[ia];
            ddy[l]=y[ib]-y[ia];
            ddz[l]=z[ib]-z[ia];

            l++;
        }

        image_array(box,ddx,ddy,ddz,parallel->nCtProc);

    }

    // total mass and center of mass

    masst=0.;
    com[0]=0.;
    com[1]=0.;
    com[2]=0.;
    for(i=0; i<param->nAtom; i++)
    {
        masst+=mass[i];
        com[0]+=mass[i]*x[i];
        com[1]+=mass[i]*y[i];
        com[2]+=mass[i]*z[i];
    }
    com[0]/=masst;
    com[1]/=masst;
    com[2]/=masst;

    //   Mass parameter for Nose-Hoover thermostat

    qmass=2.0*param->kinTemp0*X2(bath->tauT);
    pmass=2.0*param->kinTemp0*X2(bath->tauP);

    rts2=1./X2(param->timeStep);

#ifdef _OPENMP
    #pragma omp parallel for default(none) shared(param,atom) private(i)
#endif
    for(i=parallel->fAtProc; i<parallel->lAtProc; i++)
    {
        vx[i]+=0.5*param->timeStep*fx[i]*rmass[i];
        vy[i]+=0.5*param->timeStep*fy[i]*rmass[i];
        vz[i]+=0.5*param->timeStep*fz[i]*rmass[i];
    }

    /* *************************************************
     * This part has to be rework by improving kinetic()*/

    if(parallel->nProc>1)
    {
        update_double_para(param,parallel,vx,dBuffer);
        update_double_para(param,parallel,vy,dBuffer);
        update_double_para(param,parallel,vz,dBuffer);
    }

    ener->kin=kinetic(parallel,vx,vy,vz,mass,dBuffer);
    /* **************************************************/

    gammb=(2.0*ener->kin - ener->virpot - virshake - 3.0*param->press0*volume)/pmass-
          bath->chiT*bath->chiP;
    gammc=bath->chiP+param->timeStep*gammb;
    gamma=0.5*(bath->chiP+gammc);

    lambdb=(2.0*(ener->kin - param->kinTemp0 ) + pmass*X2(bath->chiP)-
            rboltzui*param->temp0)/qmass;
    lambdc=bath->chiT+param->timeStep*lambdb;
    lambda=0.5*(bath->chiT+lambdc);

    if(param->nConst>0)
        nosecycle=4;
    else
        nosecycle=5;

    for(k=0; k<nosecycle; k++)
    {

// move atoms by leapfrog algorithm

        l=0;
#ifdef _OPENMP
        #pragma omp parallel for default(none) shared(xt,yt,zt,xo,yo,zo,vxo,vyo,vzo,vxu,vyu,vzu,param,atom,lambda,gamma,com,gammc) private(i)
#endif
        for(i=parallel->fAtProc; i<parallel->lAtProc; i++)
        {

// update velocities

            vxu[l]=vxo[l]+param->timeStep*(fx[i]*rmass[i]-vx[i]*(lambda+gamma));
            vyu[l]=vyo[l]+param->timeStep*(fy[i]*rmass[i]-vy[i]*(lambda+gamma));
            vzu[l]=vzo[l]+param->timeStep*(fz[i]*rmass[i]-vz[i]*(lambda+gamma));

// update positions

            x[i]=xo[l]+param->timeStep*(vxu[l]+gammc*(0.5*(x[i]+xo[l])-com[0]));
            y[i]=yo[l]+param->timeStep*(vyu[l]+gammc*(0.5*(y[i]+yo[l])-com[1]));
            z[i]=zo[l]+param->timeStep*(vzu[l]+gammc*(0.5*(z[i]+zo[l])-com[2]));

// Temporary storage of the uncorrected positions

            if(param->nConst>0)
            {
                xt[l]=x[i];
                yt[l]=y[i];
                zt[l]=z[i];
            }

            l++;
        }

        if(param->nConst>0)
        {

            if(parallel->nProc>1)
            {
                update_double_para(param,parallel,x,dBuffer);
                update_double_para(param,parallel,y,dBuffer);
                update_double_para(param,parallel,z,dBuffer);
            }

// Apply constraint with Shake algorithm.

            cbrga=exp(3.*param->timeStep*gammc);
            cbrga=cbrt(cbrga);

            scale_box(box,cell0,cbrga);

            lf_shake(param,box,constList,parallel,x,y,z,ddx,ddy,ddz,rmass,nAtConst,strest,&virshakt,dBuffer);

            virshake+=virshakt;
            for(i=0; i<6; i++)
                stress[i]+=strest[i];

            l=0;
#ifdef _OPENMP
            #pragma omp parallel for default(none) shared(xt,yt,zt,vxu,vyu,vzu,param,atom,rts2) private(i)
#endif
            for(i=parallel->fAtProc; i<parallel->lAtProc; i++)
            {

// Corrected velocities

                vxu[l]+=(x[i]-xt[l])*param->rTimeStep;
                vyu[l]+=(y[i]-yt[l])*param->rTimeStep;
                vzu[l]+=(z[i]-zt[l])*param->rTimeStep;

// Corrected Forces

                fx[i]+=(x[i]-xt[l])*mass[i]*rts2;
                fy[i]+=(y[i]-yt[l])*mass[i]*rts2;
                fz[i]+=(z[i]-zt[l])*mass[i]*rts2;

                l++;
            }
        }

// calculate full timestep velocity

        l=0;
#ifdef _OPENMP
        #pragma omp parallel for default(none) shared(vxo,vyo,vzo,vxu,vyu,vzu,param,atom) private(i)
#endif
        for(i=parallel->fAtProc; i<parallel->lAtProc; i++)
        {

            vx[i]=0.5*(vxo[l]+vxu[l]);
            vy[i]=0.5*(vyo[l]+vyu[l]);
            vz[i]=0.5*(vzo[l]+vzu[l]);

            l++;
        }

// calculate kinetic energy

        /* *************************************************
        * This part has to be rework by improving kinetic()*/

        if(parallel->nProc>1)
        {
            update_double_para(param,parallel,vx,dBuffer);
            update_double_para(param,parallel,vy,dBuffer);
            update_double_para(param,parallel,vz,dBuffer);
        }

        ener->kin=kinetic(parallel,vx,vy,vz,mass,dBuffer);
        /* **************************************************/

        //ener->kin=kinetic(parallel,vx,vy,vz,mass,dBuffer);

        gammb=(2.0*ener->kin-ener->virpot-virshake-3.0*param->press0*volume)/pmass-
              bath->chiT*bath->chiP;
        gammc=bath->chiP+param->timeStep*gammb;
        gamma=0.5*(bath->chiP+gammc);

        lambdb=(2.0*(ener->kin-param->kinTemp0)+pmass*X2(bath->chiP)-
                rboltzui*param->temp0)/qmass;
        lambdc=bath->chiT+param->timeStep*lambdb;
        lambda=0.5*(bath->chiT+lambdc);

    }

    ener->virshake=virshake; // Has to be checked

    stress_kinetic(parallel,vx,vy,vz,mass,stresk,dBuffer);

    box->stress1+=stress[0]+stresk[0];
    box->stress2+=stress[1]+stresk[1];
    box->stress3+=stress[2]+stresk[2];
    box->stress4+=stress[1]+stresk[1];
    box->stress5+=stress[3]+stresk[3];
    box->stress6+=stress[4]+stresk[4];
    box->stress7+=stress[2]+stresk[2];
    box->stress8+=stress[4]+stresk[4];
    box->stress9+=stress[5]+stresk[5];

    cbrga=exp(3.*param->timeStep*gammc);
    cbrga=cbrt(cbrga);

    scale_box(box,cell0,cbrga);

    bath->chiT=lambdc;
    bath->chiP=gammc;

    ener->conint+=param->timeStep*lambda*(rboltzui*param->temp0+qmass/X2(bath->tauT));
    ener->consv=ener->conint+param->press0*volume+0.5*(qmass*X2(lambda)+pmass*X2(gamma));

// periodic boundary condition

    image_update(parallel,box,x,y,z);

// updated velocity
    l=0;
#ifdef _OPENMP
    #pragma omp parallel for default(none) shared(vxu,vyu,vzu,param,atom) private(i)
#endif
    for(i=parallel->fAtProc; i<parallel->lAtProc; i++)
    {

        vx[i]=vxu[l];
        vy[i]=vyu[l];
        vz[i]=vzu[l];

        l++;
    }

    vom[0]=0.;
    vom[1]=0.;
    vom[2]=0.;
    for(i=0; i<param->nAtom; i++)
    {
        vom[0]+=mass[i]*vx[i];
        vom[1]+=mass[i]*vy[i];
        vom[2]+=mass[i]*vz[i];
    }
    vom[0]/=masst;
    vom[1]/=masst;
    vom[2]/=masst;

    for(i=parallel->fAtProc; i<parallel->lAtProc; i++)
    {
        vx[i]-=vom[0];
        vy[i]-=vom[1];
        vz[i]-=vom[2];
    }

    if(parallel->nProc>1)
    {
        update_double_para(param,parallel,x,dBuffer);
        update_double_para(param,parallel,y,dBuffer);
        update_double_para(param,parallel,z,dBuffer);

        update_double_para(param,parallel,vx,dBuffer);
        update_double_para(param,parallel,vy,dBuffer);
        update_double_para(param,parallel,vz,dBuffer);

        if(param->nConst>0)
        {
            update_double_para(param,parallel,fx,dBuffer);
            update_double_para(param,parallel,fy,dBuffer);
            update_double_para(param,parallel,fz,dBuffer);
        }

    }

}

void vv_integrate(CTRL *ctrl,PARAM *param,ENERGY *ener,PBC *box,BATH *bath,
                  CONSTRAINT constList[],PARALLEL *parallel,double *x,double *y,
                  double *z,double *vx,double *vy,double *vz,
                  double *fx,double *fy,double *fz,
                  double *mass,double *rmass,int *nAtConst,double dBuffer[],int stage)
{

#ifdef TIMER
    update_timer_begin(TIMER_INTEGRATE,__func__);
#endif

    switch (ctrl->ens)
    {
    case NVE:
        vv_nve(param,ener,box,constList,parallel,x,y,z,vx,vy,vz,fx,fy,fz,mass,rmass,nAtConst,dBuffer,stage);
        break;
    case NVT_B:
        vv_nvt_b(param,ener,box,bath,constList,parallel,x,y,z,vx,vy,vz,fx,fy,fz,mass,rmass,nAtConst,dBuffer,stage);
        break;
    case NPT_B:
        vv_npt_b(param,ener,box,bath,constList,parallel,x,y,z,vx,vy,vz,fx,fy,fz,mass,rmass,nAtConst,dBuffer,stage);
        break;
    case NVT_H:
        vv_nvt_h(param,ener,box,bath,constList,parallel,x,y,z,vx,vy,vz,fx,fy,fz,mass,rmass,nAtConst,dBuffer,stage);
        break;
    case NPT_H:
        vv_npt_h(param,ener,box,bath,constList,parallel,x,y,z,vx,vy,vz,fx,fy,fz,mass,rmass,nAtConst,dBuffer,stage);
        break;
    default:
        vv_nve(param,ener,box,constList,parallel,x,y,z,vx,vy,vz,fx,fy,fz,mass,rmass,nAtConst,dBuffer,stage);
        break;
    }

#ifdef TIMER
    update_timer_end(TIMER_INTEGRATE,__func__);
#endif

}

void vv_nve(PARAM *param,ENERGY *ener,PBC *box,CONSTRAINT constList[],PARALLEL *parallel,
            double *x,double *y,double *z,double *vx,double *vy,double *vz,
            double *fx,double *fy,double *fz,
            double *mass,double *rmass,int *nAtConst,double dBuffer[],int stage)
{
    int i,l,ia,ib;
    double virshake,stress[6]= {0.},stresk[6]= {0.};

    if(param->nConst>0)
    {

        l=0;
#ifdef _OPENMP
        #pragma omp parallel for default(none) shared(constList,dd,param,atom) private(i,ia,ib)
#endif
        for(i=parallel->fCtProc; i<parallel->lCtProc; i++)
        {
            ia=constList[i].a;
            ib=constList[i].b;

            ddx[l]=x[ib]-x[ia];
            ddy[l]=y[ib]-y[ia];
            ddz[l]=z[ib]-z[ia];

            l++;
        }

        image_array(box,ddx,ddy,ddz,parallel->nCtProc);

    }

// move atoms by leapfrog algorithm

#ifdef _OPENMP
    #pragma omp parallel for default(none) shared(param,atom) private(i)
#endif
    for(i=parallel->fAtProc; i<parallel->lAtProc; i++)
    {
// update velocities

        vx[i]+=0.5*param->timeStep*fx[i]*rmass[i];
        vy[i]+=0.5*param->timeStep*fy[i]*rmass[i];
        vz[i]+=0.5*param->timeStep*fz[i]*rmass[i];
    }

    if(stage==1)
    {
#ifdef _OPENMP
        #pragma omp parallel for default(none) shared(param,atom) private(i)
#endif
        for(i=parallel->fAtProc; i<parallel->lAtProc; i++)
        {
// update positions

            x[i]+=param->timeStep*vx[i];
            y[i]+=param->timeStep*vy[i];
            z[i]+=param->timeStep*vz[i];

        }

        if(parallel->nProc>1)
        {
            update_double_para(param,parallel,x,dBuffer);
            update_double_para(param,parallel,y,dBuffer);
            update_double_para(param,parallel,z,dBuffer);
        }

        if(param->nConst>0)
        {

            if(parallel->nProc>1)
            {
                update_double_para(param,parallel,vx,dBuffer);
                update_double_para(param,parallel,vy,dBuffer);
                update_double_para(param,parallel,vz,dBuffer);
            }

// Apply constraint with Shake algorithm.

            vv_shake_r(param,box,constList,parallel,x,y,z,vx,vy,vz,ddx,ddy,ddz,rmass,nAtConst,stress,&virshake,dBuffer);
            ener->virshake=virshake;

        }

    }
    else
    {

        if(parallel->nProc>1)
        {
            update_double_para(param,parallel,vx,dBuffer);
            update_double_para(param,parallel,vy,dBuffer);
            update_double_para(param,parallel,vz,dBuffer);
        }

// calculate kinetic energy

        if(param->nConst>0)
        {

// Apply constraint with Shake algorithm.

            vv_shake_v(param,constList,parallel,vx,vy,vz,ddx,ddy,ddz,rmass,nAtConst,dBuffer);

        }

        ener->kin=kinetic(parallel,vx,vy,vz,mass,dBuffer);

        stress_kinetic(parallel,vx,vy,vz,mass,stresk,dBuffer);

        box->stress1+=stress[0]+stresk[0];
        box->stress2+=stress[1]+stresk[1];
        box->stress3+=stress[2]+stresk[2];
        box->stress4+=stress[1]+stresk[1];
        box->stress5+=stress[3]+stresk[3];
        box->stress6+=stress[4]+stresk[4];
        box->stress7+=stress[2]+stresk[2];
        box->stress8+=stress[4]+stresk[4];
        box->stress9+=stress[5]+stresk[5];

    }

    if(stage==2)
    {

// periodic boundary condition

        image_update(parallel,box,x,y,z);

        if(parallel->nProc>1)
        {
            update_double_para(param,parallel,x,dBuffer);
            update_double_para(param,parallel,y,dBuffer);
            update_double_para(param,parallel,z,dBuffer);
        }
    }

}

void vv_nvt_b(PARAM *param,ENERGY *ener,PBC *box,BATH *bath,CONSTRAINT constList[],PARALLEL *parallel,
              double *x,double *y,double *z,double *vx,double *vy,double *vz,
              double *fx,double *fy,double *fz,
              double *mass,double *rmass,int *nAtConst,double dBuffer[],int stage)
{
    int i,l,ia,ib;
    double lambda;
    double virshake,stress[6]= {0.},stresk[6]= {0.};

    l=0;
    if(param->nConst>0)
    {

#ifdef _OPENMP
        #pragma omp parallel for default(none) shared(constList,dd,param,atom) private(i,ia,ib)
#endif
        for(i=parallel->fCtProc; i<parallel->lCtProc; i++)
        {
            ia=constList[i].a;
            ib=constList[i].b;

            ddx[l]=x[ib]-x[ia];
            ddy[l]=y[ib]-y[ia];
            ddz[l]=z[ib]-z[ia];

            l++;
        }

        image_array(box,ddx,ddy,ddz,parallel->nCtProc);

    }

// move atoms by leapfrog algorithm

#ifdef _OPENMP
    #pragma omp parallel for default(none) shared(param,atom) private(i)
#endif
    for(i=parallel->fAtProc; i<parallel->lAtProc; i++)
    {
// update velocities

        vx[i]+=0.5*param->timeStep*fx[i]*rmass[i];
        vy[i]+=0.5*param->timeStep*fy[i]*rmass[i];
        vz[i]+=0.5*param->timeStep*fz[i]*rmass[i];
    }

    if(stage==1)
    {
#ifdef _OPENMP
        #pragma omp parallel for default(none) shared(param,atom) private(i)
#endif
        for(i=parallel->fAtProc; i<parallel->lAtProc; i++)
        {
// update positions

            x[i]+=param->timeStep*vx[i];
            y[i]+=param->timeStep*vy[i];
            z[i]+=param->timeStep*vz[i];

        }

        if(parallel->nProc>1)
        {
            update_double_para(param,parallel,x,dBuffer);
            update_double_para(param,parallel,y,dBuffer);
            update_double_para(param,parallel,z,dBuffer);
        }

        if(param->nConst>0)
        {

            if(parallel->nProc>1)
            {
                update_double_para(param,parallel,vx,dBuffer);
                update_double_para(param,parallel,vy,dBuffer);
                update_double_para(param,parallel,vz,dBuffer);
            }

// Apply constraint with Shake algorithm.

            vv_shake_r(param,box,constList,parallel,x,y,z,vx,vy,vz,ddx,ddy,ddz,rmass,nAtConst,stress,&virshake,dBuffer);
            ener->virshake=virshake;

        }

    }
    else
    {

        if(parallel->nProc>1)
        {
            update_double_para(param,parallel,vx,dBuffer);
            update_double_para(param,parallel,vy,dBuffer);
            update_double_para(param,parallel,vz,dBuffer);
        }

// calculate kinetic energy

        if(param->nConst>0)
        {

// Apply constraint with Shake algorithm.

            vv_shake_v(param,constList,parallel,vx,vy,vz,ddx,ddy,ddz,rmass,nAtConst,dBuffer);

        }

        ener->kin=kinetic(parallel,vx,vy,vz,mass,dBuffer);

        lambda=sqrt(1.0+param->timeStep/bath->tauT*(param->kinTemp0/ener->kin-1.0));

#ifdef _OPENMP
        #pragma omp parallel for default(none) shared(param,atom,lambda) private(i)
#endif
        for(i=parallel->fAtProc; i<parallel->lAtProc; i++)
        {
            vx[i]*=lambda;
            vy[i]*=lambda;
            vz[i]*=lambda;
        }

        ener->kin*=X2(lambda);

        if(parallel->nProc>1)
        {
            update_double_para(param,parallel,vx,dBuffer);
            update_double_para(param,parallel,vy,dBuffer);
            update_double_para(param,parallel,vz,dBuffer);
        }

        stress_kinetic(parallel,vx,vy,vz,mass,stresk,dBuffer);

        box->stress1+=stress[0]+stresk[0];
        box->stress2+=stress[1]+stresk[1];
        box->stress3+=stress[2]+stresk[2];
        box->stress4+=stress[1]+stresk[1];
        box->stress5+=stress[3]+stresk[3];
        box->stress6+=stress[4]+stresk[4];
        box->stress7+=stress[2]+stresk[2];
        box->stress8+=stress[4]+stresk[4];
        box->stress9+=stress[5]+stresk[5];

    }

    if(stage==2)
    {

// periodic boundary condition

        image_update(parallel,box,x,y,z);

        if(parallel->nProc>1)
        {
            update_double_para(param,parallel,x,dBuffer);
            update_double_para(param,parallel,y,dBuffer);
            update_double_para(param,parallel,z,dBuffer);
        }
    }

}

void vv_npt_b(PARAM *param,ENERGY *ener,PBC *box,BATH *bath,CONSTRAINT constList[],PARALLEL *parallel,
              double *x,double *y,double *z,double *vx,double *vy,double *vz,
              double *fx,double *fy,double *fz,
              double *mass,double *rmass,int *nAtConst,double dBuffer[],int stage)
{
    int i,l,ia,ib,k,nosecycle;
    double lambda,gamma,cbrga,volume,pp;
    double virshake,stress[6]= {0.},stresk[6]= {0.};

    volume=box->vol;

    if(param->nConst>0)
    {
        l=0;
#ifdef _OPENMP
        #pragma omp parallel for default(none) shared(constList,dd,param,atom) private(i,ia,ib)
#endif
        for(i=parallel->fCtProc; i<parallel->lCtProc; i++)
        {
            ia=constList[i].a;
            ib=constList[i].b;

            ddx[l]=x[ib]-x[ia];
            ddy[l]=y[ib]-y[ia];
            ddz[l]=z[ib]-z[ia];

            l++;
        }

        image_array(box,ddx,ddy,ddz,parallel->nCtProc);

    }

    if(stage==1)
    {

        ener->kin=kinetic(parallel,vx,vy,vz,mass,dBuffer);

#ifdef _OPENMP
        #pragma omp parallel for default(none) shared(param,atom) private(i)
#endif
        for(i=parallel->fAtProc; i<parallel->lAtProc; i++)
        {
//    update velocities

            vx[i]+=0.5*param->timeStep*fx[i]*rmass[i];
            vy[i]+=0.5*param->timeStep*fy[i]*rmass[i];
            vz[i]+=0.5*param->timeStep*fz[i]*rmass[i];
        }

        if(param->nConst>0)
        {
            l=0;
#ifdef _OPENMP
            #pragma omp parallel for default(none) shared(xo,yo,zo,vxo,vyo,vzo,param,atom) private(i)
#endif
            for(i=parallel->fAtProc; i<parallel->lAtProc; i++)
            {

                // Store old coordinates and old velocities.

                xo[l]=x[i];
                yo[l]=y[i];
                zo[l]=z[i];

                vxo[l]=vx[i];
                vyo[l]=vy[i];
                vzo[l]=vz[i];

                l++;

            }
        }

        if( (stage==1) && (param->nConst>0) )
            nosecycle=2;
        else
            nosecycle=1;

        for(k=0; k<nosecycle; k++)
        {
            cbrga=1.;

            if(k==nosecycle-1)
            {
                pp=(2.*ener->kin-ener->virpot-virshake)/(3.*volume);
                gamma=1.+watercomp*param->timeStep*(pp-param->press0)/bath->tauP;
                cbrga=cbrt(gamma);

                vv_scale_box(box,cbrga);

            }

#ifdef _OPENMP
            #pragma omp parallel for default(none) shared(param,atom,cbrga) private(i)
#endif
            for(i=parallel->fAtProc; i<parallel->lAtProc; i++)
            {
                // update positions

                x[i]=cbrga*x[i]+param->timeStep*vx[i];
                y[i]=cbrga*y[i]+param->timeStep*vy[i];
                z[i]=cbrga*z[i]+param->timeStep*vz[i];

            }

            if(parallel->nProc>1)
            {
                update_double_para(param,parallel,x,dBuffer);
                update_double_para(param,parallel,y,dBuffer);
                update_double_para(param,parallel,z,dBuffer);
            }

            if(param->nConst>0)
            {

                if(parallel->nProc>1)
                {
                    update_double_para(param,parallel,vx,dBuffer);
                    update_double_para(param,parallel,vy,dBuffer);
                    update_double_para(param,parallel,vz,dBuffer);
                }

                // Apply constraint with Shake algorithm.

                vv_shake_r(param,box,constList,parallel,x,y,z,vx,vy,vz,ddx,ddy,ddz,rmass,nAtConst,stress,&virshake,dBuffer);
                ener->virshake=virshake;

            }

            if(k<nosecycle-1)
            {
                l=0;
#ifdef _OPENMP
                #pragma omp parallel for default(none) shared(xo,yo,zo,vxo,vyo,vzo,param,atom) private(i)
#endif
                for(i=parallel->fAtProc; i<parallel->lAtProc; i++)
                {

                    // Store old coordinates and old velocities.

                    x[i]=xo[l];
                    y[i]=yo[l];
                    z[i]=zo[l];

                    vx[i]=vxo[l];
                    vy[i]=vyo[l];
                    vz[i]=vzo[l];

                    l++;

                }
            }
        }
    }
    else
    {

#ifdef _OPENMP
        #pragma omp parallel for default(none) shared(param,atom) private(i)
#endif
        for(i=parallel->fAtProc; i<parallel->lAtProc; i++)
        {
//    update velocities

            vx[i]+=0.5*param->timeStep*fx[i]*rmass[i];
            vy[i]+=0.5*param->timeStep*fy[i]*rmass[i];
            vz[i]+=0.5*param->timeStep*fz[i]*rmass[i];
        }

// calculate kinetic energy

        ener->kin=kinetic(parallel,vx,vy,vz,mass,dBuffer);

        lambda=sqrt(1.0+param->timeStep/bath->tauT*(param->kinTemp0/ener->kin-1.0));

        for(i=parallel->fAtProc; i<parallel->lAtProc; i++)
        {
//    update velocities

            vx[i]*=lambda;
            vy[i]*=lambda;
            vz[i]*=lambda;
        }

        if(parallel->nProc>1)
        {
            update_double_para(param,parallel,vx,dBuffer);
            update_double_para(param,parallel,vy,dBuffer);
            update_double_para(param,parallel,vz,dBuffer);
        }

        if(param->nConst>0)
        {

// Apply constraint with Shake algorithm.

            vv_shake_v(param,constList,parallel,vx,vy,vz,ddx,ddy,ddz,rmass,nAtConst,dBuffer);

        }

        ener->kin=kinetic(parallel,vx,vy,vz,mass,dBuffer);

        stress_kinetic(parallel,vx,vy,vz,mass,stresk,dBuffer);

        box->stress1+=stress[0]+stresk[0];
        box->stress2+=stress[1]+stresk[1];
        box->stress3+=stress[2]+stresk[2];
        box->stress4+=stress[1]+stresk[1];
        box->stress5+=stress[3]+stresk[3];
        box->stress6+=stress[4]+stresk[4];
        box->stress7+=stress[2]+stresk[2];
        box->stress8+=stress[4]+stresk[4];
        box->stress9+=stress[5]+stresk[5];

    }

    if(stage==2)
    {

// periodic boundary condition

        image_update(parallel,box,x,y,z);

        if(parallel->nProc>1)
        {
            update_double_para(param,parallel,x,dBuffer);
            update_double_para(param,parallel,y,dBuffer);
            update_double_para(param,parallel,z,dBuffer);
        }
    }

}

void vv_nvt_h(PARAM *param,ENERGY *ener,PBC *box,BATH *bath,CONSTRAINT constList[],PARALLEL *parallel,
              double *x,double *y,double *z,double *vx,double *vy,double *vz,
              double *fx,double *fy,double *fz,
              double *mass,double *rmass,int *nAtConst,double dBuffer[],int stage)
{
    int i,l,ia,ib;
    double lambda,qmass;
    double virshake,stress[6]= {0.},stresk[6]= {0.};

    l=0;
    if(param->nConst>0)
    {

#ifdef _OPENMP
        #pragma omp parallel for default(none) shared(constList,dd,param,atom) private(i,ia,ib)
#endif
        for(i=parallel->fCtProc; i<parallel->lCtProc; i++)
        {
            ia=constList[i].a;
            ib=constList[i].b;

            ddx[l]=x[ib]-x[ia];
            ddy[l]=y[ib]-y[ia];
            ddz[l]=z[ib]-z[ia];

            l++;
        }

        image_array(box,ddx,ddy,ddz,parallel->nCtProc);

    }

    qmass=2.0*param->kinTemp0*X2(bath->tauT);

    if(stage==1)
    {

        ener->kin=kinetic(parallel,vx,vy,vz,mass,dBuffer);

        bath->chiT+=0.5*param->timeStep*(ener->kin-param->kinTemp0)/qmass;

        lambda=exp(-0.5*param->timeStep*bath->chiT);

#ifdef _OPENMP
        #pragma omp parallel for default(none) shared(param,atom,lambda) private(i)
#endif
        for(i=parallel->fAtProc; i<parallel->lAtProc; i++)
        {
            // scale velocities

            vx[i]*=lambda;
            vy[i]*=lambda;
            vz[i]*=lambda;

            // update velocities

            vx[i]+=0.5*param->timeStep*fx[i]*rmass[i];
            vy[i]+=0.5*param->timeStep*fy[i]*rmass[i];
            vz[i]+=0.5*param->timeStep*fz[i]*rmass[i];
        }

        ener->kin*=X2(lambda);

        ener->conint+=0.5*param->timeStep*bath->chiT*qmass/X2(bath->tauT);

        bath->chiT+=0.5*param->timeStep*(ener->kin-param->kinTemp0)/qmass;

#ifdef _OPENMP
        #pragma omp parallel for default(none) shared(param,atom) private(i)
#endif
        for(i=parallel->fAtProc; i<parallel->lAtProc; i++)
        {
// update positions

            x[i]+=param->timeStep*vx[i];
            y[i]+=param->timeStep*vy[i];
            z[i]+=param->timeStep*vz[i];

        }

        if(parallel->nProc>1)
        {
            update_double_para(param,parallel,x,dBuffer);
            update_double_para(param,parallel,y,dBuffer);
            update_double_para(param,parallel,z,dBuffer);
        }

        if(param->nConst>0)
        {

            if(parallel->nProc>1)
            {
                update_double_para(param,parallel,vx,dBuffer);
                update_double_para(param,parallel,vy,dBuffer);
                update_double_para(param,parallel,vz,dBuffer);
            }

// Apply constraint with Shake algorithm.

            vv_shake_r(param,box,constList,parallel,x,y,z,vx,vy,vz,ddx,ddy,ddz,rmass,nAtConst,stress,&virshake,dBuffer);
            ener->virshake=virshake;

        }

    }
    else
    {
// calculate kinetic energy

#ifdef _OPENMP
        #pragma omp parallel for default(none) shared(param,atom) private(i)
#endif
        for(i=parallel->fAtProc; i<parallel->lAtProc; i++)
        {
            // update velocities

            vx[i]+=0.5*param->timeStep*fx[i]*rmass[i];
            vy[i]+=0.5*param->timeStep*fy[i]*rmass[i];
            vz[i]+=0.5*param->timeStep*fz[i]*rmass[i];
        }

        if(param->nConst>0)
        {

            if(parallel->nProc>1)
            {
                update_double_para(param,parallel,vx,dBuffer);
                update_double_para(param,parallel,vy,dBuffer);
                update_double_para(param,parallel,vz,dBuffer);
            }

// Apply constraint with Shake algorithm.

            vv_shake_v(param,constList,parallel,vx,vy,vz,ddx,ddy,ddz,rmass,nAtConst,dBuffer);

        }

        ener->kin=kinetic(parallel,vx,vy,vz,mass,dBuffer);

        bath->chiT+=0.5*param->timeStep*(ener->kin-param->kinTemp0)/qmass;

        lambda=exp(-0.5*param->timeStep*bath->chiT);

#ifdef _OPENMP
        #pragma omp parallel for default(none) shared(param,atom,lambda) private(i)
#endif
        for(i=parallel->fAtProc; i<parallel->lAtProc; i++)
        {
            vx[i]*=lambda;
            vy[i]*=lambda;
            vz[i]*=lambda;
        }

        ener->kin*=X2(lambda);

        ener->conint+=0.5*param->timeStep*bath->chiT*qmass/X2(bath->tauT);

        bath->chiT+=0.5*param->timeStep*(ener->kin-param->kinTemp0)/qmass;

        ener->consv=ener->conint+0.5*qmass*X2(bath->chiT);

        if(parallel->nProc>1)
        {
            update_double_para(param,parallel,vx,dBuffer);
            update_double_para(param,parallel,vy,dBuffer);
            update_double_para(param,parallel,vz,dBuffer);
        }

        stress_kinetic(parallel,vx,vy,vz,mass,stresk,dBuffer);

        box->stress1+=stress[0]+stresk[0];
        box->stress2+=stress[1]+stresk[1];
        box->stress3+=stress[2]+stresk[2];
        box->stress4+=stress[1]+stresk[1];
        box->stress5+=stress[3]+stresk[3];
        box->stress6+=stress[4]+stresk[4];
        box->stress7+=stress[2]+stresk[2];
        box->stress8+=stress[4]+stresk[4];
        box->stress9+=stress[5]+stresk[5];

    }

    if(stage==2)
    {

// periodic boundary condition

        image_update(parallel,box,x,y,z);

        if(parallel->nProc>1)
        {
            update_double_para(param,parallel,x,dBuffer);
            update_double_para(param,parallel,y,dBuffer);
            update_double_para(param,parallel,z,dBuffer);
        }
    }

}

void vv_npt_h(PARAM *param,ENERGY *ener,PBC *box,BATH *bath,CONSTRAINT constList[],PARALLEL *parallel,
              double *x,double *y,double *z,double *vx,double *vy,double *vz,
              double *fx,double *fy,double *fz,
              double *mass,double *rmass,int *nAtConst,double dBuffer[],int stage)
{
    int i,l,ia,ib,k,kk,nosecycle,hoovercycle=5;
    double hts,chts,cqts;
    double cons0,lambda,lambda0,qmass;
    double gamma,gamma0,pmass,cbrga,scale;
    double volume,volume0,cell0[9],masst=0.,com[3]= {0.},vom[3]= {0.};
    double virshake,stress[6]= {0.},stresk[6]= {0.};

    hts=0.5*param->timeStep;
    chts=hts/(double)hoovercycle;
    cqts=0.25*param->timeStep/(double)hoovercycle;

    //Store initial box parameters

    cell0[0]=box->a1;
    cell0[1]=box->a2;
    cell0[2]=box->a3;
    cell0[3]=box->b1;
    cell0[4]=box->b2;
    cell0[5]=box->b3;
    cell0[6]=box->c1;
    cell0[7]=box->c2;
    cell0[8]=box->c3;

    volume=box->vol;
    volume0=box->vol;

    masst=0.;
    for(i=0; i<param->nAtom; i++)
        masst+=mass[i];

    qmass=2.0*param->kinTemp0*X2(bath->tauT);
    pmass=2.0*param->kinTemp0*X2(bath->tauP);

    if(param->nConst>0)
    {
        l=0;
#ifdef _OPENMP
        #pragma omp parallel for default(none) shared(param,constList,dd,atom) private(i,ia,ib)
#endif
        for(i=0; i<param->nConst; i++)
        {
            ia=constList[i].a;
            ib=constList[i].b;

            ddx[l]=x[ib]-x[ia];
            ddy[l]=y[ib]-y[ia];
            ddz[l]=z[ib]-z[ia];

            l++;
        }

        image_array(box,ddx,ddy,ddz,parallel->nCtProc);

    }

    if(stage==1)
    {

        lambda0=bath->chiT;
        gamma0=bath->chiP;
        cons0=ener->conint;

        l=0;
#ifdef _OPENMP
        #pragma omp parallel for default(none) shared(param,xo,yo,zo,vxo,vyo,vzo,atom) private(i)
#endif
        for(i=parallel->fAtProc; i<parallel->lAtProc; i++)
        {

            // Store old coordinates and old velocities.

            xo[l]=x[i];
            yo[l]=y[i];
            zo[l]=z[i];

            vxo[l]=vx[i];
            vyo[l]=vy[i];
            vzo[l]=vz[i];

            l++;
        }

        if( (stage==1) && (param->nConst>0) )
            nosecycle=2;
        else
            nosecycle=1;

        for(k=0; k<nosecycle; k++)
        {

            for(kk=0; kk<hoovercycle; kk++)
            {

                // apply nvt
                ener->kin=kinetic(parallel,vx,vy,vz,mass,dBuffer);

                bath->chiT+=0.5*cqts*(2.0*(ener->kin-param->kinTemp0)+
                                      pmass*X2(bath->chiP)-rboltzui*param->temp0)/qmass;

                lambda=exp(-cqts*bath->chiT);

#ifdef _OPENMP
                #pragma omp parallel for default(none) shared(param,atom,lambda) private(i)
#endif
                for(i=parallel->fAtProc; i<parallel->lAtProc; i++)
                {
                    // scale velocities

                    vx[i]*=lambda;
                    vy[i]*=lambda;
                    vz[i]*=lambda;

                }

                ener->kin*=X2(lambda);

                ener->conint+=cqts*bath->chiT*(rboltzui*param->temp0+qmass/X2(bath->tauT));

                bath->chiT+=0.5*cqts*(2.0*(ener->kin-param->kinTemp0)+
                                      pmass*X2(bath->chiP)-rboltzui*param->temp0)/qmass;

                // apply npt

                bath->chiP+=0.5*chts*(((2.0*ener->kin-ener->virpot-virshake)-
                                       3.0*param->press0*volume)/pmass-bath->chiP*bath->chiT);

                gamma=exp(-chts*bath->chiP);

                for(i=parallel->fAtProc; i<parallel->lAtProc; i++)
                {
                    // scale velocities

                    vx[i]*=gamma;
                    vy[i]*=gamma;
                    vz[i]*=gamma;

                }

                ener->kin*=X2(gamma);

                volume*=exp(3.0*chts*bath->chiP);

                bath->chiP+=0.5*chts*(((2.0*ener->kin-ener->virpot-virshake)-
                                       3.0*param->press0*volume)/pmass-bath->chiP*bath->chiT);

                // apply nvt

                /**************** check **************/
                ener->kin=kinetic(parallel,vx,vy,vz,mass,dBuffer);

                bath->chiT+=0.5*cqts*(2.0*(ener->kin-param->kinTemp0)+
                                      pmass*X2(bath->chiP)-rboltzui*param->temp0)/qmass;

                lambda=exp(-cqts*bath->chiT);

#ifdef _OPENMP
                #pragma omp parallel for default(none) shared(param,atom,lambda) private(i)
#endif
                for(i=parallel->fAtProc; i<parallel->lAtProc; i++)
                {
                    // scale velocities

                    vx[i]*=lambda;
                    vy[i]*=lambda;
                    vz[i]*=lambda;

                }

                ener->kin*=X2(lambda);

                ener->conint+=cqts*bath->chiT*(rboltzui*param->temp0+qmass/X2(bath->tauT));

                bath->chiT+=0.5*cqts*(2.0*(ener->kin-param->kinTemp0)+
                                      pmass*X2(bath->chiP)-rboltzui*param->temp0)/qmass;
            }

            scale=cbrt(volume/volume0);
            scale_box(box,cell0,scale);

#ifdef _OPENMP
            #pragma omp parallel for default(none) shared(param,atom,lambda) private(i)
#endif
            for(i=parallel->fAtProc; i<parallel->lAtProc; i++)
            {
                // update velocities

                vx[i]+=0.5*param->timeStep*fx[i]*rmass[i];
                vy[i]+=0.5*param->timeStep*fy[i]*rmass[i];
                vz[i]+=0.5*param->timeStep*fz[i]*rmass[i];
            }

            com[0]=0.;
            com[1]=0.;
            com[2]=0.;
            for(i=0; i<param->nAtom; i++)
            {
                com[0]+=mass[i]*x[i];
                com[1]+=mass[i]*y[i];
                com[2]+=mass[i]*z[i];
            }
            com[0]/=masst;
            com[1]/=masst;
            com[2]/=masst;

            cbrga=exp(param->timeStep*bath->chiP);

#ifdef _OPENMP
            #pragma omp parallel for default(none) shared(param,atom,com,cbrga) private(i)
#endif
            for(i=parallel->fAtProc; i<parallel->lAtProc; i++)
            {
                // update positions

                x[i]=cbrga*(x[i]-com[0])+param->timeStep*vx[i]+com[0];
                y[i]=cbrga*(y[i]-com[1])+param->timeStep*vy[i]+com[1];
                z[i]=cbrga*(z[i]-com[2])+param->timeStep*vz[i]+com[2];

            }

            if(parallel->nProc>1)
            {
                update_double_para(param,parallel,x,dBuffer);
                update_double_para(param,parallel,y,dBuffer);
                update_double_para(param,parallel,z,dBuffer);
            }

            if(param->nConst>0)
            {

                if(parallel->nProc>1)
                {
                    update_double_para(param,parallel,vx,dBuffer);
                    update_double_para(param,parallel,vy,dBuffer);
                    update_double_para(param,parallel,vz,dBuffer);
                }

                // Apply constraint with Shake algorithm.

                vv_shake_r(param,box,constList,parallel,x,y,z,vx,vy,vz,ddx,ddy,ddz,rmass,nAtConst,stress,&virshake,dBuffer);
                ener->virshake=virshake;

            }

            if(k<nosecycle-1)
            {
                volume=volume0;
                bath->chiT=lambda0;
                bath->chiP=gamma0;
                ener->conint=cons0;

                l=0;
#ifdef _OPENMP
                #pragma omp parallel for default(none) shared(xo,yo,zo,vxo,vyo,vzo,atom,param) private(i)
#endif
                for(i=parallel->fAtProc; i<parallel->lAtProc; i++)
                {

                    // Store old coordinates and old velocities.

                    xo[l]=x[i];
                    yo[l]=y[i];
                    zo[l]=z[i];

                    vxo[l]=vx[i];
                    vyo[l]=vy[i];
                    vzo[l]=vz[i];

                    l++;

                }

            }
        }
    }
    else
    {
// calculate kinetic energy

#ifdef _OPENMP
        #pragma omp parallel for default(none) shared(param,atom) private(i)
#endif
        for(i=parallel->fAtProc; i<parallel->lAtProc; i++)
        {
            // update velocities

            vx[i]+=0.5*param->timeStep*fx[i]*rmass[i];
            vy[i]+=0.5*param->timeStep*fy[i]*rmass[i];
            vz[i]+=0.5*param->timeStep*fz[i]*rmass[i];
        }

        if(param->nConst>0)
        {

            if(parallel->nProc>1)
            {
                update_double_para(param,parallel,vx,dBuffer);
                update_double_para(param,parallel,vy,dBuffer);
                update_double_para(param,parallel,vz,dBuffer);
            }

// Apply constraint with Shake algorithm.

            vv_shake_v(param,constList,parallel,vx,vy,vz,ddx,ddy,ddz,rmass,nAtConst,dBuffer);

        }

        for(kk=0; kk<hoovercycle; kk++)
        {

            // apply nvt
            ener->kin=kinetic(parallel,vx,vy,vz,mass,dBuffer);

            bath->chiT+=0.5*cqts*(2.0*(ener->kin-param->kinTemp0)+
                                  pmass*X2(bath->chiP)-rboltzui*param->temp0)/qmass;

            lambda=exp(-cqts*bath->chiT);

#ifdef _OPENMP
            #pragma omp parallel for default(none) shared(param,atom,lambda) private(i)
#endif
            for(i=parallel->fAtProc; i<parallel->lAtProc; i++)
            {
                // scale velocities

                vx[i]*=lambda;
                vy[i]*=lambda;
                vz[i]*=lambda;

            }

            ener->kin*=X2(lambda);

            ener->conint+=cqts*bath->chiT*(rboltzui*param->temp0+qmass/X2(bath->tauT));

            bath->chiT+=0.5*cqts*(2.0*(ener->kin-param->kinTemp0)+
                                  pmass*X2(bath->chiP)-rboltzui*param->temp0)/qmass;

            // apply npt

            bath->chiP+=0.5*chts*(((2.0*ener->kin-ener->virpot-virshake)-
                                   3.0*param->press0*volume)/pmass-bath->chiP*bath->chiT);

            gamma=exp(-chts*bath->chiP);

            for(i=parallel->fAtProc; i<parallel->lAtProc; i++)
            {
                // scale velocities

                vx[i]*=gamma;
                vy[i]*=gamma;
                vz[i]*=gamma;

            }

            ener->kin*=X2(gamma);

            volume*=exp(3.0*chts*bath->chiP);

            bath->chiP+=0.5*chts*(((2.0*ener->kin-ener->virpot-virshake)-
                                   3.0*param->press0*volume)/pmass-bath->chiP*bath->chiT);

            // apply nvt

            /*************** check *******************/
            ener->kin=kinetic(parallel,vx,vy,vz,mass,dBuffer);

            bath->chiT+=0.5*cqts*(2.0*(ener->kin-param->kinTemp0)+
                                  pmass*X2(bath->chiP)-rboltzui*param->temp0)/qmass;

            lambda=exp(-cqts*bath->chiT);

#ifdef _OPENMP
            #pragma omp parallel for default(none) shared(param,atom,lambda) private(i)
#endif
            for(i=parallel->fAtProc; i<parallel->lAtProc; i++)
            {
                // scale velocities

                vx[i]*=lambda;
                vy[i]*=lambda;
                vz[i]*=lambda;

            }

            ener->kin*=X2(lambda);

            ener->conint+=cqts*bath->chiT*(rboltzui*param->temp0+qmass/X2(bath->tauT));

            bath->chiT+=0.5*cqts*(2.0*(ener->kin-param->kinTemp0)+
                                  pmass*X2(bath->chiP)-rboltzui*param->temp0)/qmass;
        }

        vom[0]=0.;
        vom[1]=0.;
        vom[2]=0.;
        for(i=0; i<param->nAtom; i++)
        {
            vom[0]+=mass[i]*vx[i];
            vom[1]+=mass[i]*vy[i];
            vom[2]+=mass[i]*vz[i];
        }
        vom[0]/=masst;
        vom[1]/=masst;
        vom[2]/=masst;

        for(i=parallel->fAtProc; i<parallel->lAtProc; i++)
        {
            vx[i]-=vom[0];
            vy[i]-=vom[1];
            vz[i]-=vom[2];
        }

        if(parallel->nProc>1)
        {
            update_double_para(param,parallel,vx,dBuffer);
            update_double_para(param,parallel,vy,dBuffer);
            update_double_para(param,parallel,vz,dBuffer);
        }

        scale=cbrt(volume/volume0);
        scale_box(box,cell0,scale);

        ener->consv=ener->conint+param->press0*volume+0.5*(qmass*X2(bath->chiT)+pmass*X2(bath->chiP));

        ener->kin=kinetic(parallel,vx,vy,vz,mass,dBuffer);

        stress_kinetic(parallel,vx,vy,vz,mass,stresk,dBuffer);

        box->stress1+=stress[0]+stresk[0];
        box->stress2+=stress[1]+stresk[1];
        box->stress3+=stress[2]+stresk[2];
        box->stress4+=stress[1]+stresk[1];
        box->stress5+=stress[3]+stresk[3];
        box->stress6+=stress[4]+stresk[4];
        box->stress7+=stress[2]+stresk[2];
        box->stress8+=stress[4]+stresk[4];
        box->stress9+=stress[5]+stresk[5];

    }

    if(stage==2)
    {

// periodic boundary condition

        image_update(parallel,box,x,y,z);

        if(parallel->nProc>1)
        {
            update_double_para(param,parallel,x,dBuffer);
            update_double_para(param,parallel,y,dBuffer);
            update_double_para(param,parallel,z,dBuffer);
        }
    }

}
