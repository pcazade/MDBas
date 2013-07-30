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
 * \file energy.c
 * \brief Contains highest level functions for evaluating energy and forces of a system.
 * \author Pierre-Andre Cazade and Florent Hedin
 * \version alpha-branch
 * \date 2012
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "global.h"
#include "energy.h"
#include "elec.h"
#include "vdw.h"
#include "ewald.h"
#include "spme.h"
#include "energy.h"
#include "internal.h"
#include "io.h"
#include "utils.h"
#include "errors.h"
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

/**
 * \param simulCond Pointer to structure SIMULPARAMS containing parameters of the current simulation.
 *
 * \brief Set the pointers to function used for non-bonded energy evaluation.
 */
void init_energy_ptrs(CTRL *ctrl)
{

    switch(ctrl->elecType)
    {
    case NOELEC:
        ptr_coulomb = &(coulomb_none);
        ptr_coulomb14 = &(coulomb14_none);
        break;

    case FULL:
        ptr_coulomb = &(coulomb_none);
        ptr_coulomb14 = &(coulomb14_none);
        break;

    case SHIFT1:
        ptr_coulomb = &(coulomb_shift1);
        ptr_coulomb14 = &(coulomb14_shift1);
        break;

    case SHIFT2:
        ptr_coulomb = &(coulomb_shift2);
        ptr_coulomb14 = &(coulomb14_shift2);
        break;

    case SWITCH:
        ptr_coulomb = &(coulomb_switch);
        ptr_coulomb14 = &(coulomb14_switch);
        break;

    default:
        my_error(UNKNOWN_ELEC_ERROR,__FILE__,__LINE__,0);
        break;
    }

    switch(ctrl->vdwType)
    {
    case NOVDW:
        ptr_vdw = &(vdw_none);
        ptr_vdw14 = &(vdw14_none);
        break;

    case VFULL:
        ptr_vdw = &(vdw_none);
        ptr_vdw14 = &(vdw14_none);
        break;

    case VSWITCH:
        ptr_vdw = &(vdw_switch);
        ptr_vdw14 = &(vdw14_switch);
        break;

    default:
        my_error(UNKNOWN_VDW_ERROR,__FILE__,__LINE__,0);
        break;
    }

}

/**
 * \param atom Array of structure ATOM (coordinates, forces, etc...).
 * \param ff Pointer to structure FORCEFIELD containing forcefield parameters.
 * \param ener Pointer to structure ENERGY containing values of the different energies.
 * \param simulCond Pointer to structure SIMULPARAMS containing parameters of the current simulation.
 * \param box Pointer to structure PBC containing Periodic Boundaries Conditions parameters.
 *
 * \brief Main energy function, collecting total energy by calling the required subfunctions.
 */
void energy(CTRL *ctrl,PARAM *param,PARALLEL *parallel,ENERGY *ener,EWALD *ewald,PBC *box,NEIGH *neigh,
            BOND bond[],BOND ub[],ANGLE angle[],DIHE dihe[],DIHE impr[],
            const double x[],const double y[], const double z[],
            double vx[],double vy[], double vz[],double fx[],double fy[],
            double fz[],const double q[],const double eps[],const double sig[],
            const double eps14[],const double sig14[],const int frozen[],
            int **neighList,const int neighPair[],const int neighList14[],
            int **exclList,const int exclPair[],double dBuffer[])
{

    int i;

    ener->elec=0.;
    ener->vdw=0.;
    ener->bond=0.;
    ener->ub=0.;
    ener->ang=0.;
    ener->dihe=0.;
    ener->impr=0.;

    ener->virelec=0.;
    ener->virvdw=0.;
    ener->virbond=0.;
    ener->virub=0.;
    ener->virpot=0.;

    box->stress1=0.;
    box->stress2=0.;
    box->stress3=0.;
    box->stress4=0.;
    box->stress5=0.;
    box->stress6=0.;
    box->stress7=0.;
    box->stress8=0.;
    box->stress9=0.;

    for(i=0; i<param->nAtom; i++)
    {
        fx[i]=0.;
        fy[i]=0.;
        fz[i]=0.;
    }

#ifdef TIMER
    update_timer_begin(TIMER_ENERGY_TOT,__func__);
#endif

    /* Performing non-bonding terms */

    if(ctrl->keyEwald!=0)
    {
        ewald_energy(ctrl,param,parallel,ener,ewald,box,x,y,z,fx,fy,fz,q,eps,sig,
                     neighList,neighPair,exclList,exclPair,dBuffer);
        if(ctrl->keyNb14)
            ewald14_energy(param,parallel,ener,ewald,box,neigh,x,y,z,fx,fy,fz,q,eps14,sig14,neighList14);
    }
    else
    {
        nonbond_energy(param,parallel,ener,box,x,y,z,fx,fy,fz,q,eps,sig,neighList,neighPair);
        if(ctrl->keyNb14)
            nonbond14_energy(param,parallel,ener,box,neigh,x,y,z,fx,fy,fz,q,eps14,sig14,neighList14);
    }

    /* Performing bond terms */

    if(param->nBond>0)
        bond_energy(parallel,ener,box,bond,x,y,z,fx,fy,fz);

    /* Performing angle terms */

    if(param->nAngle>0)
        angle_energy(parallel,ener,box,angle,x,y,z,fx,fy,fz);

    /* Performing Urey-Bradley terms */

    if(param->nUb>0)
        ub_energy(parallel,ener,box,ub,x,y,z,fx,fy,fz);

    /* Performing diherdral terms */

    if(param->nDihedral>0)
        dihedral_energy(parallel,ener,box,dihe,x,y,z,fx,fy,fz);

    /* Performing improper terms */

    if(param->nImproper>0)
        improper_energy(parallel,ener,box,impr,x,y,z,fx,fy,fz);

    /* Calculate potential energy */

    if(parallel->nProc>1)
    {
        dBuffer[0]=ener->elec;
        dBuffer[1]=ener->vdw;
        dBuffer[2]=ener->bond;
        dBuffer[3]=ener->ang;
        dBuffer[4]=ener->ub;
        dBuffer[5]=ener->dihe;
        dBuffer[6]=ener->impr;

        dBuffer[7]=ener->virelec;
        dBuffer[8]=ener->virvdw;
        dBuffer[9]=ener->virbond;
        dBuffer[10]=ener->virub;

        sum_double_para(dBuffer,&(dBuffer[11]),11);

        ener->elec=dBuffer[0];
        ener->vdw=dBuffer[1];
        ener->bond=dBuffer[2];
        ener->ang=dBuffer[3];
        ener->ub=dBuffer[4];
        ener->dihe=dBuffer[5];
        ener->impr=dBuffer[6];

        ener->virelec=dBuffer[7];
        ener->virvdw=dBuffer[8];
        ener->virbond=dBuffer[9];
        ener->virub=dBuffer[10];

        sum_double_para(fx,dBuffer,param->nAtom);
        sum_double_para(fy,dBuffer,param->nAtom);
        sum_double_para(fz,dBuffer,param->nAtom);
    }

    ener->pot=ener->elec+ener->vdw+ener->bond+
              ener->ang+ener->ub+ener->dihe+ener->impr;

    ener->virpot=ener->virbond+ener->virub+ener->virelec+ener->virvdw;

#ifdef TIMER
    update_timer_end(TIMER_ENERGY_TOT,__func__);
#endif

    if(param->nFrozen>0)
        freeze_atoms(param,vx,vy,vz,fx,fy,fz,frozen);

}

/**
 * \param atom Array of structure ATOM (coordinates, forces, etc...).
 * \param ff Pointer to structure FORCEFIELD containing forcefield parameters.
 * \param ener Pointer to structure ENERGY containing values of the different energies.
 * \param simulCond Pointer to structure SIMULPARAMS containing parameters of the current simulation.
 * \param box Pointer to structure PBC containing Periodic Boundaries Conditions parameters.
 *
 * \brief Energy function collecting all terms of non-bonded energies.
 */
void nonbond_energy(PARAM *param,PARALLEL *parallel,ENERGY *ener,PBC *box,const double x[],const double y[],
                    const double z[],double fx[],double fy[],double fz[],const double q[],
                    const double eps[],const double sig[],int **neighList,
                    const int neighPair[])
{

    int i,j,k,l;
    double elec=0.,evdw=0.,delec=0.,dvdw=0.,virelec=0.,virvdw=0.;
    double r,r2,rt,fxj,fyj,fzj,fxi,fyi,fzi;
    double qel,veps,vsig;
    double delta[3]/*,stress[6]={0.}*/;

#ifdef TIMER
    update_timer_begin(TIMER_ENERGY_NB,__func__);
#endif

    l=0;
#ifdef _OPENMP
    #pragma omp parallel default(none) shared(atom,param,box,vdw,neigh,ptr_coulomb,ptr_vdw) private(i,j,k,delec,dvdw,r,fx,fy,fz,fxi,fyi,fzi,delta) reduction(+:elec,evdw,virelec,virvdw)
    {
        #pragma omp for schedule(dynamic,1000) nowait
#endif
        for(i=parallel->idProc; i<param->nAtom; i+=parallel->nProc)
        {

            fxi=0.;
            fyi=0.;
            fzi=0.;

            for(k=0; k<neighPair[l]; k++)
            {
                j=neighList[l][k];

                delta[0]=x[j]-x[i];
                delta[1]=y[j]-y[i];
                delta[2]=z[j]-z[i];

                r2=dist(box,delta);
                r=sqrt(r2);
                rt=1./r;

                if(r2<=param->cutOff2)
                {
                    qel=q[i]*q[j];
                    elec+=(*ptr_coulomb)(param,&delec,qel,r2,rt);

                    veps=eps[i]*eps[j];
                    vsig=sig[i]+sig[j];
                    evdw+=(*ptr_vdw)(param,&dvdw,veps,vsig,r2,rt);

                    virelec+=delec*r;
                    virvdw+=dvdw*r;

                    fxj=(delec+dvdw)*delta[0]*rt;
                    fyj=(delec+dvdw)*delta[1]*rt;
                    fzj=(delec+dvdw)*delta[2]*rt;

                    fxi+=fxj;
                    fyi+=fyj;
                    fzi+=fzj;

                    /*stress[0]-=fx*delta[0];
                    stress[1]-=fy*delta[0];
                    stress[2]-=fz*delta[0];
                    stress[3]-=fy*delta[1];
                    stress[4]-=fz*delta[1];
                    stress[5]-=fz*delta[2];*/

#ifdef _OPENMP
                    #pragma omp atomic
#endif
                    fx[j]+=-fxj;

#ifdef _OPENMP
                    #pragma omp atomic
#endif
                    fy[j]+=-fyj;

#ifdef _OPENMP
                    #pragma omp atomic
#endif
                    fz[j]+=-fzj;

                }

            }

#ifdef _OPENMP
            #pragma omp atomic
#endif
            fx[i]+=fxi;

#ifdef _OPENMP
            #pragma omp atomic
#endif
            fy[i]+=fyi;

#ifdef _OPENMP
            #pragma omp atomic
#endif
            fz[i]+=fzi;

            l++;

        } //end of parallel for

#ifdef _OPENMP
    } //end of parallel region
#endif

    ener->elec+=elec;
    ener->vdw+=evdw;

    ener->virelec+=virelec;
    ener->virvdw+=virvdw;

    /*box->stress1+=stress[0];
    box->stress2+=stress[1];
    box->stress3+=stress[2];
    box->stress4+=stress[1];
    box->stress5+=stress[3];
    box->stress6+=stress[4];
    box->stress7+=stress[2];
    box->stress8+=stress[4];
    box->stress9+=stress[5];*/

#ifdef TIMER
    update_timer_end(TIMER_ENERGY_NB,__func__);
#endif

}

/**
 * \param atom Array of structure ATOM (coordinates, forces, etc...).
 * \param ff Pointer to structure FORCEFIELD containing forcefield parameters.
 * \param ener Pointer to structure ENERGY containing values of the different energies.
 * \param simulCond Pointer to structure SIMULPARAMS containing parameters of the current simulation.
 * \param box Pointer to structure PBC containing Periodic Boundaries Conditions parameters.
 *
 * \brief Energy function collecting all terms of 1-4 non-bonded energies.
 */
void nonbond14_energy(PARAM *param,PARALLEL *parallel,ENERGY *ener,PBC *box,NEIGH *neigh,
                      const double x[],const double y[],const double z[],double fx[],double fy[],
                      double fz[],const double q[],const double eps[],const double sig[],
                      const int neighList14[])
{

    int i,j,k,l;
    double elec=0.,evdw=0.,delec=0.,dvdw=0.,virelec=0.,virvdw=0.;
    double r,r2,rt,fxj,fyj,fzj;
    double qel,veps,vsig;
    double delta[3]/*,stress[6]={0.}*/;

#ifdef TIMER
    update_timer_begin(TIMER_ENERGY_NB14,__func__);
#endif

    l=0;
    for(k=parallel->idProc; k<neigh->nPair14; k+=parallel->nProc)
    {
        i=neighList14[2*l];
        j=neighList14[2*l+1];

        delec=0.;
        dvdw=0.;

        delta[0]=x[j]-x[i];
        delta[1]=y[j]-y[i];
        delta[2]=z[j]-z[i];

        r2=dist(box,delta);
        r=sqrt(r2);
        rt=1./r;

        qel=q[i]*q[j];
        elec+=(*ptr_coulomb14)(param,&delec,qel,r2,rt);

        veps=eps[i]*eps[j];
        vsig=sig[i]+sig[j];
        evdw+=(*ptr_vdw14)(param,&dvdw,veps,vsig,r2,rt);

        virelec+=delec*r;
        virvdw+=dvdw*r;

        fxj=(delec+dvdw)*delta[0]*rt;
        fyj=(delec+dvdw)*delta[1]*rt;
        fzj=(delec+dvdw)*delta[2]*rt;

        /*stress[0]-=fx*delta[0];
        stress[1]-=fy*delta[0];
        stress[2]-=fz*delta[0];
        stress[3]-=fy*delta[1];
        stress[4]-=fz*delta[1];
        stress[5]-=fz*delta[2];*/

        fx[i]+=fxj;
        fy[i]+=fyj;
        fz[i]+=fzj;

        fx[j]+=-fxj;
        fy[j]+=-fyj;
        fz[j]+=-fzj;

        l++;
    }

    ener->elec+=elec;
    ener->vdw+=evdw;

    ener->virelec+=virelec;
    ener->virvdw+=virvdw;

    /*box->stress1+=stress[0];
    box->stress2+=stress[1];
    box->stress3+=stress[2];
    box->stress4+=stress[1];
    box->stress5+=stress[3];
    box->stress6+=stress[4];
    box->stress7+=stress[2];
    box->stress8+=stress[4];
    box->stress9+=stress[5];*/

#ifdef TIMER
    update_timer_end(TIMER_ENERGY_NB14,__func__);
#endif

}

/**
 * \param atom Array of structure ATOM (coordinates, forces, etc...).
 * \param ff Pointer to structure FORCEFIELD containing forcefield parameters.
 * \param ener Pointer to structure ENERGY containing values of the different energies.
 * \param simulCond Pointer to structure SIMULPARAMS containing parameters of the current simulation.
 * \param box Pointer to structure PBC containing Periodic Boundaries Conditions parameters.
 *
 * \brief Energy function collecting all terms of non-bonded energies.
 */
void ewald_energy(CTRL *ctrl,PARAM *param,PARALLEL *parallel,ENERGY *ener,EWALD *ewald,PBC *box,const double x[],
                  const double y[],const double z[],double fx[],double fy[],
                  double fz[],const double q[],const double eps[],const double sig[],
                  int **neighList,const int neighPair[],int *exclList[],
                  const int exclPair[],double dBuffer[])
{

    int i,j,k,l;
    double eEwaldRec=0.,virEwaldRec=0.,eEwaldDir=0.,dEwaldDir=0.,virEwaldDir=0.;
    double eEwaldCorr=0.,dEwaldCorr=0.,virEwaldCorr=0.;
    double evdw=0.,dvdw=0.,virvdw=0.;
    double r,r2,rt,fxj,fyj,fzj,fxi,fyi,fzi;
    double qel,veps,vsig;
    double delta[3],stress1[6]= {0.};

#ifdef TIMER
    update_timer_begin(TIMER_ENERGY_NB,__func__);
#endif

    if(ctrl->keyEwald==1)
        eEwaldRec=ewald_rec(param,parallel,ewald,box,x,y,z,fx,fy,fz,q,stress1,&virEwaldRec,dBuffer);
    else if(ctrl->keyEwald==2)
        eEwaldRec=spme_energy(param,parallel,ewald,box,x,y,z,fx,fy,fz,q,stress1,&virEwaldRec,dBuffer);
    
    l=0;
#ifdef _OPENMP
    #pragma omp parallel default(none) shared(atom,param,box,vdw,neigh,ptr_coulomb,ptr_vdw) private(i,j,k,delec,dvdw,r,fx,fy,fz,fxi,fyi,fzi,delta) reduction(+:elec,evdw,virelec,virvdw)
    {
        #pragma omp for schedule(dynamic,1000) nowait
#endif
        for(i=parallel->idProc; i<param->nAtom; i+=parallel->nProc)
        {

            fxi=0.;
            fyi=0.;
            fzi=0.;

            for(k=0; k<neighPair[l]; k++)
            {
                j=neighList[l][k];

                delta[0]=x[j]-x[i];
                delta[1]=y[j]-y[i];
                delta[2]=z[j]-z[i];

                r2=dist(box,delta);
                r=sqrt(r2);
                rt=1./r;

                if(r2<=param->cutOff2)
                {
                    qel=param->chargeConst*q[i]*q[j];
                    eEwaldDir+=ewald_dir(ewald,&dEwaldDir,qel,r,rt);

                    veps=eps[i]*eps[j];
                    vsig=sig[i]+sig[j];
                    evdw+=(*ptr_vdw)(param,&dvdw,veps,vsig,r2,rt);

                    virEwaldDir+=dEwaldDir*r;
                    virvdw+=dvdw*r;

                    fxj=(dEwaldDir+dvdw)*delta[0]*rt;
                    fyj=(dEwaldDir+dvdw)*delta[1]*rt;
                    fzj=(dEwaldDir+dvdw)*delta[2]*rt;

                    fxi+=fxj;
                    fyi+=fyj;
                    fzi+=fzj;

                    /*stress[0]-=fx*delta[0];
                    stress[1]-=fy*delta[0];
                    stress[2]-=fz*delta[0];
                    stress[3]-=fy*delta[1];
                    stress[4]-=fz*delta[1];
                    stress[5]-=fz*delta[2];*/

#ifdef _OPENMP
                    #pragma omp atomic
#endif
                    fx[j]+=-fxj;

#ifdef _OPENMP
                    #pragma omp atomic
#endif
                    fy[j]+=-fyj;

#ifdef _OPENMP
                    #pragma omp atomic
#endif
                    fz[j]+=-fzj;

                }
            }

#ifdef _OPENMP
            #pragma omp atomic
#endif
            fx[i]+=fxi;

#ifdef _OPENMP
            #pragma omp atomic
#endif
            fy[i]+=fyi;

#ifdef _OPENMP
            #pragma omp atomic
#endif
            fz[i]+=fzi;

            l++;

        } //end of parallel for

#ifdef _OPENMP
    } //end of parallel region
#endif

    l=0;
    for(i=parallel->idProc; i<param->nAtom; i+=parallel->nProc)
    {

        fxi=0.;
        fyi=0.;
        fzi=0.;

        for(k=0; k<exclPair[l]; k++)
        {
            j=exclList[l][k];

            if(j>i)
            {

                delta[0]=x[j]-x[i];
                delta[1]=y[j]-y[i];
                delta[2]=z[j]-z[i];

                r2=dist(box,delta);
                r=sqrt(r2);
                rt=1./r;

                qel=param->chargeConst*q[i]*q[j];
                eEwaldCorr+=ewald_corr(ewald,&dEwaldCorr,qel,r,rt);

                virEwaldCorr+=dEwaldCorr*r;

                fxj=dEwaldCorr*delta[0]*rt;
                fyj=dEwaldCorr*delta[1]*rt;
                fzj=dEwaldCorr*delta[2]*rt;

                fxi+=fxj;
                fyi+=fyj;
                fzi+=fzj;

                /*stress[0]-=fx*delta[0];
                stress[1]-=fy*delta[0];
                stress[2]-=fz*delta[0];
                stress[3]-=fy*delta[1];
                stress[4]-=fz*delta[1];
                stress[5]-=fz*delta[2];*/

                fx[j]+=-fxj;
                fy[j]+=-fyj;
                fz[j]+=-fzj;

            }

        }

        fx[i]+=fxi;
        fy[i]+=fyi;
        fz[i]+=fzi;

        l++;
    }

    ener->elec+=eEwaldDir+eEwaldRec+eEwaldCorr;
    ener->vdw+=evdw;

    ener->virelec+=virEwaldDir+virEwaldRec+virEwaldCorr;
    ener->virvdw+=virvdw;

    /*box->stress1+=stress[0];
    box->stress2+=stress[1];
    box->stress3+=stress[2];
    box->stress4+=stress[1];
    box->stress5+=stress[3];
    box->stress6+=stress[4];
    box->stress7+=stress[2];
    box->stress8+=stress[4];
    box->stress9+=stress[5];*/

#ifdef TIMER
    update_timer_end(TIMER_ENERGY_NB,__func__);
#endif

}

/**
 * \param atom Array of structure ATOM (coordinates, forces, etc...).
 * \param ff Pointer to structure FORCEFIELD containing forcefield parameters.
 * \param ener Pointer to structure ENERGY containing values of the different energies.
 * \param simulCond Pointer to structure SIMULPARAMS containing parameters of the current simulation.
 * \param box Pointer to structure PBC containing Periodic Boundaries Conditions parameters.
 *
 * \brief Energy function collecting all terms of 1-4 non-bonded energies.
 */
void ewald14_energy(PARAM *param,PARALLEL *parallel,ENERGY *ener,EWALD *ewald,PBC *box,NEIGH *neigh,const double x[],
                    const double y[],const double z[],double fx[],double fy[],
                    double fz[],const double q[],const double eps[],const double sig[],
                    const int neighList14[])
{

    int i,j,k,l;
    double eEwaldDir=0.,dEwaldDir=0.,virEwaldDir=0.;
    double eEwaldCorr=0.,dEwaldCorr=0.,virEwaldCorr=0.;
    double evdw=0.,dvdw=0.,virvdw=0.;
    double r,r2,rt,fxj,fyj,fzj;
    double qel,veps,vsig;
    double delta[3]/*,stress[6]={0.}*/;

#ifdef TIMER
    update_timer_begin(TIMER_ENERGY_NB14,__func__);
#endif

    l=0;
    for(k=parallel->idProc; k<neigh->nPair14; k+=parallel->nProc)
    {
        i=neighList14[2*l];
        j=neighList14[2*l+1];

        dEwaldCorr=0.;
        dvdw=0.;

        delta[0]=x[j]-x[i];
        delta[1]=y[j]-y[i];
        delta[2]=z[j]-z[i];

        r2=dist(box,delta);
        r=sqrt(r2);
        rt=1./r;

        qel=param->chargeConst*q[i]*q[j];

        eEwaldDir+=ewald_dir14(param,ewald,&dEwaldDir,qel,r,rt);
        eEwaldCorr+=ewald_corr14(param,ewald,&dEwaldCorr,qel,r,rt);

        veps=eps[i]*eps[j];
        vsig=sig[i]+sig[j];
        evdw+=(*ptr_vdw14)(param,&dvdw,veps,vsig,r2,rt);

        virEwaldDir+=dEwaldDir*r;
        virEwaldCorr+=dEwaldCorr*r;
        virvdw+=dvdw*r;

        fxj=(dEwaldDir+dEwaldCorr+dvdw)*delta[0]*rt;
        fyj=(dEwaldDir+dEwaldCorr+dvdw)*delta[1]*rt;
        fzj=(dEwaldDir+dEwaldCorr+dvdw)*delta[2]*rt;

        /*stress[0]-=fx*delta[0];
        stress[1]-=fy*delta[0];
        stress[2]-=fz*delta[0];
        stress[3]-=fy*delta[1];
        stress[4]-=fz*delta[1];
        stress[5]-=fz*delta[2];*/

        fx[i]+=fxj;
        fy[i]+=fyj;
        fz[i]+=fzj;

        fx[j]+=-fxj;
        fy[j]+=-fyj;
        fz[j]+=-fzj;

        l++;

    }

    ener->elec+=eEwaldDir+eEwaldCorr;
    ener->vdw+=evdw;

    ener->virelec+=virEwaldDir+virEwaldCorr;
    ener->virvdw+=virvdw;

    /*box->stress1+=stress[0];
    box->stress2+=stress[1];
    box->stress3+=stress[2];
    box->stress4+=stress[1];
    box->stress5+=stress[3];
    box->stress6+=stress[4];
    box->stress7+=stress[2];
    box->stress8+=stress[4];
    box->stress9+=stress[5];*/

#ifdef TIMER
    update_timer_end(TIMER_ENERGY_NB14,__func__);
#endif

}
