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
 * \file internal.c
 * \brief Contains functions evaluating the internal energy.
 * \author Pierre-Andre Cazade and Florent Hedin
 * \version alpha-branch
 * \date 2012
 */

#include <stdio.h>
#include <float.h>
#include <math.h>

#include "global.h"
#include "utils.h"

#if (defined TIMING && defined __unix__ && !defined __STRICT_ANSI__)
#define TIMER
#include "timing.h"
#endif

void bond_energy(const PARAM *param,const PARALLEL *parallel,ENERGY *ener,const PBC *box,
                 const BOND bond[],const double *x,const double *y,const double *z,
                 double *fx,double *fy,double *fz)
{
    int i,j,ll;
    double morsea,morseb;
    double r,tfx,tfy,tfz,dbond,virbond=0.;
    double delta[3]/*,stress[6]={0.}*/;

#ifdef TIMER
    update_timer_begin(TIMER_ENERGY_BOND,__func__);
#endif

    for(ll=parallel->fBdProc; ll<parallel->lBdProc; ll++)
    {

        i=bond[ll].a;
        j=bond[ll].b;

        delta[0]=x[j]-x[i];
        delta[1]=y[j]-y[i];
        delta[2]=z[j]-z[i];

        r=sqrt(dist(box,delta));

        switch(bond[ll].type)
        {
        case BHARM:
            ener->bond+=0.5*bond[ll].k*X2(r-bond[ll].r0);
            dbond=bond[ll].k*(r-bond[ll].r0);
            break;

        case BMORSE:
            morsea=exp(-bond[ll].beta*(r-bond[ll].r0));
            morseb=X2(morsea);
            ener->bond+=bond[ll].k*(morseb-2.*morsea)+bond[ll].k;
            dbond=2.*bond[ll].k*bond[ll].beta*(morsea-morseb);
            break;

        default:
            ener->bond+=0.5*bond[ll].k*X2(r-bond[ll].r0);
            dbond=bond[ll].k*(r-bond[ll].r0);
            break;

        }

        virbond+=dbond*r;

        tfx=dbond*delta[0]/r;
        tfy=dbond*delta[1]/r;
        tfz=dbond*delta[2]/r;

        fx[i]+=tfx;
        fy[i]+=tfy;
        fz[i]+=tfz;

        fx[j]+=-tfx;
        fy[j]+=-tfy;
        fz[j]+=-tfz;

        /*stress[0]-=tfx*delta[0];
        stress[1]-=tfy*delta[0];
        stress[2]-=tfz*delta[0];
        stress[3]-=tfy*delta[1];
        stress[4]-=tfz*delta[1];
        stress[5]-=tfz*delta[2];*/

    }

    ener->virbond+=virbond;

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
    update_timer_end(TIMER_ENERGY_BOND,__func__);
#endif

}

void ub_energy(const PARAM *param,const PARALLEL *parallel,ENERGY *ener,const PBC *box,
               const BOND ub[],const double *x,const double *y,const double *z,
               double *fx,double *fy,double *fz)
{
    int i,j,ll;
    double r,tfx,tfy,tfz,dub,virub=0.;
    double delta[3]/*,stress[6]={0.}*/;

#ifdef TIMER
    update_timer_begin(TIMER_ENERGY_UB,__func__);
#endif

    for(ll=parallel->fUbProc; ll<parallel->lUbProc; ll++)
    {
        i=ub[ll].a;
        j=ub[ll].b;

        delta[0]=x[j]-x[i];
        delta[1]=y[j]-y[i];
        delta[2]=z[j]-z[i];

        r=sqrt(dist(box,delta));

        ener->ub+=0.5*ub[ll].k*X2(r-ub[ll].r0);
        dub=ub[ll].k*(r-ub[ll].r0);

        virub+=dub*r;

        tfx=dub*delta[0]/r;
        tfy=dub*delta[1]/r;
        tfz=dub*delta[2]/r;

        fx[i]+=tfx;
        fy[i]+=tfy;
        fz[i]+=tfz;

        fx[j]+=-tfx;
        fy[j]+=-tfy;
        fz[j]+=-tfz;

        /*stress[0]-=tfx*delta[0];
        stress[1]-=tfy*delta[0];
        stress[2]-=tfz*delta[0];
        stress[3]-=tfy*delta[1];
        stress[4]-=tfz*delta[1];
        stress[5]-=tfz*delta[2];*/

    }

    ener->virub+=virub;

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
    update_timer_end(TIMER_ENERGY_UB,__func__);
#endif

}

void angle_energy(const PARAM *param,const PARALLEL *parallel,ENERGY *ener,const PBC *box,
                  const ANGLE angle[],const double *x,const double *y,const double *z,
                  double *fx,double *fy,double *fz)
{
    int i,j,k,ll;
    double dangle,rab,rbc,rabt,rbct,cost,sint,theta;
    double dab[3],dbc[3]/*,stress[6]*/;
    double fxa,fya,fza,fxc,fyc,fzc;

#ifdef TIMER
    update_timer_begin(TIMER_ENERGY_ANGL,__func__);
#endif

    for(ll=parallel->fAgProc; ll<parallel->lAgProc; ll++)
    {

        i=angle[ll].a;
        j=angle[ll].b;
        k=angle[ll].c;

        dab[0]=x[i]-x[j];
        dab[1]=y[i]-y[j];
        dab[2]=z[i]-z[j];

        rab=sqrt(dist(box,dab));
        rabt=1./rab;

        dbc[0]=x[k]-x[j];
        dbc[1]=y[k]-y[j];
        dbc[2]=z[k]-z[j];

        rbc=sqrt(dist(box,dbc));
        rbct=1./rbc;

        cost=(dab[0]*dbc[0]+dab[1]*dbc[1]+dab[2]*dbc[2])/(rab*rbc);
        sint=MAX(1.e-8,sqrt(1.-(cost*cost)));
        theta=acos(cost);

        ener->ang+=0.5*angle[ll].k*X2(theta-angle[ll].theta0);
        dangle=angle[ll].k*(theta-angle[ll].theta0)/sint;

        fxa=dangle*(dbc[0]*rbct-dab[0]*cost*rabt)*rabt;
        fya=dangle*(dbc[1]*rbct-dab[1]*cost*rabt)*rabt;
        fza=dangle*(dbc[2]*rbct-dab[2]*cost*rabt)*rabt;

        fxc=dangle*(dab[0]*rabt-dbc[0]*cost*rbct)*rbct;
        fyc=dangle*(dab[1]*rabt-dbc[1]*cost*rbct)*rbct;
        fzc=dangle*(dab[2]*rabt-dbc[2]*cost*rbct)*rbct;

        fx[i]+=fxa;
        fy[i]+=fya;
        fz[i]+=fza;

        fx[j]+=-fxa-fxc;
        fy[j]+=-fya-fyc;
        fz[j]+=-fza-fzc;

        fx[k]+=fxc;
        fy[k]+=fyc;
        fz[k]+=fzc;

        /*stress[0]=rab*fxa*dab[0]+rbc*fxc*dbc[0];
        stress[1]=rab*fya*dab[0]+rbc*fyc*dbc[0];
        stress[2]=rab*fza*dab[0]+rbc*fzc*dbc[0];
        stress[3]=rab*fya*dab[1]+rbc*fyc*dbc[1];
        stress[4]=rab*fza*dab[1]+rbc*fzc*dbc[1];
        stress[5]=rab*fza*dab[2]+rbc*fzc*dbc[2];*/

    }

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
    update_timer_end(TIMER_ENERGY_ANGL,__func__);
#endif

}

void dihedral_energy(const PARAM *param,const PARALLEL *parallel,ENERGY *ener,const PBC *box,
                     const DIHE dihe[],const double *x,const double *y,const double *z,
                     double *fx,double *fy,double *fz)
{
    int i,j,k,l,ll;
    double edihe=0.,ddihe=0.;
    double cosp,sinp,phi;
    double /*rab,*/rbc,/*rcd,*/rpb,rpc,r2pb,r2pc,pbpc;
    double dab[3],dbc[3],dcd[3],/*dac[3],*/pb[3],pc[3]/*,stress[6]={0.}*/;
    double fax,fay,faz,fbx,fby,fbz,fcx,fcy,fcz,fdx,fdy,fdz;

#ifdef TIMER
    update_timer_begin(TIMER_ENERGY_DIHE,__func__);
#endif

    for(ll=parallel->fDhProc; ll<parallel->lDhProc; ll++)
    {
        i=dihe[ll].a;
        j=dihe[ll].b;
        k=dihe[ll].c;
        l=dihe[ll].d;

        dab[0]=x[i]-x[j];
        dab[1]=y[i]-y[j];
        dab[2]=z[i]-z[j];

        dist(box,dab);

        dbc[0]=x[j]-x[k];
        dbc[1]=y[j]-y[k];
        dbc[2]=z[j]-z[k];

        rbc=sqrt(dist(box,dbc));

        dcd[0]=x[k]-x[l];
        dcd[1]=y[k]-y[l];
        dcd[2]=z[k]-z[l];

        dist(box,dcd);

// construct first dihedral vector

        pb[0]=dab[1]*dbc[2]-dab[2]*dbc[1];
        pb[1]=dab[2]*dbc[0]-dab[0]*dbc[2];
        pb[2]=dab[0]*dbc[1]-dab[1]*dbc[0];

        r2pb=X2(pb[0])+X2(pb[1])+X2(pb[2]);
        rpb=sqrt(r2pb);

// construct second dihedral vector

        pc[0]=dbc[1]*dcd[2]-dbc[2]*dcd[1];
        pc[1]=dbc[2]*dcd[0]-dbc[0]*dcd[2];
        pc[2]=dbc[0]*dcd[1]-dbc[1]*dcd[0];

        r2pc=X2(pc[0])+X2(pc[1])+X2(pc[2]);
        rpc=sqrt(r2pc);

// determine dihedral angle

        pbpc=pb[0]*pc[0]+pb[1]*pc[1]+pb[2]*pc[2];
        cosp=pbpc/(rpb*rpc);

        sinp=(dbc[0]*(pc[1]*pb[2]-pc[2]*pb[1])+dbc[1]*(pb[0]*pc[2]-pb[2]*pc[0])+
              dbc[2]*(pc[0]*pb[1]-pc[1]*pb[0]))/(rpb*rpc*rbc);

        phi=atan2(sinp,cosp);

// avoid singularity in sinp

        if(sinp>=0.)
        {
            sinp=MAX(DBL_EPSILON,fabs(sinp));
        }
        else
        {
            sinp=-(MAX(DBL_EPSILON,fabs(sinp)));
        }

// calculate potential energy and scalar force term

        switch(dihe[ll].type)
        {
        case DCOS: // cosine dihedral

            edihe=dihe[ll].k*(1.+cos(dihe[ll].mult*phi-dihe[ll].phi0));
            ddihe=-dihe[ll].k*dihe[ll].mult*
                  sin(dihe[ll].mult*phi-dihe[ll].phi0)/(rpb*rpc*sinp);

            break;

        case DHARM: // harmonic dihedral

            phi=phi-dihe[ll].phi0;
            phi=phi-nint(phi/TWOPI)*TWOPI;

            edihe=0.5*dihe[ll].k*(phi*phi);

            ddihe=dihe[ll].k*phi/(rpb*rpc*sinp);

            break;

        default:

            edihe=dihe[ll].k*(1.+cos(dihe[ll].mult*phi-dihe[ll].phi0));
            ddihe=-dihe[ll].k*dihe[ll].mult*
                  sin(dihe[ll].mult*phi-dihe[ll].phi0)/(rpb*rpc*sinp);

            break;
        }

// calculate potential energy

        ener->dihe+=edihe;

        fax=ddihe*((-pc[1]*dbc[2]+pc[2]*dbc[1])-pbpc*(-pb[1]*dbc[2]+pb[2]*dbc[1])/r2pb);
        fay=ddihe*(( pc[0]*dbc[2]-pc[2]*dbc[0])-pbpc*( pb[0]*dbc[2]-pb[2]*dbc[0])/r2pb);
        faz=ddihe*((-pc[0]*dbc[1]+pc[1]*dbc[0])-pbpc*(-pb[0]*dbc[1]+pb[1]*dbc[0])/r2pb);

        fcx=ddihe*((-pc[1]*dab[2]+pc[2]*dab[1])-pbpc*(-pb[1]*dab[2]+pb[2]*dab[1])/r2pb);
        fcy=ddihe*(( pc[0]*dab[2]-pc[2]*dab[0])-pbpc*( pb[0]*dab[2]-pb[2]*dab[0])/r2pb);
        fcz=ddihe*((-pc[0]*dab[1]+pc[1]*dab[0])-pbpc*(-pb[0]*dab[1]+pb[1]*dab[0])/r2pb);

        fbx=ddihe*((-pb[1]*dcd[2]+pb[2]*dcd[1])-pbpc*(-pc[1]*dcd[2]+pc[2]*dcd[1])/r2pc);
        fby=ddihe*(( pb[0]*dcd[2]-pb[2]*dcd[0])-pbpc*( pc[0]*dcd[2]-pc[2]*dcd[0])/r2pc);
        fbz=ddihe*((-pb[0]*dcd[1]+pb[1]*dcd[0])-pbpc*(-pc[0]*dcd[1]+pc[1]*dcd[0])/r2pc);

        fdx=ddihe*((-pb[1]*dbc[2]+pb[2]*dbc[1])-pbpc*(-pc[1]*dbc[2]+pc[2]*dbc[1])/r2pc);
        fdy=ddihe*(( pb[0]*dbc[2]-pb[2]*dbc[0])-pbpc*( pc[0]*dbc[2]-pc[2]*dbc[0])/r2pc);
        fdz=ddihe*((-pb[0]*dbc[1]+pb[1]*dbc[0])-pbpc*(-pc[0]*dbc[1]+pc[1]*dbc[0])/r2pc);

        fx[i]+=fax;
        fy[i]+=fay;
        fz[i]+=faz;

        fx[j]+=-fax-fcx+fbx;
        fy[j]+=-fay-fcy+fby;
        fz[j]+=-faz-fcz+fbz;

        fx[k]+=fcx-fbx-fdx;
        fy[k]+=fcy-fby-fdy;
        fz[k]+=fcz-fbz-fdz;

        fx[l]+=fdx;
        fy[l]+=fdy;
        fz[l]+=fdz;

        /*stress[0]+=dab[0]*fax+dbc[0]*(fbx-fcx)-dcd[0]*fdx;
        stress[1]+=dab[1]*fax+dbc[1]*(fbx-fcx)-dcd[1]*fdx;
        stress[2]+=dab[2]*fax+dbc[2]*(fbx-fcx)-dcd[2]*fdx;
        stress[3]+=dab[1]*fay+dbc[1]*(fby-fcy)-dcd[1]*fdy;
        stress[4]+=dab[1]*faz+dbc[1]*(fbz-fcz)-dcd[1]*fdz;
        stress[5]+=dab[2]*faz+dbc[2]*(fbz-fcz)-dcd[2]*fdz; */

    }

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
    update_timer_end(TIMER_ENERGY_DIHE,__func__);
#endif

}

void improper_energy(const PARAM *param,const PARALLEL *parallel,ENERGY *ener,const PBC *box,
                     const DIHE impr[],const double *x,const double *y,const double *z,
                     double *fx,double *fy,double *fz)
{
    int i,j,k,l,ll;
    double edihe=0.,ddihe=0.;
    double cosp,sinp,phi;
    double /*rab,*/rbc,/*rcd,*/rpb,rpc,r2pb,r2pc,pbpc;
    double dab[3],dbc[3],dcd[3],/*dac[3],*/pb[3],pc[3]/*,stress[6]={0.}*/;
    double fax,fay,faz,fbx,fby,fbz,fcx,fcy,fcz,fdx,fdy,fdz;

#ifdef TIMER
    update_timer_begin(TIMER_ENERGY_UB,__func__);
#endif

    for(ll=parallel->fIpProc; ll<parallel->lIpProc; ll++)
    {
        i=impr[ll].a;
        j=impr[ll].b;
        k=impr[ll].c;
        l=impr[ll].d;

        dab[0]=x[i]-x[j];
        dab[1]=y[i]-y[j];
        dab[2]=z[i]-z[j];

        dist(box,dab);

        dbc[0]=x[j]-x[k];
        dbc[1]=y[j]-y[k];
        dbc[2]=z[j]-z[k];

        rbc=sqrt(dist(box,dbc));

        dcd[0]=x[k]-x[l];
        dcd[1]=y[k]-y[l];
        dcd[2]=z[k]-z[l];

        dist(box,dcd);

// construct first dihedral vector

        pb[0]=dab[1]*dbc[2]-dab[2]*dbc[1];
        pb[1]=dab[2]*dbc[0]-dab[0]*dbc[2];
        pb[2]=dab[0]*dbc[1]-dab[1]*dbc[0];

        r2pb=X2(pb[0])+X2(pb[1])+X2(pb[2]);
        rpb=sqrt(r2pb);

// construct second dihedral vector

        pc[0]=dbc[1]*dcd[2]-dbc[2]*dcd[1];
        pc[1]=dbc[2]*dcd[0]-dbc[0]*dcd[2];
        pc[2]=dbc[0]*dcd[1]-dbc[1]*dcd[0];

        r2pc=X2(pc[0])+X2(pc[1])+X2(pc[2]);
        rpc=sqrt(r2pc);

// determine dihedral angle

        pbpc=pb[0]*pc[0]+pb[1]*pc[1]+pb[2]*pc[2];
        cosp=pbpc/(rpb*rpc);

        sinp=(dbc[0]*(pc[1]*pb[2]-pc[2]*pb[1])+dbc[1]*(pb[0]*pc[2]-pb[2]*pc[0])+
              dbc[2]*(pc[0]*pb[1]-pc[1]*pb[0]))/(rpb*rpc*rbc);

        phi=atan2(sinp,cosp);

// avoid singularity in sinp

        if(sinp>=0.)
        {
            sinp=MAX(DBL_EPSILON,fabs(sinp));
        }
        else
        {
            sinp=-(MAX(DBL_EPSILON,fabs(sinp)));
        }

// calculate potential energy and scalar force term

        switch(impr[ll].type)
        {
        case DCOS: // cosine improper dihedral

            edihe=impr[ll].k*(1.+cos(impr[ll].mult*phi-impr[ll].phi0));
            ddihe=-impr[ll].k*impr[ll].mult*
                  sin(impr[ll].mult*phi-impr[ll].phi0)/(rpb*rpc*sinp);

            break;

        case DHARM: // harmonic improper dihedral

            phi=phi-impr[ll].phi0;
            phi=phi-nint(phi/TWOPI)*TWOPI;
            edihe=0.5*impr[ll].k*(phi*phi);
            ddihe=impr[ll].k*phi/(rpb*rpc*sinp);

            break;

        default:

            edihe=impr[ll].k*(1.+cos(impr[ll].mult*phi-impr[ll].phi0));
            ddihe=-impr[ll].k*impr[ll].mult*
                  sin(impr[ll].mult*phi-impr[ll].phi0)/(rpb*rpc*sinp);

            break;
        }

// calculate potential energy

        ener->impr+=edihe;

        fax=ddihe*((-pc[1]*dbc[2]+pc[2]*dbc[1])-pbpc*(-pb[1]*dbc[2]+pb[2]*dbc[1])/r2pb);
        fay=ddihe*(( pc[0]*dbc[2]-pc[2]*dbc[0])-pbpc*( pb[0]*dbc[2]-pb[2]*dbc[0])/r2pb);
        faz=ddihe*((-pc[0]*dbc[1]+pc[1]*dbc[0])-pbpc*(-pb[0]*dbc[1]+pb[1]*dbc[0])/r2pb);

        fcx=ddihe*((-pc[1]*dab[2]+pc[2]*dab[1])-pbpc*(-pb[1]*dab[2]+pb[2]*dab[1])/r2pb);
        fcy=ddihe*(( pc[0]*dab[2]-pc[2]*dab[0])-pbpc*( pb[0]*dab[2]-pb[2]*dab[0])/r2pb);
        fcz=ddihe*((-pc[0]*dab[1]+pc[1]*dab[0])-pbpc*(-pb[0]*dab[1]+pb[1]*dab[0])/r2pb);

        fbx=ddihe*((-pb[1]*dcd[2]+pb[2]*dcd[1])-pbpc*(-pc[1]*dcd[2]+pc[2]*dcd[1])/r2pc);
        fby=ddihe*(( pb[0]*dcd[2]-pb[2]*dcd[0])-pbpc*( pc[0]*dcd[2]-pc[2]*dcd[0])/r2pc);
        fbz=ddihe*((-pb[0]*dcd[1]+pb[1]*dcd[0])-pbpc*(-pc[0]*dcd[1]+pc[1]*dcd[0])/r2pc);

        fdx=ddihe*((-pb[1]*dbc[2]+pb[2]*dbc[1])-pbpc*(-pc[1]*dbc[2]+pc[2]*dbc[1])/r2pc);
        fdy=ddihe*(( pb[0]*dbc[2]-pb[2]*dbc[0])-pbpc*( pc[0]*dbc[2]-pc[2]*dbc[0])/r2pc);
        fdz=ddihe*((-pb[0]*dbc[1]+pb[1]*dbc[0])-pbpc*(-pc[0]*dbc[1]+pc[1]*dbc[0])/r2pc);

        fx[i]+=fax;
        fy[i]+=fay;
        fz[i]+=faz;

        fx[j]+=-fax-fcx+fbx;
        fy[j]+=-fay-fcy+fby;
        fz[j]+=-faz-fcz+fbz;

        fx[k]+=fcx-fbx-fdx;
        fy[k]+=fcy-fby-fdy;
        fz[k]+=fcz-fbz-fdz;

        fx[l]+=fdx;
        fy[l]+=fdy;
        fz[l]+=fdz;

        /*stress[0]+=dab[0]*fax+dbc[0]*(fbx-fcx)-dcd[0]*fdx;
        stress[1]+=dab[1]*fax+dbc[1]*(fbx-fcx)-dcd[1]*fdx;
        stress[2]+=dab[2]*fax+dbc[2]*(fbx-fcx)-dcd[2]*fdx;
        stress[3]+=dab[1]*fay+dbc[1]*(fby-fcy)-dcd[1]*fdy;
        stress[4]+=dab[1]*faz+dbc[1]*(fbz-fcz)-dcd[1]*fdz;
        stress[5]+=dab[2]*faz+dbc[2]*(fbz-fcz)-dcd[2]*fdz;*/

    }

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
    update_timer_end(TIMER_ENERGY_UB,__func__);
#endif

}
