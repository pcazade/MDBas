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

#include <stdlib.h>
#include <float.h>
#include <math.h>

#include <cuda.h>
#include <cuda_runtime.h>
#include <cufft.h>
#include <cuComplex.h>

#include "global.h"
#include "utils.h"
#include "memory.h"
#include "spme.h"
#include "kernels_utils.h"

#ifdef USING_MPI
#include "parallel.h"
#else
#include "serial.h"
#endif

static int newJob;

static real eEwaldself,systq;

static real *sx,*sy,*sz;

static real *bsp,*qsp;
static real **bsp1,**bsp2,**bsp3;
static real **bsd1,**bsd2,**bsd3;

static cplx *bspc1,*bspc2,*bspc3;
static cplx *epl1,*epl2,*epl3;

static cufftComplex *d_ftqsp;
static cufftHandle d_fft3d;

void init_spme(CTRL *ctrl,PARAM *param,PARALLEL *parallel,EWALD *ewald,PBC *box)
{
    int i,m1maxp2,m2maxp2,m3maxp2;

    newJob=1;

    ewald->prec=fmin( fabs( ewald->prec ) , 0.5 );
    ewald->tol=sqrt( fabs( log( ewald->prec * param->cutOff ) ) );

    if(!ctrl->keyAlpha)
        ewald->alpha=sqrt( fabs( log( ewald->prec * param->cutOff * ewald->tol ) ) ) / param->cutOff;

    ewald->tol1=sqrt( -log( ewald->prec * param->cutOff * X2( 2.0 * ewald->tol * ewald->alpha ) ) );

    if(!ctrl->keyMmax)
    {
        ewald->m1max=nint(0.25+box->pa*ewald->alpha*ewald->tol1/PI);
        ewald->m2max=nint(0.25+box->pb*ewald->alpha*ewald->tol1/PI);
        ewald->m3max=nint(0.25+box->pc*ewald->alpha*ewald->tol1/PI);
    }

    m1maxp2=1;
    while( (ewald->m1max>m1maxp2) && (m1maxp2<256) )
    {
        m1maxp2*=2;
    }
    ewald->m1max=2*m1maxp2;

    m2maxp2=1;
    while( (ewald->m2max>m2maxp2) && (m2maxp2<256) )
    {
        m2maxp2*=2;
    }
    ewald->m2max=2*m2maxp2;

    m3maxp2=1;
    while( (ewald->m3max>m3maxp2) && (m3maxp2<256) )
    {
        m3maxp2*=2;
    }
    ewald->m3max=2*m3maxp2;

    epl1=(cplx*)my_malloc(ewald->m1max*sizeof(*epl1));
    epl2=(cplx*)my_malloc(ewald->m2max*sizeof(*epl2));
    epl3=(cplx*)my_malloc(ewald->m3max*sizeof(*epl3));

    bspc1=(cplx*)my_malloc(ewald->m1max*sizeof(*bspc1));
    bspc2=(cplx*)my_malloc(ewald->m2max*sizeof(*bspc2));
    bspc3=(cplx*)my_malloc(ewald->m3max*sizeof(*bspc3));

    sx=(real*)my_malloc(parallel->maxAtProc*sizeof(*sx));
    sy=(real*)my_malloc(parallel->maxAtProc*sizeof(*sy));
    sz=(real*)my_malloc(parallel->maxAtProc*sizeof(*sz));

    bsp=(real*)my_malloc(ewald->nbsp*sizeof(*bsp));

    bsp1=(real**)my_malloc(parallel->maxAtProc*sizeof(*bsp1));
    bsp2=(real**)my_malloc(parallel->maxAtProc*sizeof(*bsp2));
    bsp3=(real**)my_malloc(parallel->maxAtProc*sizeof(*bsp3));

    bsd1=(real**)my_malloc(parallel->maxAtProc*sizeof(*bsd1));
    bsd2=(real**)my_malloc(parallel->maxAtProc*sizeof(*bsd2));
    bsd3=(real**)my_malloc(parallel->maxAtProc*sizeof(*bsd3));

    for(i=0; i<parallel->maxAtProc; i++)
    {
        bsp1[i]=(real*)my_malloc(ewald->nbsp*sizeof(**bsp1));
        bsp2[i]=(real*)my_malloc(ewald->nbsp*sizeof(**bsp2));
        bsp3[i]=(real*)my_malloc(ewald->nbsp*sizeof(**bsp3));

        bsd1[i]=(real*)my_malloc(ewald->nbsp*sizeof(**bsd1));
        bsd2[i]=(real*)my_malloc(ewald->nbsp*sizeof(**bsd2));
        bsd3[i]=(real*)my_malloc(ewald->nbsp*sizeof(**bsd3));
    }

    ewald->mmax=ewald->m1max*ewald->m2max*ewald->m3max;

    qsp=(real*)my_malloc(ewald->mmax*sizeof(*qsp));
    d_ftqsp=(cplx*)cudaMalloc(ewald->mmax*sizeof(*d_ftqsp));

}

void spme_free(PARALLEL *parallel)
{
    int i;

    free(epl1);
    free(epl2);
    free(epl3);

    free(bspc1);
    free(bspc2);
    free(bspc3);

    free(sx);
    free(sy);
    free(sz);

    free(bsp);

    for(i=0; i<parallel->maxAtProc; i++)
    {
        free(bsp1[i]);
        free(bsp2[i]);
        free(bsp3[i]);

        free(bsd1[i]);
        free(bsd2[i]);
        free(bsd3[i]);
    }

    free(bsp1);
    free(bsp2);
    free(bsp3);

    free(bsd1);
    free(bsd2);
    free(bsd3);

    free(qsp);

    cudaFree(d_ftqsp);
    cufftDestroy(d_fft3d);

}

void epl_cplx(EWALD *ewald)
{

    int i,hm1max,hm2max,hm3max;

    hm1max=ewald->m1max/2;
    hm2max=ewald->m2max/2;
    hm3max=ewald->m3max/2;

    epl1[0]=cudaComp(1.0);

    for(i=1; i<=hm1max; i++)
    {
        epl1[i]=cudaExpI(TWOPI*(real)i/(real)ewald->m1max);
        epl1[ewald->m1max-i]=cudaConj(epl1[i]);
    }

    epl2[0]=cudaComp(1.0);

    for(i=1; i<=hm2max; i++)
    {
        epl2[i]=cudaExpI(TWOPI*(real)i/(real)ewald->m2max);
        epl2[ewald->m2max-i]=cudaConj(epl2[i]);
    }

    epl3[0]=cudaComp(1.0);

    for(i=1; i<=hm3max; i++)
    {
        epl3[i]=cudaExpI(TWOPI*(real)i/(real)ewald->m3max);
        epl3[ewald->m3max-i]=cudaConj(epl3[i]);
    }

}

void bspcoef(EWALD *ewald)
{

    int i,j,k;
    cplx coeff;

    bsp[0]=0.0;
    bsp[1]=1.0;

    for(k=2; k<ewald->nbsp; k++)
    {

        bsp[k]=0.0;

        for(j=k; j>0; j--)
        {
            bsp[j] = ( (real)j*bsp[j]+(real)(k+1-j)*bsp[j-1] ) / ( (real)k );
        }
    }

    for(i=0; i<ewald->m1max; i++)
    {

        coeff=cudaComp(0.0);

        for(k=0; k<ewald->nbsp-1; k++)
        {
            coeff=cudaAdd(coeff,cudaMul(cudaComp(bsp[k+1]),epl1[( (i*k) % ewald->m1max )]));
        }

        bspc1[i]=cudaDiv(epl1[( ( i* (ewald->nbsp-1) ) % ewald->m1max )],coeff);

    }

    for(i=0; i<ewald->m2max; i++)
    {

        coeff=cudaComp(0.0);

        for(k=0; k<ewald->nbsp-1; k++)
        {
            coeff=cudaAdd(coeff,cudaMul(cudaComp(bsp[k+1]),epl2[( (i*k) % ewald->m2max )]));
        }

        bspc2[i]=cudaDiv(epl2[( ( i* (ewald->nbsp-1) ) % ewald->m2max )],coeff);

    }

    for(i=0; i<ewald->m3max; i++)
    {

        coeff=0.+I*0.0;

        for(k=0; k<ewald->nbsp-1; k++)
        {
            coeff=cudaAdd(coeff,cudaMul(cudaComp(bsp[k+1]),epl3[( (i*k) % ewald->m3max )]));
        }

        bspc3[i]=cudaDiv(epl3[( ( i* (ewald->nbsp-1) ) % ewald->m3max )],coeff);

    }

}

void bspgen(PARALLEL *parallel,EWALD *ewald)
{

    int i,j,k;
    real tsx,tsy,tsz;

    for(i=0; i<parallel->nAtProc; i++)
    {

        bsd1[i][0]=1.0;
        bsd2[i][0]=1.0;
        bsd3[i][0]=1.0;

        bsd1[i][1]=-1.0;
        bsd2[i][1]=-1.0;
        bsd3[i][1]=-1.0;

        bsp1[i][0]=sx[i]-(int)sx[i];
        bsp2[i][0]=sy[i]-(int)sy[i];
        bsp3[i][0]=sz[i]-(int)sz[i];

        bsp1[i][1]=1.0-bsp1[i][0];
        bsp2[i][1]=1.0-bsp2[i][0];
        bsp3[i][1]=1.0-bsp3[i][0];

    }

    for(k=2; k<ewald->nbsp; k++)
    {

        for(i=0; i<parallel->nAtProc; i++)
        {
            bsp1[i][k]=0.0;
            bsp2[i][k]=0.0;
            bsp3[i][k]=0.0;
        }

        for(j=k; j>0; j--)
        {

            if( k == (ewald->nbsp-1) )
            {
                for(i=0; i<parallel->nAtProc; i++)
                {
                    bsd1[i][j]=bsp1[i][j]-bsp1[i][j-1];
                    bsd2[i][j]=bsp2[i][j]-bsp2[i][j-1];
                    bsd3[i][j]=bsp3[i][j]-bsp3[i][j-1];
                }
            }

            for(i=0; i<parallel->nAtProc; i++)
            {
                tsx=sx[i]+(real)j-(int)sx[i];
                tsy=sy[i]+(real)j-(int)sy[i];
                tsz=sz[i]+(real)j-(int)sz[i];

                bsp1[i][j]=(tsx*bsp1[i][j]+((real)(k+1)-tsx)*bsp1[i][j-1])/((real)k);
                bsp2[i][j]=(tsy*bsp2[i][j]+((real)(k+1)-tsy)*bsp2[i][j-1])/((real)k);
                bsp3[i][j]=(tsz*bsp3[i][j]+((real)(k+1)-tsz)*bsp3[i][j-1])/((real)k);
            }

        }

        if( k == (ewald->nbsp-1) )
        {
            for(i=0; i<parallel->nAtProc; i++)
            {
                bsd1[i][0]=bsp1[i][0];
                bsd2[i][0]=bsp2[i][0];
                bsd3[i][0]=bsp3[i][0];
            }
        }

        for(i=0; i<parallel->nAtProc; i++)
        {
            tsx=sx[i]-(int)sx[i];
            tsy=sy[i]-(int)sy[i];
            tsz=sz[i]-(int)sz[i];

            bsp1[i][0]=tsx*bsp1[i][0]/((real)k);
            bsp2[i][0]=tsy*bsp2[i][0]/((real)k);
            bsp3[i][0]=tsz*bsp3[i][0]/((real)k);
        }

    }

}

__global__ void mesh(real d_sx[],real d_x[],real d_y[],real d_z[],
		     real u1,real u2,real u3,int fAtProc,int lAtProc,int mmax)
{
  int ll=blockIdx.x * blockDim.x + threadIdx.x;
  int i=ll+fAtProc;
  
  if(i<lAtProc)
    d_sx[ll]=(real)mmax*(d_x[i]*u1+d_y[i]*u2+d_z[i]*u3+0.5);
}

real spme_energy(PARAM *param,PARALLEL *parallel,EWALD *ewald,PBC *box,const real x[],
                   const real y[],const real z[],real fx[],real fy[],real fz[],
                   const real q[],real stress[6],real *virEwaldRec,real dBuffer[])
{

    int i,ii,j,jj,k,kk,l,ll;
    int m1,m2,m3;
    int hm1max,hm2max,hm3max;

    cplx cam,etmp;

    real vam,qtmp,fact1;
    real tt,bm1,bm2,bm3;
    real rm,rrm,rmx,rmy,rmz;
    real rm1x,rm1y,rm1z,rm2x,rm2y,rm2z;
    real recCutOff,recCutOff2,rAlpha2,rVol;
    real eEwaldRec,eNonNeutral;
    real fbx,fby,fbz;
    real fm[3];
    
    real d_x,d_y,d_z;
    cplx *d_bspc1,*d_bspc2,*d_bspc3;
    cplx *d_bsp1,*d_bsp2,*d_bsp3;
    cplx *d_bsd1,*d_bsd2,*d_bsd3;
    
    dim3 threadsPerBlock(nThrdsX,nThrdsY);
    dim3 nAtBlocks(parallel->nAtProc/threadsPerBlock.x,parallel->nAtProc/threadsPerBlock.y);

    if(newJob)
    {
        newJob=0;

        systq=0.;
        eEwaldself=0.;

        for(i=parallel->fAtProc; i<parallel->lAtProc; i++)
        {
            systq+=q[i];
            eEwaldself+=X2(q[i]);
        }

        if(parallel->nProc>1)
        {
            dBuffer[0]=systq;
            dBuffer[1]=eEwaldself;
            sum_double_para(dBuffer,&(dBuffer[3]),2);
            systq=dBuffer[0];
            eEwaldself=dBuffer[1];
        }

        eEwaldself=-param->chargeConst*ewald->alpha*eEwaldself/SQRTPI;

//     Calculate complex exponentials
        epl_cplx(ewald);

//     Calculate the coefficient of the B-spline
        bspcoef(ewald);
	
# ifdef DOUBLE_CUDA
	cufftPlan3d(&d_fft3d,ewald->m1max,ewald->m2max,ewald->m3max,CUFFT_Z2Z);
#else
	cufftPlan3d(&d_fft3d,ewald->m1max,ewald->m2max,ewald->m3max,CUFFT_C2C);
#endif

    } //   End if(newJob)

    hm1max=ewald->m1max/2;
    hm2max=ewald->m2max/2;
    hm3max=ewald->m3max/2;

    eEwaldRec=0.;
    *virEwaldRec=0.;

    stress[0]=0.;
    stress[1]=0.;
    stress[2]=0.;
    stress[3]=0.;
    stress[4]=0.;
    stress[5]=0.;

    rVol=TWOPI/box->vol;
    rAlpha2=-0.25/(X2(ewald->alpha));

    //   Set the cutoff in reciprocal space
    recCutOff=fmin( ( (real)ewald->m1max*box->u ) , ( (real)ewald->m2max*box->v ) );
    recCutOff=fmin( recCutOff , (real)ewald->m3max*box->w );
    recCutOff=recCutOff*1.05*TWOPI;
    recCutOff2=X2(recCutOff);
    
//  Memory allocation on the GPU
    
    cudaMalloc((void**) &d_x,param->nAtom*sizeof(real));
    cudaMalloc((void**) &d_y,param->nAtom*sizeof(real));
    cudaMalloc((void**) &d_z,param->nAtom*sizeof(real));
    
    cudaMalloc((void**) &d_bspc1,ewald->m1max*sizeof(cplx));
    cudaMalloc((void**) &d_bspc2,ewald->m2max*sizeof(cplx));
    cudaMalloc((void**) &d_bspc3,ewald->m3max*sizeof(cplx));
    
    cudaMalloc((void**) &d_bsp1,parallel->maxAtProc*ewald->nbsp*sizeof(cplx));
    cudaMalloc((void**) &d_bsp2,parallel->maxAtProc*ewald->nbsp*sizeof(cplx));
    cudaMalloc((void**) &d_bsp3,parallel->maxAtProc*ewald->nbsp*sizeof(cplx));
    
    cudaMalloc((void**) &d_bsd1,parallel->maxAtProc*ewald->nbsp*sizeof(cplx));
    cudaMalloc((void**) &d_bsd2,parallel->maxAtProc*ewald->nbsp*sizeof(cplx));
    cudaMalloc((void**) &d_bsd3,parallel->maxAtProc*ewald->nbsp*sizeof(cplx));
    
// Copy arrays from host to device
  
    cudaMemcpy(d_x,x,param->nAtom*sizeof(real),cudaMemcpyHostToDevice);
    cudaMemcpy(d_y,y,param->nAtom*sizeof(real),cudaMemcpyHostToDevice);
    cudaMemcpy(d_z,z,param->nAtom*sizeof(real),cudaMemcpyHostToDevice);
    
    cudaMemcpy(d_bspc1,bspc1,param->nAtom*sizeof(cplx),cudaMemcpyHostToDevice);
    cudaMemcpy(d_bspc2,bspc2,param->nAtom*sizeof(cplx),cudaMemcpyHostToDevice);
    cudaMemcpy(d_bspc3,bspc3,param->nAtom*sizeof(cplx),cudaMemcpyHostToDevice);
    
//   Address atoms to the cells of the mesh [0...mimax]
    
    mesh<<<(parallel->nAtProc/nThrds),nThrds>>>(d_sx,d_x,d_y,d_z,box->u1,box->u2,box->u3,parallel->fAtProc,parallel->lAtProc,ewald->m1max);
    mesh<<<(parallel->nAtProc/nThrds),nThrds>>>(d_sy,d_x,d_y,d_z,box->v1,box->v2,box->v3,parallel->fAtProc,parallel->lAtProc,ewald->m2max);
    mesh<<<(parallel->nAtProc/nThrds),nThrds>>>(d_sz,d_x,d_y,d_z,box->w1,box->w2,box->w3,parallel->fAtProc,parallel->lAtProc,ewald->m3max);

//   Construct the B-splines
    bspgen(parallel,ewald);
    

//   Initialise charge array Q(k1,k2,k3)
    for(i=0; i<ewald->mmax; i++)
        qsp[i]=0.0;

//   Fill charge array Q(k1,k2,k3)
    ll=0;
    for(l=parallel->fAtProc; l<parallel->lAtProc; l++)
    {
        for(i=0; i<ewald->nbsp; i++)
        {
            ii=(int)sx[ll]-i;

            if(ii>=ewald->m1max)
                ii=0;

            if(ii<0)
                ii+=ewald->m1max;

            m1=ewald->m2max*ii;

            for(j=0; j<ewald->nbsp; j++)
            {
                jj=(int)sy[ll]-j;

                if(jj>=ewald->m2max)
                    jj=0;

                if(jj<0)
                    jj+=ewald->m2max;

                m2=ewald->m3max*(jj+m1);

                for(k=0; k<ewald->nbsp; k++)
                {
                    kk=(int)sz[ll]-k;

                    if(kk>=ewald->m3max)
                        kk=0;

                    if(kk<0)
                        kk+=ewald->m3max;

                    m3=kk+m2;
                    qsp[m3]+=q[l]*bsp1[ll][i]*bsp2[ll][j]*bsp3[ll][k];

                }
            }
        }
        ll++;
    }

    if(parallel->nProc>1)
        sum_double_para(qsp,dBuffer,ewald->mmax);

    for(m3=0; m3<ewald->mmax; m3++)
    {
        d_ftqsp[m3]=cudaComp(qsp[m3]);
    }

//   Perform the Fourier Transform of Q(k1,k2,k3)
# ifdef DOUBLE_CUDA
    cufftExecZ2Z(d_fft3d,d_ftqsp,d_ftqsp,CUFFT_INVERSE);
#else
    cufftExecC2C(d_fft3d,d_ftqsp,d_ftqsp,CUFFT_INVERSE);
#endif
    
    for(i=0; i<ewald->m1max; i++)
    {
        ii=i;

        if(i>hm1max)
            ii=i-ewald->m1max;

        m1=ewald->m2max*i;

        tt=TWOPI*(real)ii;

        rm1x=tt*box->u1;
        rm1y=tt*box->u2;
        rm1z=tt*box->u3;

        bm1=creal(bspc1[i]*conj(bspc1[i]));

        for(j=0; j<ewald->m2max; j++)
        {
            jj=j;

            if(j>hm2max)
                jj=j-ewald->m2max;

            m2=ewald->m3max*(j+m1);

            tt=TWOPI*(real)jj;

            rm2x=rm1x+(tt*box->v1);
            rm2y=rm1y+(tt*box->v2);
            rm2z=rm1z+(tt*box->v3);

            bm2=bm1*creal(bspc2[j]*conj(bspc2[j]));

            for(k=0; k<ewald->m3max; k++)
            {

                kk=k;

                if(k>hm3max)
                    kk=k-ewald->m3max;

                m3=k+m2;

                tt=TWOPI*(real)kk;

                rmx=rm2x+(tt*box->w1);
                rmy=rm2y+(tt*box->w2);
                rmz=rm2z+(tt*box->w3);

                bm3=bm2*creal(bspc3[k]*conj(bspc3[k]));

                rm=X2(rmx)+X2(rmy)+X2(rmz);

                if( (rm>DBL_EPSILON) && (rm<=recCutOff2 ) )
                {
                    rrm=1.0/rm;

                    cam=bm3*exp(rAlpha2*rm)*rrm*ftqsp[m3];

                    vam=2.0*(rrm-rAlpha2)*creal(cam*conj(ftqsp[m3]));

                    stress[0]-=vam*rmx*rmx;
                    stress[1]-=vam*rmx*rmy;
                    stress[2]-=vam*rmx*rmz;
                    stress[3]-=vam*rmy*rmy;
                    stress[4]-=vam*rmy*rmz;
                    stress[5]-=vam*rmz*rmz;

                    ftqsp[m3]=cam;

                }
                else
                {
                    ftqsp[m3]=0.0+I*0.0;
                }
            }
        }
    }

    /**   Beginning of the forces calculation section   */

# ifdef DOUBLE_CUDA
    cufftExecZ2Z(d_fft3d,d_ftqsp,d_ftqsp,CUFFT_FORWARD);
#else
    cufftExecC2C(d_fft3d,d_ftqsp,d_ftqsp,CUFFT_FORWARD);
#endif

    fact1=-2.0*rVol*param->chargeConst;

    ll=0;
    for(l=parallel->fAtProc; l<parallel->lAtProc; l++)
    {
        for(i=0; i<ewald->nbsp; i++)
        {
            ii=(int)sx[ll]-i;

            if(ii>=ewald->m1max)
                ii=0;

            if(ii<0)
                ii+=ewald->m1max;

            m1=ewald->m2max*ii;

            for(j=0; j<ewald->nbsp; j++)
            {
                jj=(int)sy[ll]-j;

                if(jj>=ewald->m2max)
                    jj=0;

                if(jj<0)
                    jj+=ewald->m2max;

                m2=ewald->m3max*(jj+m1);

                for(k=0; k<ewald->nbsp; k++)
                {
                    kk=(int)sz[ll]-k;

                    if(kk>=ewald->m3max)
                        kk=0;

                    if(kk<0)
                        kk+=ewald->m3max;

                    m3=kk+m2;

                    qtmp=creal(ftqsp[m3]);

                    fbx=qtmp*bsd1[ll][i]*bsp2[ll][j]*bsp3[ll][k]*(real)ewald->m1max;
                    fby=qtmp*bsp1[ll][i]*bsd2[ll][j]*bsp3[ll][k]*(real)ewald->m2max;
                    fbz=qtmp*bsp1[ll][i]*bsp2[ll][j]*bsd3[ll][k]*(real)ewald->m3max;

                    fx[l]+=fact1*q[l]*(fbx*box->u1+fby*box->v1+fbz*box->w1);
                    fy[l]+=fact1*q[l]*(fbx*box->u2+fby*box->v2+fbz*box->w2);
                    fz[l]+=fact1*q[l]*(fbx*box->u3+fby*box->v3+fbz*box->w3);

                }
            }
        }
        ll++;
    }

//   Set the sum of the forces to 0

    fm[0]=0.0;
    fm[1]=0.0;
    fm[2]=0.0;

    for(l=parallel->fAtProc; l<parallel->lAtProc; l++)
    {
        fm[0]+=fx[l];
        fm[1]+=fy[l];
        fm[2]+=fz[l];
    }

    if(parallel->nProc>1)
        sum_real_para(fm,dBuffer,3);

    fm[0]/=(real)param->nAtom;
    fm[1]/=(real)param->nAtom;
    fm[2]/=(real)param->nAtom;

    for(l=parallel->fAtProc; l<parallel->lAtProc; l++)
    {
        fx[l]-=fm[0];
        fy[l]-=fm[1];
        fz[l]-=fm[2];
    }

    /**   End of the forces calculation section   */

    etmp=0.0+I*0.0;
    for(m3=0; m3<ewald->mmax; m3++)
    {
        etmp+=ftqsp[m3]*qsp[m3];
    }

    eNonNeutral=-(0.5*PI*param->chargeConst)*(X2(systq/ewald->alpha)/box->vol)/(real)parallel->nProc;

    fact1=rVol*param->chargeConst/(real)parallel->nProc;

    eEwaldRec=creal(etmp);

    stress[0]=fact1*(stress[0]+eEwaldRec)+eNonNeutral;
    stress[1]=fact1*stress[1];
    stress[2]=fact1*stress[2];
    stress[3]=fact1*(stress[3]+eEwaldRec)+eNonNeutral;
    stress[4]=fact1*stress[4];
    stress[5]=fact1*(stress[5]+eEwaldRec)+eNonNeutral;

    eEwaldRec=fact1*eEwaldRec+eNonNeutral+eEwaldself/(real)parallel->nProc;

    *virEwaldRec=-(stress[0]+stress[3]+stress[5]);

    return eEwaldRec;

}

