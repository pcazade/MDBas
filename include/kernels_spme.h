#ifndef SPMEH_INCLUDED
#define SPMEH_INCLUDED

#ifdef	__cplusplus
extern "C" {
#endif

    void init_spme(CTRL *ctrl,PARAM *param,PARALLEL *parallel,EWALD *ewald,PBC *box);

    void spme_free(PARALLEL *parallel);

    void epl_cplx(EWALD *ewald);

    void bspcoef(EWALD *ewald);

    void bspgen(PARALLEL *parallel,EWALD *ewald);
    
    __global__ void mesh(real d_sx[],real d_x[],real d_y[],real d_z[],
			 real u1,real u2,real u3,int fAtProc,int lAtProc,int mmax);

    real spme_energy(PARAM *param,PARALLEL *parallel,EWALD *ewald,PBC *box,const real x[],
                       const real y[],const real z[],real fx[],real fy[],real fz[],
                       const real q[],real stress[6],real *virEwaldRec,real dBuffer[]);
#ifdef	__cplusplus
}
#endif

#endif
