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

// real spme_energy(PARAM *param,EWALD *ewald,PBC *box,const real x[],const real y[],const real z[],
// 	           real fx[],real fy[],real fz[],const real q[],real stress[6],real *virEwaldRec);

// real spme_energy(PARAM *param,EWALD *ewald,PBC *box,real x[],real y[],real z[],
// 	           real fx[],real fy[],real fz[],real q[],real stress[6],real *virEwaldRec);

    real spme_energy(PARAM *param,PARALLEL *parallel,EWALD *ewald,PBC *box,const real x[],
                       const real y[],const real z[],real fx[],real fy[],real fz[],
                       const real q[],real stress[6],real *virEwaldRec,real dBuffer[]);
#ifdef	__cplusplus
}
#endif

#endif
