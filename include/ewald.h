#ifndef EWALDH_INCLUDED
#define EWALDH_INCLUDED

#ifdef	__cplusplus
extern "C" {
#endif

    void init_ewald(CTRL *ctrl,PARAM *param,PARALLEL *parallel,EWALD *ewald,PBC *box);

    void ewald_free(EWALD *ewald);

    real ewald_rec(PARAM *param,PARALLEL *parallel,EWALD *ewald,PBC *box,const real x[],
                     const real y[],const real z[],real fx[],real fy[],real fz[],
                     const real q[],real stress[6],real *virEwaldRec,real dBuffer[]);

    real ewald_dir(EWALD *ewald,real *dEwaldDir,const real qel,
                     const real r,const real rt);

    real ewald_corr(EWALD *ewald,real *dEwaldCorr,const real qel,
                      const real r,const real rt);

    real ewald_dir14(PARAM *param,EWALD *ewald,real *dEwaldDir,const real qel,
                       const real r,const real rt);

    real ewald_corr14(PARAM *param,EWALD *ewald,real *dEwaldCorr,
                        const real qel,const real r,const real rt);

#ifdef	__cplusplus
}
#endif

#endif
