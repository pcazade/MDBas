#ifndef VDWH_INCLUDED
#define VDWH_INCLUDED

#ifdef	__cplusplus
extern "C" {
#endif

    real vdw_none(const PARAM *param,real *dvdw,const real veps,
                    const real vsig,const real r2, const real rt);

    void vdw_full(const PARAM *param,PARALLEL *parallel, ENERGY *ener, const PBC *box,real *x,
                  real *y,real *z,real *fx, real *fy, real *fz,real *eps,real *sig,
                  int **exclList,int *exclPair);

    real vdw_switch(const PARAM *param,real *dvdw,const real veps,
                      const real vsig,const real r2, const real rt);

    real vdw14_none(const PARAM *param,real *dvdw,const real veps,
                      const real vsig,const real r2, const real rt);

    real vdw14_full(const PARAM *param,real *dvdw,const real veps,
                      const real vsig,const real r2, const real rt);

    real vdw14_switch(const PARAM *param,real *dvdw,const real veps,
                        const real vsig,const real r2, const real rt);

#ifdef	__cplusplus
}
#endif

#endif
