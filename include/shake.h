#ifndef SHAKEH_INCLUDED
#define SHAKEH_INCLUDED

#ifdef	__cplusplus
extern "C" {
#endif

    void shake_allocate_arrays(const PARAM *param,const PARALLEL *parallel);

    void shake_free_arrays();

    void lf_shake(PARAM *param,PBC *box,CONSTRAINT constList[],PARALLEL *parallel,
                  real x[],real y[],real z[],
                  real ddx[],real ddy[],real ddz[],real rmass[],
                  int *nAtConst,real stress[6],real *virshake,real dBuffer[]);

    void vv_shake_r(PARAM *param,PBC *box,CONSTRAINT constList[],PARALLEL *parallel,
                    real x[],real y[],real z[],
                    real vx[],real vy[],real vz[],
                    real ddx[],real ddy[],real ddz[],real rmass[],
                    int *nAtConst,real stress[6],real *virshake,real dBuffer[]);

    void vv_shake_v(PARAM *param,CONSTRAINT constList[],PARALLEL *parallel,
                    real vx[],real vy[],real vz[],real ddx[],
                    real ddy[],real ddz[],real rmass[],int *nAtConst,real dBuffer[]);

#ifdef	__cplusplus
}
#endif

#endif
