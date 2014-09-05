#ifndef INTEGRATEH_INCLUDED
#define INTEGRATEH_INCLUDED

#ifdef	__cplusplus
extern "C" {
#endif

// called once at the beginning of simulation for allocating working arrays used by integration
    void integrators_allocate_arrays(CTRL *ctrl, PARAM *param, PARALLEL *parallel);

// called once at the end of simulation for freeing working arrays used by integration
    void integrators_free_arrays(CTRL *ctrl, PARAM *param);

    void lf_integrate(CTRL *ctrl,PARAM *param,ENERGY *ener,PBC *box,
                      BATH *bath,CONSTRAINT constList[],PARALLEL *parallel,
                      real *x,real *y,real *z,
                      real *vx,real *vy,real *vz,
                      real *fx,real *fy,real *fz,
                      real *mass,real *rmass,int *nAtConst,real dBuffer[]);

    void lf_nve(PARAM *param,ENERGY *ener,PBC *box,CONSTRAINT constList[],PARALLEL *parallel,
                real *x,real *y,real *z,
                real *vx,real *vy,real *vz,
                real *fx,real *fy,real *fz,
                real *mass,real *rmass,int *nAtConst,real dBuffer[]);

    void lf_nvt_b(PARAM *param,ENERGY *ener,PBC *box,BATH *bath,CONSTRAINT constList[],PARALLEL *parallel,
                  real *x,real *y,real *z,
                  real *vx,real *vy,real *vz,
                  real *fx,real *fy,real *fz,
                  real *mass,real *rmass,int *nAtConst,real dBuffer[]);

    void lf_npt_b(PARAM *param,ENERGY *ener,PBC *box,BATH *bath,CONSTRAINT constList[],PARALLEL *parallel,
                  real *x,real *y,real *z,
                  real *vx,real *vy,real *vz,
                  real *fx,real *fy,real *fz,
                  real *mass,real *rmass,int *nAtConst,real dBuffer[]);

    void lf_nvt_h(PARAM *param,ENERGY *ener,PBC *box,BATH *bath,CONSTRAINT constList[],PARALLEL *parallel,
                  real *x,real *y,real *z,
                  real *vx,real *vy,real *vz,
                  real *fx,real *fy,real *fz,
                  real *mass,real *rmass,int *nAtConst,real dBuffer[]);

    void lf_npt_h(PARAM *param,ENERGY *ener,PBC *box,BATH *bath,CONSTRAINT constList[],PARALLEL *parallel,
                  real *x,real *y,real *z,
                  real *vx,real *vy,real *vz,
                  real *fx,real *fy,real *fz,
                  real *mass,real *rmass,int *nAtConst,real dBuffer[]);


    void vv_integrate(CTRL *ctrl,PARAM *param,ENERGY *ener,PBC *box,BATH *bath,
                      CONSTRAINT constList[],PARALLEL *parallel,real *x,real *y,
                      real *z,real *vx,real *vy,real *vz,
                      real *fx,real *fy,real *fz,
                      real *mass,real *rmass,int *nAtConst,real dBuffer[],int stage);

    void vv_nve(PARAM *param,ENERGY *ener,PBC *box,CONSTRAINT constList[],PARALLEL *parallel,
                real *x,real *y,real *z,real *vx,real *vy,real *vz,
                real *fx,real *fy,real *fz,
                real *mass,real *rmass,int *nAtConst,real dBuffer[],int stage);

    void vv_nvt_b(PARAM *param,ENERGY *ener,PBC *box,BATH *bath,CONSTRAINT constList[],PARALLEL *parallel,
                  real *x,real *y,real *z,real *vx,real *vy,real *vz,
                  real *fx,real *fy,real *fz,
                  real *mass,real *rmass,int *nAtConst,real dBuffer[],int stage);

    void vv_npt_b(PARAM *param,ENERGY *ener,PBC *box,BATH *bath,CONSTRAINT constList[],PARALLEL *parallel,
                  real *x,real *y,real *z,real *vx,real *vy,real *vz,
                  real *fx,real *fy,real *fz,
                  real *mass,real *rmass,int *nAtConst,real dBuffer[],int stage);

    void vv_nvt_h(PARAM *param,ENERGY *ener,PBC *box,BATH *bath,CONSTRAINT constList[],PARALLEL *parallel,
                  real *x,real *y,real *z,real *vx,real *vy,real *vz,
                  real *fx,real *fy,real *fz,
                  real *mass,real *rmass,int *nAtConst,real dBuffer[],int stage);

    void vv_npt_h(PARAM *param,ENERGY *ener,PBC *box,BATH *bath,CONSTRAINT constList[],PARALLEL *parallel,
                  real *x,real *y,real *z,real *vx,real *vy,real *vz,
                  real *fx,real *fy,real *fz,
                  real *mass,real *rmass,int *nAtConst,real dBuffer[],int stage);

#ifdef	__cplusplus
}
#endif

#endif
