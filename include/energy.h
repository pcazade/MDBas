#ifndef ENERGYH_INCLUDED
#define ENERGYH_INCLUDED

#ifdef	__cplusplus
extern "C" {
#endif

    void init_energy_ptrs(CTRL *ctrl);

    void energy(CTRL *ctrl,PARAM *param,PARALLEL *parallel,ENERGY *ener,EWALD *ewald,PBC *box,NEIGH *neigh,
                BOND bond[],BOND ub[],ANGLE angle[],DIHE dihe[],DIHE impr[],
                const real x[],const real y[], const real z[],
                real vx[],real vy[], real vz[],real fx[],real fy[],
                real fz[],const real q[],const real eps[],const real sig[],
                const real eps14[],const real sig14[],const int frozen[],
                int **neighList,const int neighPair[],const int neighList14[],
                int **exclList,const int exclPair[],real dBuffer[]);

    void nonbond_energy(PARAM *param,PARALLEL *parallel,ENERGY *ener,PBC *box,
                        const real x[],const real y[],const real z[],real fx[],real fy[],
                        real fz[],const real q[],const real eps[],const real sig[],
                        int **neighList,const int neighPair[]);

    void nonbond14_energy(PARAM *param,PARALLEL *parallel,ENERGY *ener,PBC *box,NEIGH *neigh,
                          const real x[],const real y[],const real z[],real fx[],real fy[],
                          real fz[],const real q[],const real eps[],const real sig[],
                          const int neighList14[]);

    void ewald_energy(CTRL *ctrl,PARAM *param,PARALLEL *parallel,ENERGY *ener,EWALD *ewald,PBC *box,const real x[],
                      const real y[],const real z[],real fx[],real fy[],
                      real fz[],const real q[],const real eps[],const real sig[],
                      int **neighList,const int neighPair[],int **exclList,
                      const int exclPair[],real dBuffer[]);

    void ewald14_energy(PARAM *param,PARALLEL *parallel,ENERGY *ener,EWALD *ewald,PBC *box,NEIGH *neigh,const real x[],
                        const real y[],const real z[],real fx[],real fy[],
                        real fz[],const real q[],const real eps[],const real sig[],
                        const int neighList14[]);

    /* pointers to function for electrostatic energy functions */
    real (*ptr_coulomb)(const PARAM *param,real *delec,const real qel,
                          const real r2,const real rt);

    real (*ptr_coulomb14)(const PARAM *param,real *delec,const real qel,
                            const real r2,const real rt);

    /* pointers to function for vdw energy functions */
    real (*ptr_vdw)(const PARAM *param,real *dvdw,const real veps,
                      const real vsig,const real r2, const real rt);

    real (*ptr_vdw14)(const PARAM *param,real *dvdw,const real veps,
                        const real vsig,const real r2, const real rt);

#ifdef	__cplusplus
}
#endif

#endif
