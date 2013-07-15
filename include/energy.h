#ifndef ENERGYH_INCLUDED
#define ENERGYH_INCLUDED

#ifdef	__cplusplus
extern "C" {
#endif

    void init_energy_ptrs(CTRL *ctrl);

    void energy(CTRL *ctrl,PARAM *param,PARALLEL *parallel,ENERGY *ener,EWALD *ewald,PBC *box,NEIGH *neigh,
                BOND bond[],BOND ub[],ANGLE angle[],DIHE dihe[],DIHE impr[],
                const double restrict x[],const double restrict y[], const double restrict z[],
                double restrict vx[],double restrict vy[], double restrict vz[],double restrict fx[],double restrict fy[],
                double restrict fz[],const double restrict q[],const double restrict eps[],const double restrict sig[],
                const double restrict eps14[],const double restrict sig14[],const int restrict frozen[],
                const int restrict **neighList,const int restrict neighPair[],const int neighList14[],
                const int **exclList,const int restrict exclPair[],const double restrict dBuffer[]);

    void nonbond_energy(PARAM *param,PARALLEL *parallel,ENERGY *ener,PBC *box,
                        const double x[],const double y[],const double z[],double fx[],double fy[],
                        double fz[],const double q[],const double eps[],const double sig[],
                        int **neighList,const int neighPair[]);

    void nonbond14_energy(PARAM *param,PARALLEL *parallel,ENERGY *ener,PBC *box,NEIGH *neigh,
                          const double x[],const double y[],const double z[],double fx[],double fy[],
                          double fz[],const double q[],const double eps[],const double sig[],
                          const int neighList14[]);

    void ewald_energy(CTRL *ctrl,PARAM *param,PARALLEL *parallel,ENERGY *ener,EWALD *ewald,PBC *box,const double x[],
                      const double y[],const double z[],double fx[],double fy[],
                      double fz[],const double q[],const double eps[],const double sig[],
                      int **neighList,const int neighPair[],int **exclList,
                      const int exclPair[],double dBuffer[]);

    void ewald14_energy(PARAM *param,PARALLEL *parallel,ENERGY *ener,EWALD *ewald,PBC *box,NEIGH *neigh,const double x[],
                        const double y[],const double z[],double fx[],double fy[],
                        double fz[],const double q[],const double eps[],const double sig[],
                        const int neighList14[]);

    /* pointers to function for electrostatic energy functions */
    double (*ptr_coulomb)(const PARAM *param,double *delec,const double qel,
                          const double r2,const double rt);

    double (*ptr_coulomb14)(const PARAM *param,double *delec,const double qel,
                            const double r2,const double rt);

    /* pointers to function for vdw energy functions */
    double (*ptr_vdw)(const PARAM *param,double *dvdw,const double veps,
                      const double vsig,const double r2, const double rt);

    double (*ptr_vdw14)(const PARAM *param,double *dvdw,const double veps,
                        const double vsig,const double r2, const double rt);

#ifdef	__cplusplus
}
#endif

#endif
