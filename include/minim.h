#ifndef MINIMH_INCLUDED
#define MINIMH_INCLUDED

#ifdef	__cplusplus
extern "C" {
#endif

    void minimise(CTRL *ctrl,PARAM *param,ENERGY *ener,PBC *box,NEIGH *neigh,
                  ATOM atom[],BOND bond[],BOND ub[],ANGLE angle[],DIHE dihe[],
                  DIHE impr[]);

    void steepestDescent(CTRL *ctrl,PARAM *param,ENERGY *ener,PBC *box,NEIGH *neigh,
                         ATOM atom[],BOND bond[],BOND ub[],ANGLE angle[],DIHE dihe[],
                         DIHE impr[],double x[],double y[],double z[],double fx[],
                         double fy[],double fz[]);

    void conjugateGradients(CTRL *ctrl,PARAM *param,ENERGY *ener,PBC *box,NEIGH *neigh,
                            ATOM atom[],BOND bond[],BOND ub[],ANGLE angle[],DIHE dihe[],
                            DIHE impr[]);

#ifdef	__cplusplus
}
#endif

#endif
