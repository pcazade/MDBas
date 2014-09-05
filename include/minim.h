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
                         DIHE impr[],real x[],real y[],real z[],real fx[],
                         real fy[],real fz[]);

    void conjugateGradients(CTRL *ctrl,PARAM *param,ENERGY *ener,PBC *box,NEIGH *neigh,
                            ATOM atom[],BOND bond[],BOND ub[],ANGLE angle[],DIHE dihe[],
                            DIHE impr[]);

#ifdef	__cplusplus
}
#endif

#endif
