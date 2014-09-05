#ifndef INITH_INCLUDED
#define INITH_INCLUDED

#ifdef	__cplusplus
extern "C" {
#endif

    void init_system(int *argc, char ***argv,IO *inout,CTRL *ctrl,PARAM *param,
                     PARALLEL *parallel,ENERGY *ener,BATH *bath,NEIGH *neigh,EWALD *ewald,
                     PBC *box,ATOM **atom,CONSTRAINT **constList,BOND **bond,ANGLE **angle,
                     DIHE **dihe,DIHE **impr,BOND **ub,real **x,real **y, real **z,
                     real **vx,real **vy,real **vz,real **fx,real **fy, real **fz,
                     real **mass,real **rmass,real **q,real **eps,real **sig,
                     real **eps14,real **sig14,int **frozen,int **nAtConst,int ***neighList,
                     int **neighPair,int **neighList14,int ***exclList,int **exclPair,
                     real **dBuffer,int **iBuffer);

    void init_variables(CTRL *ctrl,PARAM *param,PARALLEL *parallel,BATH *bath,NEIGH *neigh,
                        EWALD *ewald,PBC *box);

    void setup(CTRL *ctrl,PARAM *param,ATOM atom[],CONSTRAINT **constList,
               BOND **bond,ANGLE **angle,DIHE **dihe,DIHE **impr,BOND **ub,
               real mass[],real rmass[],int frozen[],int nAtConst[]);

    void init_vel(PARAM *param,PARALLEL *parallel,PBC *box,CONSTRAINT constList[],
                  real x[],real y[],real z[],real vx[],real vy[],
                  real vz[],real mass[],real rmass[],int frozen[],
                  int nAtConst[], real dBuffer[]);

    void init_constvel(PARAM *param,PBC *box,CONSTRAINT constList[],real x[],
                       real y[],real z[],real vx[],real vy[],real vz[],
                       real mass[],int nAtConst[]);

    void init_box(PBC *box);

    void free_all(CTRL *ctrl,PARAM *param, PARALLEL *parallel,EWALD *ewald,ATOM **atom,
                  CONSTRAINT **constList,BOND **bond,ANGLE **angle,DIHE **dihe,DIHE **impr,
                  BOND **ub,real **x,real **y, real **z,real **vx,real **vy,
                  real **vz,real **fx,real **fy, real **fz,real **mass,real **rmass,
                  real **q,real **eps,real **sig,real **eps14,real **sig14,int **frozen,
                  int **nAtConst,int ***neighList,int **neighPair,int **neighList14,
                  int ***exclList,int **exclPair,real **dBuffer,int **iBuffer);

#ifdef	__cplusplus
}
#endif

#endif
