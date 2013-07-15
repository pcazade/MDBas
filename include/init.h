#ifndef INITH_INCLUDED
#define INITH_INCLUDED

#ifdef	__cplusplus
extern "C" {
#endif

    void init_system(int *argc, char ***argv,IO *inout,CTRL *ctrl,PARAM *param,
                     PARALLEL *parallel,ENERGY *ener,BATH *bath,NEIGH *neigh,EWALD *ewald,
                     PBC *box,ATOM **atom,CONSTRAINT **constList,BOND **bond,ANGLE **angle,
                     DIHE **dihe,DIHE **impr,BOND **ub,double **x,double **y, double **z,
                     double **vx,double **vy,double **vz,double **fx,double **fy, double **fz,
                     double **mass,double **rmass,double **q,double **eps,double **sig,
                     double **eps14,double **sig14,int **frozen,int **nAtConst,int ***neighList,
                     int **neighPair,int **neighList14,int ***exclList,int **exclPair,
                     double **dBuffer,int **iBuffer);

    void init_variables(CTRL *ctrl,PARAM *param,PARALLEL *parallel,BATH *bath,NEIGH *neigh,
                        EWALD *ewald,PBC *box);

    void setup(CTRL *ctrl,PARAM *param,ATOM atom[],CONSTRAINT **constList,
               BOND **bond,ANGLE **angle,DIHE **dihe,DIHE **impr,BOND **ub,
               double mass[],double rmass[],int frozen[],int nAtConst[]);

    void init_vel(PARAM *param,PARALLEL *parallel,PBC *box,CONSTRAINT constList[],
                  double x[],double y[],double z[],double vx[],double vy[],
                  double vz[],double mass[],double rmass[],int frozen[],
                  int nAtConst[], double dBuffer[]);

    void init_constvel(PARAM *param,PBC *box,CONSTRAINT constList[],double x[],
                       double y[],double z[],double vx[],double vy[],double vz[],
                       double mass[],int nAtConst[]);

    void init_box(PBC *box);

    void free_all(CTRL *ctrl,PARAM *param, PARALLEL *parallel,EWALD *ewald,ATOM **atom,
                  CONSTRAINT **constList,BOND **bond,ANGLE **angle,DIHE **dihe,DIHE **impr,
                  BOND **ub,double **x,double **y, double **z,double **vx,double **vy,
                  double **vz,double **fx,double **fy, double **fz,double **mass,double **rmass,
                  double **q,double **eps,double **sig,double **eps14,double **sig14,int **frozen,
                  int **nAtConst,int ***neighList,int **neighPair,int **neighList14,
                  int ***exclList,int **exclPair,double **dBuffer,int **iBuffer);

#ifdef	__cplusplus
}
#endif

#endif
