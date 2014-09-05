#ifndef IOH_INCLUDED
#define IOH_INCLUDED

#ifdef	__cplusplus
extern "C" {
#endif

    void read_command_line(int *argc, char ***argv,IO *inout,PARALLEL *parallel);

    void read_SIMU(IO *inout,CTRL *ctrl,PARAM *param,BATH *bath,NEIGH *neigh,EWALD *ewald,PBC *box);

    void read_CONF(IO *inout,PARAM *param,ATOM **atom,real **x,real **y, real **z);

    void read_rest(IO *inout,PARAM *param,ENERGY *ener,BATH *bath,ATOM **atom,
                   real **x,real **y, real **z,real **vx,real **vy,
                   real **vz,real **fx,real **fy, real **fz);

    void read_FORF(IO *inout,PARAM *param,ATOM atom[],CONSTRAINT **constList,BOND **bond,
                   ANGLE **angle,DIHE **dihe,DIHE **impr,BOND **ub,real *eps,real *sig,
                   real *eps14,real *sig14,real *mass,real *q,int *frozen,int *nAtConst);

    void write_CONF(IO *inout,PARAM *param,ATOM atom[],real *x,real *y, real *z);

    void write_prop(IO *inout,PARAM *param,ENERGY *ener,PBC *box);

    void write_rest(IO *inout,PARAM *param,ENERGY *ener,BATH *bath,ATOM atom[],
                    real *x,real *y, real *z,real *vx,real *vy,real *vz,
                    real *fx,real *fy, real *fz);

    void write_DCD_header(IO *inout,CTRL *ctrl,PARAM *param, PBC *box,int frozen[]);

    void write_DCD_traj(IO *inout,PARAM *param,PBC *box,ATOM *atom,real x[],real y[],real z[],int frozen[]);

#ifdef	__cplusplus
}
#endif

#endif
