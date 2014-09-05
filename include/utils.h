#ifndef UTILSH_INCLUDED
#define UTILSH_INCLUDED

#ifdef	__cplusplus
extern "C" {
#endif

    real dist(const PBC *box, real delta[3]);

    real dist2(const PBC *box, real *dx, real *dy,real *dz);

    void freeze_atoms(const PARAM *param, real *vx, real *vy,real *vz,
                      real *fx,real *fy,real *fz,const int *frozen);

    void image_update(const PARALLEL *parallel, const PBC *box,real x[],real y[],real z[]);

    void image_array(const PBC *box, real dx[],real dy[],real dz[],const int size_array);

    void traj_rebuild(const PARAM *param,const PBC *box,const ATOM *atom,real x[],real y[], real z[]);
    
    void scale_box(PBC *box, const real cell0[9], const real scale);

    void vv_scale_box(PBC *box, const real scale);

    void box_to_lattice(const PBC *box,real lattice[6]);

    void box_to_crystal(const PBC *box,real crystal[6]);

    real getKin(const PARALLEL *parallel, const real vx[],const real vy[],
                   const real vz[],const real mass[],real dBuffer[]);

    void getKinStress(const PARALLEL *parallel,const real vx[],const real vy[],
                        const real vz[],const real mass[],real stress[6],
                        real dBuffer[]);

    void getKin0(PARAM *param, const PBC *box);

    void getDegFree(PARAM *param, const PBC *box);
    
    void getCom(const PARALLEL *parallel,const real mass[],
	    const real x[],const real y[],const real z[],
	    real com[3],real dBuffer[]);
    
    void getVom(const PARALLEL *parallel,const real mass[],
	    const real vx[],const real vy[],const real vz[],
	    real vom[3],real dBuffer[]);

    void nocase(char *str);

    int nint(const real x);

#ifdef	__cplusplus
}
#endif

#endif
