#ifndef UTILSH_INCLUDED
#define UTILSH_INCLUDED

#ifdef	__cplusplus
extern "C" {
#endif

    double dist(const PBC *box, double delta[3]);

    double dist2(const PBC *box, double *dx, double *dy,double *dz);

    void freeze_atoms(const PARAM *param, double *vx, double *vy,double *vz,
                      double *fx,double *fy,double *fz,const int *frozen);

    void image_update(const PARALLEL *parallel, const PBC *box,double x[],double y[],double z[]);

    void image_array(const PBC *box, double dx[],double dy[],double dz[],const int size_array);

    void scale_box(PBC *box, const double cell0[9], const double scale);

    void vv_scale_box(PBC *box, const double scale);

    void box_to_lattice(const PBC *box,double lattice[6]);

    void box_to_crystal(const PBC *box,double crystal[6]);

    double getKin(const PARALLEL *parallel, const double vx[],const double vy[],
                   const double vz[],const double mass[],double dBuffer[]);

    void getKinStress(const PARALLEL *parallel,const double vx[],const double vy[],
                        const double vz[],const double mass[],double stress[6],
                        double dBuffer[]);

    void getKin0(PARAM *param, const PBC *box);

    void getDegFree(PARAM *param, const PBC *box);
    
    void getCom(const PARALLEL *parallel,const double mass[],
	    const double x[],const double y[],const double z[],
	    double com[3],double dBuffer[]);
    
    void getVom(const PARALLEL *parallel,const double mass[],
	    const double vx[],const double vy[],const double vz[],
	    double vom[3],double dBuffer[]);

    void nocase(char *str);

    int nint(const double x);

#ifdef	__cplusplus
}
#endif

#endif
