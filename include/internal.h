#ifndef INTERNALH_INCLUDED
#define INTERNALH_INCLUDED

#ifdef	__cplusplus
extern "C" {
#endif

    void bond_energy(const PARALLEL *parallel,ENERGY *ener,const PBC *box,
                     const BOND bond[],const double *x,
                     const double *y,const double *z,double *fx,double *fy,double *fz);

    void ub_energy(const PARALLEL *parallel,ENERGY *ener,const PBC *box,
                   const BOND ub[],const double *x,
                   const double *y,const double *z,double *fx,double *fy,double *fz);

    void angle_energy(const PARALLEL *parallel,ENERGY *ener,const PBC *box,
                      const ANGLE angle[],const double *x,
                      const double *y,const double *z,double *fx,double *fy,double *fz);

    void dihedral_energy(const PARALLEL *parallel,ENERGY *ener,
                         const PBC *box,const DIHE dihe[],const double *x,
                         const double *y,const double *z,double *fx,double *fy,double *fz);

    void improper_energy(const PARALLEL *parallel,ENERGY *ener,
                         const PBC *box,const DIHE impr[],const double *x,
                         const double *y,const double *z,double *fx,double *fy,double *fz);

#ifdef	__cplusplus
}
#endif

#endif
