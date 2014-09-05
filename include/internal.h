#ifndef INTERNALH_INCLUDED
#define INTERNALH_INCLUDED

#ifdef	__cplusplus
extern "C" {
#endif

    void bond_energy(const PARALLEL *parallel,ENERGY *ener,const PBC *box,
                     const BOND bond[],const real *x,
                     const real *y,const real *z,real *fx,real *fy,real *fz);

    void ub_energy(const PARALLEL *parallel,ENERGY *ener,const PBC *box,
                   const BOND ub[],const real *x,
                   const real *y,const real *z,real *fx,real *fy,real *fz);

    void angle_energy(const PARALLEL *parallel,ENERGY *ener,const PBC *box,
                      const ANGLE angle[],const real *x,
                      const real *y,const real *z,real *fx,real *fy,real *fz);

    void dihedral_energy(const PARALLEL *parallel,ENERGY *ener,
                         const PBC *box,const DIHE dihe[],const real *x,
                         const real *y,const real *z,real *fx,real *fy,real *fz);

    void improper_energy(const PARALLEL *parallel,ENERGY *ener,
                         const PBC *box,const DIHE impr[],const real *x,
                         const real *y,const real *z,real *fx,real *fy,real *fz);

#ifdef	__cplusplus
}
#endif

#endif
