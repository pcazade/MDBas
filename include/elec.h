#ifndef ELECH_INCLUDED
#define ELECH_INCLUDED

#ifdef	__cplusplus
extern "C" {
#endif

    real coulomb_none(const PARAM *param,real *delec,const real qel,
                        const real r2,const real rt);

    void coulomb_full(ENERGY *ener,PARAM *param,PARALLEL *parallel,PBC *box,real *x,real *y,
                      real *z,real *fx, real *fy, real *fz,real *q,
                      int **exclList,int *exclPair);

    real coulomb_shift1(const PARAM *param,real *delec,const real qel,
                          const real r2,const real rt);

    real coulomb_shift2(const PARAM *param,real *delec,const real qel,
                          const real r2,const real rt);

    real coulomb_switch(const PARAM *param,real *delec,const real qel,
                          const real r2,const real rt);
    
    real coulomb_damp(const PARAM *param,real *delec,const real qel,
			const real r2,const real rt);

    real coulomb14_none(const PARAM *param,real *delec,const real qel,
                          const real r2,const real rt);

    real coulomb14_full(const PARAM *param,real *delec,const real qel,
                          const real r2,const real rt);

    real coulomb14_shift1(const PARAM *param,real *delec,const real qel,
                            const real r2,const real rt);

    real coulomb14_shift2(const PARAM *param,real *delec,const real qel,
                            const real r2,const real rt);

    real coulomb14_switch(const PARAM *param,real *delec,const real qel,
                            const real r2,const real rt);
    
    real coulomb14_damp(const PARAM *param,real *delec,const real qel,
			  const real r2,const real rt);

    /** Pointer to the output file. **/
    extern FILE *outFile;

#ifdef	__cplusplus
}
#endif

#endif
