#ifndef VDWH_INCLUDED
#define VDWH_INCLUDED

#ifdef	__cplusplus
extern "C" {
#endif

double vdw_none(const PARAM *param,double *dvdw,const double veps,
		const double vsig,const double r2, const double rt);

void vdw_full(const PARAM *param, ENERGY *ener, const PBC *box,double *x,
	      double *y,double *z,double *fx, double *fy, double *fz,double *eps,double *sig,
	      int **exclList,int *exclPair);

double vdw_switch(const PARAM *param,double *dvdw,const double veps,
		  const double vsig,const double r2, const double rt);

double vdw14_none(const PARAM *param,double *dvdw,const double veps,
		  const double vsig,const double r2, const double rt);

double vdw14_full(const PARAM *param,double *dvdw,const double veps,
		  const double vsig,const double r2, const double rt);

double vdw14_switch(const PARAM *param,double *dvdw,const double veps,
		    const double vsig,const double r2, const double rt);

#ifdef	__cplusplus
}
#endif

#endif