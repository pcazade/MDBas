#ifndef EWALDH_INCLUDED
#define EWALDH_INCLUDED

#ifdef	__cplusplus
extern "C" {
#endif

void init_ewald(CTRL *ctrl,PARAM *param,PARALLEL *parallel,EWALD *ewald,PBC *box);

void ewald_free(EWALD *ewald);

double ewald_rec(PARAM *param,PARALLEL *parallel,EWALD *ewald,PBC *box,const double x[],
		 const double y[],const double z[],double fx[],double fy[],double fz[],
		 const double q[],double stress[6],double *virEwaldRec,double dBuffer[]);

double ewald_dir(EWALD *ewald,double *dEwaldDir,const double qel,
		 const double r,const double rt);

double ewald_corr(EWALD *ewald,double *dEwaldCorr,const double qel,
	         const double r,const double rt);

double ewald_dir14(PARAM *param,EWALD *ewald,double *dEwaldDir,const double qel,
	         const double r,const double rt);

double ewald_corr14(PARAM *param,EWALD *ewald,double *dEwaldCorr,
		    const double qel,const double r,const double rt);

#ifdef	__cplusplus
}
#endif

#endif