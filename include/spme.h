#ifndef SPMEH_INCLUDED
#define SPMEH_INCLUDED

#ifdef	__cplusplus
extern "C" {
#endif

void init_spme(CTRL *ctrl,PARAM *param,PARALLEL *parallel,EWALD *ewald,PBC *box);

void spme_free(PARALLEL *parallel);

void epl_cplx(EWALD *ewald);

void bspcoef(EWALD *ewald);

void bspgen(PARALLEL *parallel,EWALD *ewald);

// double spme_energy(PARAM *param,EWALD *ewald,PBC *box,const double x[],const double y[],const double z[],
// 	           double fx[],double fy[],double fz[],const double q[],double stress[6],double *virEwaldRec);

// double spme_energy(PARAM *param,EWALD *ewald,PBC *box,double x[],double y[],double z[],
// 	           double fx[],double fy[],double fz[],double q[],double stress[6],double *virEwaldRec);

double spme_energy(PARAM *param,PARALLEL *parallel,EWALD *ewald,PBC *box,const double x[],
		   const double y[],const double z[],double fx[],double fy[],double fz[],
		   const double q[],double stress[6],double *virEwaldRec,double dBuffer[]);
#ifdef	__cplusplus
}
#endif

#endif