#ifndef SPMEH_INCLUDED
#define SPMEH_INCLUDED

#ifdef	__cplusplus
extern "C" {
#endif

void init_spme(CTRL *ctrl,PARAM *param,EWALD *ewald,PBC *box);

void spme_free(PARAM *param);

void epl_cplx(EWALD *ewald);

void bspcoef(EWALD *ewald);

void bspgen(PARAM *param,EWALD *ewald);

// double spme_energy(PARAM *param,EWALD *ewald,PBC *box,const double x[],const double y[],const double z[],
// 	           double fx[],double fy[],double fz[],const double q[],double stress[6],double *virEwaldRec);

// double spme_energy(PARAM *param,EWALD *ewald,PBC *box,double x[],double y[],double z[],
// 	           double fx[],double fy[],double fz[],double q[],double stress[6],double *virEwaldRec);

double spme_energy(PARAM *param,EWALD *ewald,PBC *box,const double x[],const double y[],const double z[],
	           double fx[],double fy[],double fz[],const double q[],double stress[6],double *virEwaldRec);
#ifdef	__cplusplus
}
#endif

#endif