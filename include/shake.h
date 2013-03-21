#ifndef SHAKEH_INCLUDED
#define SHAKEH_INCLUDED

#ifdef	__cplusplus
extern "C" {
#endif
    
void shake_allocate_arrays(const PARAM *param);

void shake_free_arrays();

void lf_shake(PARAM *param,PBC *box,CONSTRAINT constList[],
	      double x[],double y[],double z[],
	      double ddx[],double ddy[],double ddz[],double rmass[],
	      int *nAtConst,double stress[6],double *virshake);

void vv_shake_r(PARAM *param,PBC *box,CONSTRAINT constList[],
		double x[],double y[],double z[],
		double vx[],double vy[],double vz[],
		double ddx[],double ddy[],double ddz[],double rmass[],
		int *nAtConst,double stress[6],double *virshake);

void vv_shake_v(PARAM *param,CONSTRAINT constList[],
		double vx[],double vy[],double vz[],double ddx[],
		double ddy[],double ddz[],double rmass[],int *nAtConst);

#ifdef	__cplusplus
}
#endif

#endif
