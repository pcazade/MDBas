#ifndef POLARH_INCLUDED
#define POLARH_INCLUDED

#ifdef	__cplusplus
extern "C" {
#endif

void init_polar(CTRL *ctrl,PARAM *param,POLAR *polar,const double alPol[]);

void free_polar(CTRL *ctrl);

void static_field(PARAM *param,PARALLEL *parallel,PBC *box,
		  const double x[],const double y[],const double z[],
		  const double q[],int **neighList,const int neighPair[],
		  double dBuffer[]);

double polar_ener_iter(PARAM *param,PARALLEL *parallel,PBC *box, POLAR *polar,
		       const double x[],const double y[],const double z[],double fx[],
		       double fy[],double fz[],const double q[],const double alPol[],
		       double *virpol,int **neighList,const int neighPair[],double dBuffer[]);

double polar_ener_inv(PARAM *param,PARALLEL *parallel,PBC *box, POLAR *polar,
		      const double x[],const double y[],const double z[],double fx[],
		      double fy[],double fz[],const double q[],const double alPol[],
		      double *virpol,int **neighList,const int neighPair[],double dBuffer[]);

void polar_forces(PARAM *param,PARALLEL *parallel,PBC *box,const double x[],
		  const double y[],const double z[],double fx[],double fy[],
		  double fz[],const double q[],double *virpol,int **neighList,
		  const int neighPair[]);
  
#ifdef	__cplusplus
}
#endif

#endif