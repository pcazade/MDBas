#ifndef POLARH_INCLUDED
#define POLARH_INCLUDED

#ifdef	__cplusplus
extern "C" {
#endif

void init_polar(CTRL *ctrl,PARAM *param,POLAR *polar,const real alPol[]);

void free_polar(CTRL *ctrl);

void static_field(PARAM *param,PARALLEL *parallel,PBC *box,
		  const real x[],const real y[],const real z[],
		  const real q[],int **neighList,const int neighPair[],
		  real dBuffer[]);

real polar_ener_iter(PARAM *param,PARALLEL *parallel,PBC *box, POLAR *polar,
		       const real x[],const real y[],const real z[],real fx[],
		       real fy[],real fz[],const real q[],const real alPol[],
		       real *virpol,int **neighList,const int neighPair[],real dBuffer[]);

real polar_ener_inv(PARAM *param,PARALLEL *parallel,PBC *box, POLAR *polar,
		      const real x[],const real y[],const real z[],real fx[],
		      real fy[],real fz[],const real q[],const real alPol[],
		      real *virpol,int **neighList,const int neighPair[],real dBuffer[]);

void polar_forces(PARAM *param,PARALLEL *parallel,PBC *box,const real x[],
		  const real y[],const real z[],real fx[],real fy[],
		  real fz[],const real q[],real *virpol,int **neighList,
		  const int neighPair[]);
  
#ifdef	__cplusplus
}
#endif

#endif