#ifndef INTEGRATEH_INCLUDED
#define INTEGRATEH_INCLUDED

#ifdef	__cplusplus
extern "C" {
#endif

// called once at the beginning of simulation for allocating working arrays used by integration
void integrators_allocate_arrays(CTRL *ctrl, PARAM *param);

// called once at the end of simulation for freeing working arrays used by integration
void integrators_free_arrays(CTRL *ctrl, PARAM *param);

void lf_integrate(CTRL *ctrl,PARAM *param,ENERGY *ener,PBC *box,
		  BATH *bath,CONSTRAINT constList[],PARALLEL *parallel,
		  double *x,double *y,double *z,
		  double *vx,double *vy,double *vz,
		  double *fx,double *fy,double *fz,
		  double *mass,double *rmass,int *nAtConst);

void lf_nve(PARAM *param,ENERGY *ener,PBC *box,CONSTRAINT constList[],PARALLEL *parallel,
	    double *x,double *y,double *z,
	    double *vx,double *vy,double *vz,
	    double *fx,double *fy,double *fz,
	    double *mass,double *rmass,int *nAtConst);

void lf_nvt_b(PARAM *param,ENERGY *ener,PBC *box,BATH *bath,CONSTRAINT constList[],PARALLEL *parallel,
	      double *x,double *y,double *z,
	      double *vx,double *vy,double *vz,
	      double *fx,double *fy,double *fz,
	      double *mass,double *rmass,int *nAtConst);

void lf_npt_b(PARAM *param,ENERGY *ener,PBC *box,BATH *bath,CONSTRAINT constList[],PARALLEL *parallel,
	      double *x,double *y,double *z,
	      double *vx,double *vy,double *vz,
	      double *fx,double *fy,double *fz,
	      double *mass,double *rmass,int *nAtConst);

void lf_nvt_h(PARAM *param,ENERGY *ener,PBC *box,BATH *bath,CONSTRAINT constList[],PARALLEL *parallel,
	      double *x,double *y,double *z,
	      double *vx,double *vy,double *vz,
	      double *fx,double *fy,double *fz,
	      double *mass,double *rmass,int *nAtConst);

void lf_npt_h(PARAM *param,ENERGY *ener,PBC *box,BATH *bath,CONSTRAINT constList[],PARALLEL *parallel,
	      double *x,double *y,double *z,
	      double *vx,double *vy,double *vz,
	      double *fx,double *fy,double *fz,
	      double *mass,double *rmass,int *nAtConst);


void vv_integrate(CTRL *ctrl,PARAM *param,ENERGY *ener,PBC *box,BATH *bath,
		  CONSTRAINT constList[],PARALLEL *parallel,double *x,double *y,
		  double *z,double *vx,double *vy,double *vz,
		  double *fx,double *fy,double *fz,
		  double *mass,double *rmass,int *nAtConst,int stage);

void vv_nve(PARAM *param,ENERGY *ener,PBC *box,CONSTRAINT constList[],PARALLEL *parallel,
	    double *x,double *y,double *z,double *vx,double *vy,double *vz,
	    double *fx,double *fy,double *fz,
	    double *mass,double *rmass,int *nAtConst,int stage);

void vv_nvt_b(PARAM *param,ENERGY *ener,PBC *box,BATH *bath,CONSTRAINT constList[],PARALLEL *parallel,
	      double *x,double *y,double *z,double *vx,double *vy,double *vz,
	      double *fx,double *fy,double *fz,
	      double *mass,double *rmass,int *nAtConst,int stage);

void vv_npt_b(PARAM *param,ENERGY *ener,PBC *box,BATH *bath,CONSTRAINT constList[],PARALLEL *parallel,
	      double *x,double *y,double *z,double *vx,double *vy,double *vz,
	      double *fx,double *fy,double *fz,
	      double *mass,double *rmass,int *nAtConst,int stage);

void vv_nvt_h(PARAM *param,ENERGY *ener,PBC *box,BATH *bath,CONSTRAINT constList[],PARALLEL *parallel,
	      double *x,double *y,double *z,double *vx,double *vy,double *vz,
	      double *fx,double *fy,double *fz,
	      double *mass,double *rmass,int *nAtConst,int stage);

void vv_npt_h(PARAM *param,ENERGY *ener,PBC *box,BATH *bath,CONSTRAINT constList[],PARALLEL *parallel,
	      double *x,double *y,double *z,double *vx,double *vy,double *vz,
	      double *fx,double *fy,double *fz,
	      double *mass,double *rmass,int *nAtConst,int stage);

#ifdef	__cplusplus
}
#endif

#endif