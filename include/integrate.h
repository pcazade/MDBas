/*
 * Copyright (c) 2013 Pierre-Andre Cazade
 * Copyright (c) 2013 Florent hedin
 * 
 * This file is part of MDBas.
 *
 * MDBas is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MDBas is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with MDBas.  If not, see <http://www.gnu.org/licenses/>.
 */

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
		  BATH *bath,CONSTRAINT constList[],
		  double *x,double *y,double *z,
		  double *vx,double *vy,double *vz,
		  double *fx,double *fy,double *fz,
		  double *mass,double *rmass,int *nAtConst);

void lf_nve(PARAM *param,ENERGY *ener,PBC *box,CONSTRAINT constList[],
	    double *x,double *y,double *z,
	    double *vx,double *vy,double *vz,
	    double *fx,double *fy,double *fz,
	    double *mass,double *rmass,int *nAtConst);

void lf_nvt_b(PARAM *param,ENERGY *ener,PBC *box,BATH *bath,CONSTRAINT constList[],
	      double *x,double *y,double *z,
	      double *vx,double *vy,double *vz,
	      double *fx,double *fy,double *fz,
	      double *mass,double *rmass,int *nAtConst);

void lf_npt_b(PARAM *param,ENERGY *ener,PBC *box,BATH *bath,CONSTRAINT constList[],
	      double *x,double *y,double *z,
	      double *vx,double *vy,double *vz,
	      double *fx,double *fy,double *fz,
	      double *mass,double *rmass,int *nAtConst);

void lf_nvt_h(PARAM *param,ENERGY *ener,PBC *box,BATH *bath,CONSTRAINT constList[],
	      double *x,double *y,double *z,
	      double *vx,double *vy,double *vz,
	      double *fx,double *fy,double *fz,
	      double *mass,double *rmass,int *nAtConst);

void lf_npt_h(PARAM *param,ENERGY *ener,PBC *box,BATH *bath,CONSTRAINT constList[],
	      double *x,double *y,double *z,
	      double *vx,double *vy,double *vz,
	      double *fx,double *fy,double *fz,
	      double *mass,double *rmass,int *nAtConst);


void vv_integrate(CTRL *ctrl,PARAM *param,ENERGY *ener,PBC *box,BATH *bath,
		  CONSTRAINT constList[],double *x,double *y,
		  double *z,double *vx,double *vy,double *vz,
		  double *fx,double *fy,double *fz,
		  double *mass,double *rmass,int *nAtConst,int stage);

void vv_nve(PARAM *param,ENERGY *ener,PBC *box,CONSTRAINT constList[],
	    double *x,double *y,double *z,double *vx,double *vy,double *vz,
	    double *fx,double *fy,double *fz,
	    double *mass,double *rmass,int *nAtConst,int stage);

void vv_nvt_b(PARAM *param,ENERGY *ener,PBC *box,BATH *bath,CONSTRAINT constList[],
	      double *x,double *y,double *z,double *vx,double *vy,double *vz,
	      double *fx,double *fy,double *fz,
	      double *mass,double *rmass,int *nAtConst,int stage);

void vv_npt_b(PARAM *param,ENERGY *ener,PBC *box,BATH *bath,CONSTRAINT constList[],
	      double *x,double *y,double *z,double *vx,double *vy,double *vz,
	      double *fx,double *fy,double *fz,
	      double *mass,double *rmass,int *nAtConst,int stage);

void vv_nvt_h(PARAM *param,ENERGY *ener,PBC *box,BATH *bath,CONSTRAINT constList[],
	      double *x,double *y,double *z,double *vx,double *vy,double *vz,
	      double *fx,double *fy,double *fz,
	      double *mass,double *rmass,int *nAtConst,int stage);

void vv_npt_h(PARAM *param,ENERGY *ener,PBC *box,BATH *bath,CONSTRAINT constList[],
	      double *x,double *y,double *z,double *vx,double *vy,double *vz,
	      double *fx,double *fy,double *fz,
	      double *mass,double *rmass,int *nAtConst,int stage);

#ifdef	__cplusplus
}
#endif

#endif