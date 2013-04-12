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

#ifndef IOH_INCLUDED
#define IOH_INCLUDED

#ifdef	__cplusplus
extern "C" {
#endif

void read_SIMU(IO *inout,CTRL *ctrl,PARAM *param,BATH *bath,NEIGH *neigh,EWALD *ewald,PBC *box);

void read_CONF(IO *inout,PARAM *param,ATOM **atom,double **x,double **y, double **z);

void read_rest(IO *inout,PARAM *param,ENERGY *ener,BATH *bath,ATOM **atom,
	       double **x,double **y, double **z,double **vx,double **vy,
	       double **vz,double **fx,double **fy, double **fz);

void read_FORF(IO *inout,PARAM *param,ATOM atom[],CONSTRAINT **constList,BOND **bond,
	       ANGLE **angle,DIHE **dihe,DIHE **impr,BOND **ub,double *eps,double *sig,
	       double *eps14,double *sig14,double *mass,double *q,int *frozen,int *nAtConst);

void write_CONF(IO *inout,PARAM *param,ATOM atom[],double *x,double *y, double *z);

void write_prop(IO *inout,PARAM *param,ENERGY *ener,PBC *box);

void write_rest(IO *inout,PARAM *param,ENERGY *ener,BATH *bath,ATOM atom[],
	        double *x,double *y, double *z,double *vx,double *vy,double *vz,
	        double *fx,double *fy, double *fz);

void write_DCD_header(IO *inout,CTRL *ctrl,PARAM *param, PBC *box,int frozen[]);

void write_DCD_traj(IO *inout,PARAM *param,PBC *box,double x[],double y[],double z[],int frozen[]);

#ifdef	__cplusplus
}
#endif

#endif
