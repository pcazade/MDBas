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

#ifndef INITH_INCLUDED
#define INITH_INCLUDED

#ifdef	__cplusplus
extern "C" {
#endif

void init_system(int argc, char* argv[],IO *inout,CTRL *ctrl,PARAM *param,ENERGY *ener,BATH *bath,
		 NEIGH *neigh,EWALD *ewald,PBC *box,ATOM **atom,CONSTRAINT **constList,
		 BOND **bond,ANGLE **angle,DIHE **dihe,DIHE **impr,BOND **ub,double **x,
		 double **y, double **z,double **vx,double **vy,double **vz,double **fx,
		 double **fy, double **fz,double **mass,double **rmass,double **q,
		 double **eps,double **sig,double **eps14,double **sig14,int **frozen,
		 int **nAtConst,int **neighList,int **neighPair,int **neighOrder,
		 int **neighList14,int ***exclList,int **exclPair);

void setup(CTRL *ctrl,PARAM *param,ATOM atom[],CONSTRAINT **constList,
	   BOND **bond,ANGLE **angle,DIHE **dihe,DIHE **impr,BOND **ub,
	   double mass[],double rmass[],int frozen[],int nAtConst[]);

void init_vel(PARAM *param,PBC *box,CONSTRAINT constList[],double x[],
	      double y[],double z[],double vx[],double vy[],double vz[],
	      double mass[],double rmass[],int frozen[],int nAtConst[]);

void init_constvel(PARAM *param,PBC *box,CONSTRAINT constList[],double x[],
		   double y[],double z[],double vx[],double vy[],double vz[],
		   double mass[],int nAtConst[]);

void init_box(PBC *box);

void free_all(CTRL *ctrl,PARAM *param,EWALD *ewald,ATOM **atom,CONSTRAINT **constList,
	      BOND **bond,ANGLE **angle,DIHE **dihe,DIHE **impr,BOND **ub,double **x,
	      double **y, double **z,double **vx,double **vy,double **vz,double **fx,
	      double **fy, double **fz,double **mass,double **rmass,double **q,
	      double **eps,double **sig,double **eps14,double **sig14,int **frozen,
	      int **nAtConst,int **neighList,int **neighPair,int **neighOrder,
	      int **neighList14,int ***exclList,int **exclPair);

#ifdef	__cplusplus
}
#endif

#endif