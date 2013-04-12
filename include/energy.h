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

#ifndef ENERGYH_INCLUDED
#define ENERGYH_INCLUDED

#ifdef	__cplusplus
extern "C" {
#endif

void init_energy_ptrs(CTRL *ctrl);

void energy(CTRL *ctrl,PARAM *param,ENERGY *ener,EWALD *ewald,PBC *box,NEIGH *neigh,
	    BOND bond[],BOND ub[],ANGLE angle[],DIHE dihe[],DIHE impr[],
	    const double x[],const double y[], const double z[],
	    double vx[],double vy[], double vz[],double fx[],double fy[],
	    double fz[],const double q[],const double eps[],const double sig[],
	    const double eps14[],const double sig14[],const int frozen[],
	    const int neighList[],const int neighPair[],const int neighOrder[],
	    const int neighList14[],int **exclList,const int exclPair[]);

void nonbond_energy(PARAM *param,ENERGY *ener,PBC *box,
		    const double x[],const double y[],const double z[],double fx[],double fy[],
		    double fz[],const double q[],const double eps[],const double sig[],
		    const int neighList[],const int neighPair[],const int neighOrder[]);

void nonbond14_energy(PARAM *param,ENERGY *ener,PBC *box,NEIGH *neigh,
		      const double x[],const double y[],const double z[],double fx[],double fy[],
		      double fz[],const double q[],const double eps[],const double sig[],
		      const int neighList14[]);

void ewald_energy(CTRL *ctrl,PARAM *param,ENERGY *ener,EWALD *ewald,PBC *box,const double x[],
		  const double y[],const double z[],double fx[],double fy[],
		  double fz[],const double q[],const double eps[],const double sig[],
		  const int neighList[],const int neighPair[],const int neighOrder[],
		  int **exclList,const int exclPair[]);

void ewald14_energy(PARAM *param,ENERGY *ener,EWALD *ewald,PBC *box,NEIGH *neigh,const double x[],
		    const double y[],const double z[],double fx[],double fy[],
		    double fz[],const double q[],const double eps[],const double sig[],
		    const int neighList14[]);

/* pointers to function for electrostatic energy functions */
double (*ptr_coulomb)(const PARAM *param,double *delec,const double qel,
		      const double r2,const double rt);

double (*ptr_coulomb14)(const PARAM *param,double *delec,const double qel,
			const double r2,const double rt);

/* pointers to function for vdw energy functions */
double (*ptr_vdw)(const PARAM *param,double *dvdw,const double veps,
		  const double vsig,const double r2, const double rt);

double (*ptr_vdw14)(const PARAM *param,double *dvdw,const double veps,
		    const double vsig,const double r2, const double rt);

#ifdef	__cplusplus
}
#endif

#endif