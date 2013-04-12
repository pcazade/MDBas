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
