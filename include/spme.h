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