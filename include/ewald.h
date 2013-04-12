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

#ifndef EWALDH_INCLUDED
#define EWALDH_INCLUDED

#ifdef	__cplusplus
extern "C" {
#endif

void init_ewald(CTRL *ctrl,PARAM *param,EWALD *ewald,PBC *box);

void ewald_free(EWALD *ewald);

double ewald_rec(PARAM *param,EWALD *ewald,PBC *box,const double x[],const double y[],const double z[],
	       double fx[],double fy[],double fz[],const double q[],double stress[6],double *virEwaldRec);

double ewald_dir(EWALD *ewald,double *dEwaldDir,const double qel,
		 const double r,const double rt);

double ewald_corr(EWALD *ewald,double *dEwaldCorr,const double qel,
	         const double r,const double rt);

double ewald_dir14(PARAM *param,EWALD *ewald,double *dEwaldDir,const double qel,
	         const double r,const double rt);

double ewald_corr14(PARAM *param,EWALD *ewald,double *dEwaldCorr,
		    const double qel,const double r,const double rt);

#ifdef	__cplusplus
}
#endif

#endif