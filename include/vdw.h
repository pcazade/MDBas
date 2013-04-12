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

#ifndef VDWH_INCLUDED
#define VDWH_INCLUDED

#ifdef	__cplusplus
extern "C" {
#endif

double vdw_none(const PARAM *param,double *dvdw,const double veps,
		const double vsig,const double r2, const double rt);

void vdw_full(const PARAM *param, ENERGY *ener, const PBC *box,double *x,
	      double *y,double *z,double *fx, double *fy, double *fz,double *eps,double *sig,
	      int **exclList,int *exclPair);

double vdw_switch(const PARAM *param,double *dvdw,const double veps,
		  const double vsig,const double r2, const double rt);

double vdw14_none(const PARAM *param,double *dvdw,const double veps,
		  const double vsig,const double r2, const double rt);

double vdw14_full(const PARAM *param,double *dvdw,const double veps,
		  const double vsig,const double r2, const double rt);

double vdw14_switch(const PARAM *param,double *dvdw,const double veps,
		    const double vsig,const double r2, const double rt);

#ifdef	__cplusplus
}
#endif

#endif