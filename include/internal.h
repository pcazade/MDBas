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

#ifndef INTERNALH_INCLUDED
#define INTERNALH_INCLUDED

#ifdef	__cplusplus
extern "C" {
#endif
    
void bond_energy(const PARAM *param,ENERGY *ener,const PBC *box,const BOND bond[],const double *x,
		 const double *y,const double *z,double *fx,double *fy,double *fz);

void ub_energy(const PARAM *param,ENERGY *ener,const PBC *box,const BOND ub[],const double *x,
	       const double *y,const double *z,double *fx,double *fy,double *fz);

void angle_energy(const PARAM *param,ENERGY *ener,const PBC *box,const ANGLE angle[],const double *x,
		  const double *y,const double *z,double *fx,double *fy,double *fz);

void dihedral_energy(const PARAM *param,ENERGY *ener,const PBC *box,const DIHE dihe[],const double *x,
		     const double *y,const double *z,double *fx,double *fy,double *fz);

void improper_energy(const PARAM *param,ENERGY *ener,const PBC *box,const DIHE impr[],const double *x,
		     const double *y,const double *z,double *fx,double *fy,double *fz);

#ifdef	__cplusplus
}
#endif

#endif