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

#ifndef MINIMH_INCLUDED
#define MINIMH_INCLUDED

#ifdef	__cplusplus
extern "C" {
#endif
    
void minimise(CTRL *ctrl,PARAM *param,ENERGY *ener,PBC *box,NEIGH *neigh,
	      ATOM atom[],BOND bond[],BOND ub[],ANGLE angle[],DIHE dihe[],
	      DIHE impr[]);

void steepestDescent(CTRL *ctrl,PARAM *param,ENERGY *ener,PBC *box,NEIGH *neigh,
		     ATOM atom[],BOND bond[],BOND ub[],ANGLE angle[],DIHE dihe[],
		     DIHE impr[],double x[],double y[],double z[],double fx[],
		     double fy[],double fz[]);

void conjugateGradients(CTRL *ctrl,PARAM *param,ENERGY *ener,PBC *box,NEIGH *neigh,
			ATOM atom[],BOND bond[],BOND ub[],ANGLE angle[],DIHE dihe[],
			DIHE impr[]);

#ifdef	__cplusplus
}
#endif

#endif