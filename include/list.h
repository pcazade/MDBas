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

#ifndef LISTH_INCLUDED
#define LISTH_INCLUDED

#ifdef	__cplusplus
extern "C" {
#endif

void makelist(CTRL *ctrl,PARAM *param,PBC *box,NEIGH *neigh,CONSTRAINT constList[],
	      BOND bond[],ANGLE angle[],DIHE dihe[],DIHE impr[],double x[], double y[],
	      double z[],int frozen[],int **neighList,int **neighPair,int **neighOrder,
	      int **neighList14,int ***exclList,int **exclPair);

void exclude_list(CTRL *ctrl,PARAM *param,NEIGH *neigh,CONSTRAINT constList[],
		  BOND bond[],ANGLE angle[],DIHE dihe[],DIHE impr[],
		  int **neighList14,int ***exclList,int **exclPair);

void verlet_list(PARAM *param,PBC *box,NEIGH *neigh,double x[],double y[],double z[],
		 int frozen[],int **neighList,int **neighPair,int **neighOrder,
		 int **exclList,int exclPair[]);

void verlet_list_update(PARAM *param,PBC *box,NEIGH *neigh,double x[],double y[],double z[],
			int frozen[],int **neighList,int neighPair[],int **exclList,
			int exclPair[]);

void link_cell_exclude_list(CTRL *ctrl,PARAM *param,NEIGH *neigh,CONSTRAINT constList[],
			    BOND bond[],ANGLE angle[],DIHE dihe[],DIHE impr[],
			    int **neighList14,int ***exclList,int **exclPair);

void link_cell_verlet_list(PARAM *param,PBC *box,NEIGH *neigh,double x[],double y[],double z[],
			   int frozen[],int **neighList,int **neighPair,int **neighOrder,
			   int **exclList,int exclPair[]);

void link_cell_verlet_list_update(PARAM *param,PBC *box,NEIGH *neigh,double x[],double y[],double z[],
				  int frozen[],int **neighList,int neighPair[],int neighOrder[],
				  int **exclList,int exclPair[]);

void fast_verlet_list(PARAM *param,PBC *box,NEIGH *neigh,double x[],double y[],double z[],
		      int frozen[],int **neighList,int **neighPair,int **neighOrder,
		      int **exclList,int exclPair[]);

#ifdef	__cplusplus
}
#endif

#endif