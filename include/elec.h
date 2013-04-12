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

#ifndef ELECH_INCLUDED
#define ELECH_INCLUDED

#ifdef	__cplusplus
extern "C" {
#endif

double coulomb_none(const PARAM *param,double *delec,const double qel,
		      const double r2,const double rt);

void coulomb_full(ENERGY *ener,PARAM *param,PBC *box,double *x,double *y,
		  double *z,double *fx, double *fy, double *fz,double *q,
		  int **exclList,int *exclPair);

double coulomb_shift1(const PARAM *param,double *delec,const double qel,
		      const double r2,const double rt);

double coulomb_shift2(const PARAM *param,double *delec,const double qel,
		      const double r2,const double rt);

double coulomb_switch(const PARAM *param,double *delec,const double qel,
		      const double r2,const double rt);

double coulomb14_none(const PARAM *param,double *delec,const double qel,
		      const double r2,const double rt);

double coulomb14_full(const PARAM *param,double *delec,const double qel,
		      const double r2,const double rt);

double coulomb14_shift1(const PARAM *param,double *delec,const double qel,
			const double r2,const double rt);

double coulomb14_shift2(const PARAM *param,double *delec,const double qel,
			const double r2,const double rt);

double coulomb14_switch(const PARAM *param,double *delec,const double qel,
			const double r2,const double rt);

/** Pointer to the output file. **/
extern FILE *outFile;

#ifdef	__cplusplus
}
#endif

#endif