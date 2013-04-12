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

#ifndef UTILSH_INCLUDED
#define UTILSH_INCLUDED

#ifdef	__cplusplus
extern "C" {
#endif

double dist(const PBC *box, double delta[3]);

double dist2(const PBC *box, double *dx, double *dy,double *dz);

void freeze_atoms(const PARAM *param, double *vx, double *vy,double *vz,
		  double *fx,double *fy,double *fz,const int *frozen);

void image_update(const PARAM *param, const PBC *box,double x[],double y[],double z[]);

void image_array(const PBC *box, double dx[],double dy[],double dz[],const int size_array);

void scale_box(PBC *box, const double cell0[9], const double scale);

void vv_scale_box(PBC *box, const double scale);

void box_to_lattice(const PBC *box,double lattice[6]);

void box_to_crystal(const PBC *box,double crystal[6]);

double kinetic(const PARAM *param, const double vx[],const double vy[],
	       const double vz[],const double mass[]);

void stress_kinetic(const PARAM *param,const double vx[],const double vy[],
		    const double vz[],const double mass[],double stress[6]);

void get_kinfromtemp(PARAM *param, const PBC *box);

void get_degfree(PARAM *param, const PBC *box);

void nocase(char *str);

int nint(const double x);

#ifdef	__cplusplus
}
#endif

#endif
