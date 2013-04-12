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

#ifndef USER_H
#define	USER_H

#ifdef	__cplusplus
extern "C" {
#endif

/* pointers for user-defined functions */
typedef double (*UserEnergyPtr) (const PARAM *param, const PBC *box,
                const double x[], const double y[], const double z[],
                double fx[], double fy[], double fz[],
                const int neighList[], const int neighPair[], const int neighOrder[],
                const int neighList14[]);

/* Functions from user.c ; for managing user extensions */
UserEnergyPtr loadUserPlugin(const char pluginName[], const char funcName[]);
void closeUserPlugin(void);

#ifdef	__cplusplus
}
#endif

#endif	/* USER_H */

