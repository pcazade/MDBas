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

/**
 * \file minim.c
 * \brief Contains functions performing energy minimisation.
 * \author Pierre-Andre Cazade and Florent Hedin
 * \version alpha-branch
 * \date 2012
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "global.h"
#include "energy.h"

/** Pointer to the output file. **/
extern FILE *outFile;

void minimise(CTRL *ctrl,PARAM *param,ENERGY *ener,PBC *box,NEIGH *neigh,
              ATOM atom[],BOND bond[],BOND ub[],ANGLE angle[],DIHE dihe[],
              DIHE impr[])
{

}

void steepestDescent(CTRL *ctrl,PARAM *param,ENERGY *ener,PBC *box,NEIGH *neigh,
                     ATOM atom[],BOND bond[],BOND ub[],ANGLE angle[],DIHE dihe[],
                     DIHE impr[],double x[],double y[],double z[],double fx[],
                     double fy[],double fz[])
{
    double step = 1.0e-7 ;
    double prec = 1.0e-3 ;
    int maxSteps = 10000 ;

    double diff;
    double eprev = 0. , enow = 0. ;

    int i, currSt=0 ;

    //energy(ctrl,param,ener,box,neigh,atom,bond,ub,angle,dihe,impr,vdw);
    eprev = ener->pot;

//   fprintf(outFile,"eprev : %lf\n",eprev);

    do
    {
        for (i=0 ; i < param->nAtom ; i++)
        {
            x[i] += step * fx[i] ;
            y[i] += step * fy[i] ;
            z[i] += step * fz[i] ;
        }

        //energy(ctrl,param,ener,box,neigh,atom,bond,ub,angle,dihe,impr,vdw);
        enow = ener->pot;
//     fprintf(outFile,"enow : %lf\n",enow);

        diff = fabs(enow-eprev);
        eprev=enow;

        currSt++;

//     fprintf(outFile,"Steepest Descent : after step %d : EDiff = %lf \n", currSt, diff );

    } while ( (diff >= prec) && (currSt<=maxSteps) ) ;

}

void conjugateGradients(CTRL *ctrl,PARAM *param,ENERGY *ener,PBC *box,NEIGH *neigh,
                        ATOM atom[],BOND bond[],BOND ub[],ANGLE angle[],DIHE dihe[],
                        DIHE impr[])
{

}

