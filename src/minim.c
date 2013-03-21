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

