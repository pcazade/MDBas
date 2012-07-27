#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "global.h"
#include "energy.h"

void minimise(ATOM *atom,FORCEFIELD *ff,ENERGYFORCE *enerFor,SIMULPARAMS *simulCond)
{
  
}

void steepestDescent(ATOM *atom,FORCEFIELD *ff,ENERGYFORCE *enerFor,SIMULPARAMS *simulCond)
{
  double step = 1.0e-7 ;
  double prec = 1.0e-3 ;
  int maxSteps = 10000 ;
  
  double diff;
  double eprev = 0. , enow = 0. ;
  
  int i, currSt=0 ;
  
  energy(atom, ff, enerFor, simulCond);
  eprev = enerFor->energyPot;
  
//   printf("eprev : %lf\n",eprev);

  do
  {
    for (i=0 ; i < atom->natom ; i++)
    {
      atom->x[i] += step * atom->fx[i] ;
      atom->y[i] += step * atom->fy[i] ;
      atom->z[i] += step * atom->fz[i] ;
    }
    
    energy(atom, ff, enerFor, simulCond);
    enow = enerFor->energyPot;
//     printf("enow : %lf\n",enow);
    
    diff = fabs(enow-eprev);
    eprev=enow;
    
    currSt++;
    
//     printf("Steepest Descent : after step %d : EDiff = %lf \n", currSt, diff );
    
  } while ( (diff >= prec) && (currSt<=maxSteps) ) ;
  
}

void conjugateGradients(ATOM *atom,FORCEFIELD *ff,ENERGYFORCE *enerFor,SIMULPARAMS *simulCond)
{
  
}

