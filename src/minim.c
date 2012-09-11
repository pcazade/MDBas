#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "global.h"
#include "energy.h"

void minimise(ATOM atom[],FORCEFIELD *ff,ENERGY *ener,SIMULPARAMS *simulCond)
{
  
}

void steepestDescent(ATOM atom[],FORCEFIELD *ff,ENERGY *ener,SIMULPARAMS *simulCond,PBC *box)
{
  double step = 1.0e-7 ;
  double prec = 1.0e-3 ;
  int maxSteps = 10000 ;
  
  double diff;
  double eprev = 0. , enow = 0. ;
  
  int i, currSt=0 ;
  
  energy(atom, ff, ener, simulCond, box);
  eprev = ener->pot;
  
//   printf("eprev : %lf\n",eprev);

  do
  {
    for (i=0 ; i < simulCond->natom ; i++)
    {
      atom[i].x += step * atom[i].fx ;
      atom[i].y += step * atom[i].fy ;
      atom[i].z += step * atom[i].fz ;
    }
    
    energy(atom, ff, ener, simulCond, box);
    enow = ener->pot;
//     printf("enow : %lf\n",enow);
    
    diff = fabs(enow-eprev);
    eprev=enow;
    
    currSt++;
    
//     printf("Steepest Descent : after step %d : EDiff = %lf \n", currSt, diff );
    
  } while ( (diff >= prec) && (currSt<=maxSteps) ) ;
  
}

void conjugateGradients(ATOM atom[],FORCEFIELD *ff,ENERGY *ener,SIMULPARAMS *simulCond)
{
  
}

