#include <stdlib.h>
#include <stdio.h>

#include "global.h"
#include "user.h"


// A Dummy energy function
double MyEnergyFunction (const PARAM *param, const PBC *box,
                const double x[], const double y[], const double z[],
                double fx[], double fy[], double fz[],
                const int neighList[], const int neighPair[], const int neighOrder[],
                const int neighList14[])
{
    puts("Hello from my Plugin !");
    
    return 0.0;
}

