#include <stdlib.h>
#include <stdio.h>

#include "global.h"
#include "user.h"


// A Dummy energy function
real MyEnergyFunction (const PARAM *param, const PBC *box,
                         const real x[], const real y[], const real z[],
                         real fx[], real fy[], real fz[],
                         const int neighList[], const int neighPair[], const int neighOrder[],
                         const int neighList14[])
{
    puts("Hello from my Plugin !");

    return 0.0;
}

