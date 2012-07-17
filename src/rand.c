#include <math.h>
#include "global.h"
#include "dSFMT.h"

static dsfmt_t dsfmt;

void init_rand(unsigned int seed)
{
  dsfmt_init_gen_rand(&dsfmt,seed);
}

double get_rand()
{
  return dsfmt_genrand_open_open(&dsfmt);
}

void get_BoxMuller(double *u, double *v)
{
  double a,b,s;
  
  do
  {
    
      a = 2. * get_rand() - 1.;
      b = 2. * get_rand() - 1.;
      s = a*a + b*b;
      
  }while (s >= 1.);
  
  *u= a*sqrt(-2.*log(s)/s);
  *v= b*sqrt(-2.*log(s)/s);

}
