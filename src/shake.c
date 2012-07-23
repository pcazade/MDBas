#include <stdlib.h>
#include <math.h>

void lf_shake(ATOM *atom,SIMULPARAMS *simulCond,FORCEFIELD *ff,CONSTRAINT constList)
{
  int i,ia,ib,icycle,converged;
  double *xo,*yo,*zo,*xt,*yt,*zt,*rd,*rt,ts2,maxdist;
  DELTA *dd,*dt;
  CONSTRAINT constList;
  
  for(i=0;i<ff->nconst;i++)
  {
    ia=constList[i]->x;
    ib=constList[i]->y;
    
    rd[i]=distance2(ia,ib,atom,&(dd[i]),simulCond);
  }
  
  icycle=0;
  converged=0;
  
  ts2=X2(simulCond->timeStep);

  do
  {
    maxdist=0.;
    
    for(i=0;i<ff->nconst;i++)
    {
      ia=constList[i]->x;
      ib=constList[i]->y;
      
      rt[i]=distance2(ia,ib,atom,&(dt[i]),simulCond);
      
      
    }
    
    
    
    
  }while( (!converged) && (icycle<simulCond->maxcycle) );
}
