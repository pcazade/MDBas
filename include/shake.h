#ifndef SHAKEH_INCLUDED
#define SHAKEH_INCLUDED

void lf_shake(ATOM *atom,SIMULPARAMS *simulCond,CONSTRAINT *constList,DELTA *dd,PBC *box);
void vv_shake_r(ATOM *atom,SIMULPARAMS *simulCond,CONSTRAINT *constList,DELTA *dd,PBC *box);
void vv_shake_v(ATOM *atom,SIMULPARAMS *simulCond,CONSTRAINT *constList,DELTA *dd);

#endif