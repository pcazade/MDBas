#ifndef IOH_INCLUDED
#define IOH_INCLUDED

void read_TOP(INPUTS *inp);
void read_PSF(INPUTS *inp,ATOM *atom,FORCEFIELD *ff,ENERGYFORCE *enerFor,SIMULPARAMS *simulCond);
void read_PAR(INPUTS *inp);
void read_CONF(INPUTS *inp, ATOM *atom);
void setup(INPUTS *inp,ATOM *atom,FORCEFIELD *ff,ENERGYFORCE *enerFor,SIMULPARAMS *simulCond);
void free_temp_array(INPUTS *inp);
void error(int errorNumber);

#endif
