#ifndef IOH_INCLUDED
#define IOH_INCLUDED

void read_SIMU(SIMULPARAMS *simulCond,FORCEFIELD *ff,PBC *box);
void read_TOP(INPUTS *inp);
void read_PSF(INPUTS *inp,ATOM **atom,FORCEFIELD *ff,SIMULPARAMS *simulCond,CONSTRAINT **constList);
void read_PAR(INPUTS *inp);
void read_CONF(ATOM atom[],SIMULPARAMS *simulCond);
void setup(INPUTS *inp,ATOM atom[],FORCEFIELD *ff,SIMULPARAMS *simulCond,CONSTRAINT *constList);
void write_CONF(ATOM atom[],SIMULPARAMS *simulCond);
void write_prop(SIMULPARAMS *simulCond,ENERGY *ener,PBC *box);
void write_rest(SIMULPARAMS *simulCond,ENERGY *ener,ATOM *atom);
void read_rest(SIMULPARAMS *simulCond,ENERGY *ener,ATOM *atom);
void write_FORF(INPUTS *inp,ATOM atom[],FORCEFIELD *ff,SIMULPARAMS *simulCond);
void write_DCD_header(SIMULPARAMS *simulCond, PBC *box);
void write_DCD_traj(ATOM atom[], SIMULPARAMS *simulCond, PBC *box);
void free_temp_array(INPUTS *inp);
void error(int errorNumber);

/** Pointer to the output file. **/
extern FILE *outFile;

#endif
