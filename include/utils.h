#ifndef UTILSH_INCLUDED
#define UTILSH_INCLUDED

double distance(int i,int j, ATOM atom[],double *delta,SIMULPARAMS *simulCond,PBC *box);
double distance2(int i,int j, ATOM atom[],DELTA *d,SIMULPARAMS *simulCond,PBC *box);
void init_vel(ATOM atom[],SIMULPARAMS *simulCond,CONSTRAINT *constList,PBC *box);
void init_constvel(ATOM atom[],SIMULPARAMS *simulCond,CONSTRAINT *constList,PBC *box);
void image_update(ATOM atom[],SIMULPARAMS *simulCond,PBC *box);
void image_array(int size_array,DELTA *d,SIMULPARAMS *simulCond,PBC *box);
void init_box(PBC *box);
<<<<<<< .mine
double kinetic(ATOM *atom,SIMULPARAMS *simulCond);
double stress_kinetic(ATOM *atom,SIMULPARAMS *simulCond,double *stress);
void get_kinfromtemp(ATOM *atom,SIMULPARAMS *simulCond,PBC *box);
void get_degfree(ATOM *atom,SIMULPARAMS *simulCond,PBC *box);
=======
double kinetic(ATOM atom[],SIMULPARAMS *simulCond);
void get_kinfromtemp(ATOM atom[],SIMULPARAMS *simulCond,PBC *box);
void get_degfree(ATOM atom[],SIMULPARAMS *simulCond,PBC *box);
>>>>>>> .r66
void nocase(char *str);
int nint(double x);

#endif
