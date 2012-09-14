#ifndef UTILSH_INCLUDED
#define UTILSH_INCLUDED

double distance(int i,int j, ATOM atom[],double delta[3],SIMULPARAMS *simulCond,PBC *box);
double distance2(int i,int j, ATOM atom[],DELTA *d,SIMULPARAMS *simulCond,PBC *box);

void init_vel(ATOM atom[],SIMULPARAMS *simulCond,CONSTRAINT *constList,PBC *box);
void init_constvel(ATOM atom[],SIMULPARAMS *simulCond,CONSTRAINT *constList,PBC *box);

void image_update(ATOM atom[],SIMULPARAMS *simulCond,PBC *box);
void image_array(int size_array,DELTA d[],SIMULPARAMS *simulCond,PBC *box);
void init_box(PBC *box);
void scale_box(PBC *box,double scale,double cell0[9]);

double kinetic(ATOM atom[],SIMULPARAMS *simulCond);
void stress_kinetic(ATOM atom[],SIMULPARAMS *simulCond,double stress[6]);
void get_kinfromtemp(ATOM atom[],SIMULPARAMS *simulCond,PBC *box);
void get_degfree(SIMULPARAMS *simulCond,PBC *box);

void nocase(char *str);
int nint(double x);

#endif
