#ifndef UTILSH_INCLUDED
#define UTILSH_INCLUDED

double distance(int i,int j, ATOM *atom,double *delta,SIMULPARAMS *simulCond);
void init_vel(ATOM *atom,SIMULPARAMS *simulCond);
void image_update(ATOM *atom,SIMULPARAMS *simulCond);
double kinetic(ATOM *atom);
void get_degfree(ATOM *atom,SIMULPARAMS *simulCond);
void nocase(char *str);
int nint(double x);

#endif
