#ifndef RANDH_INCLUDED
#define RANDH_INCLUDED

#ifdef	__cplusplus
extern "C" {
#endif

    void init_rand(unsigned int seed);
    double get_rand();
    void get_BoxMuller(double *u, double *v);

#ifdef	__cplusplus
}
#endif

#endif
