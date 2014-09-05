#ifndef RANDH_INCLUDED
#define RANDH_INCLUDED

#ifdef	__cplusplus
extern "C" {
#endif

    void init_rand(unsigned int seed);
    real get_rand();
    void get_BoxMuller(real *u, real *v);

#ifdef	__cplusplus
}
#endif

#endif
