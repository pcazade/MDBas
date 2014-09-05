#ifndef USER_H
#define	USER_H

#ifdef	__cplusplus
extern "C" {
#endif

    /* pointers for user-defined functions */
    typedef real (*UserEnergyPtr) (const PARAM *param, const PBC *box,
                                     const real x[], const real y[], const real z[],
                                     real fx[], real fy[], real fz[],
                                     const int neighList[], const int neighPair[], const int neighOrder[],
                                     const int neighList14[]);

    /* Functions from user.c ; for managing user extensions */
    UserEnergyPtr loadUserPlugin(const char pluginName[], const char funcName[]);
    void closeUserPlugin(void);

#ifdef	__cplusplus
}
#endif

#endif	/* USER_H */

