#ifndef USER_H
#define	USER_H

#ifdef	__cplusplus
extern "C" {
#endif

    /* pointers for user-defined functions */
    typedef double (*UserEnergyPtr) (const PARAM *param, const PBC *box,
                                     const double x[], const double y[], const double z[],
                                     double fx[], double fy[], double fz[],
                                     const int neighList[], const int neighPair[], const int neighOrder[],
                                     const int neighList14[]);

    /* Functions from user.c ; for managing user extensions */
    UserEnergyPtr loadUserPlugin(const char pluginName[], const char funcName[]);
    void closeUserPlugin(void);

#ifdef	__cplusplus
}
#endif

#endif	/* USER_H */

