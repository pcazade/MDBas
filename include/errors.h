/**
 * \file errors.h
 * \brief Contains the enumeration where error codes are defined, plus prototypes of errors.c
 * \author Pierre-Andre Cazade and Florent Hedin
 * \version alpha-branch
 * \date 2012-2013
 */

#ifndef ERRORS_INCLUDED
#define ERRORS_INCLUDED

#ifdef	__cplusplus
extern "C" {
#endif

    /*!
     * \enum ERROR_TYPE
     * \brief Enumeration of the error codes used when calling the function my_error(...)
     */
    enum ERROR_TYPE {
        MEM_ALLOC_ERROR = 1, /*!< There was an error when attempting to allocate memory : used by my_malloc(), my_realloc() or my_calloc() ; this should never happen excepted if the host computer does not have enough memory ... */

        TOP_FILE_ERROR = 10, /*!< There was an error when opening the topology file ; check the name of the file and the path to it. */

        PSF_FILE_ERROR = 20, /*!< There was an error when opening the PSF file ; check the name of the file and the path to it. */
        PSF_BADLINE_ERROR = 22,
        PSF_BOND_SEQ_ERROR = 23,
        PSF_ANGL_SEQ_ERROR = 24,
        PSF_DIHE_SEQ_ERROR = 25,
        PSF_IMPR_SEQ_ERROR = 26,
        PSF_UNDEF_BOND_ERROR = 71,
        PSF_UNDEF_ANGL_ERROR = 72,
        PSF_UNDEF_DIHE_ERROR = 73,
        PSF_UNDEF_IMPR_ERROR = 74,

        PAR_FILE_ERROR = 30, /*!< There was an error when opening the PAR file ; check the name of the file and the path to it. */
        PAR_DIHE_FOURIER_HARMONIC_ERROR = 50,

        CONF_FILE_ERROR = 40, /*!< There was an error when opening the CONF file ; check the name of the file and the path to it. */
        CONF_ATNUM_ERROR = 41,

        RESCONF_FILE_ERROR = 47, /*!< There was an error when opening the restart CONF file ; check the name of the file and the path to it. */

        SIMU_FILE_ERROR = 60, /*!< There was an error when opening the main simulation file (SIMU) ; check the name of the file and the path to it. */
        SIMU_KEYWORD_ERROR = 61,
        SIMU_PARAM_ERROR = 62,
        SIMU_NOTFOUND_ERROR = 63,

        DIHE_NONPARAM_ERROR = 110,

        EXCLLIST_LASTATOM_ERROR = 111,
        EXCLLIST_SUMATOM_ERROR = 112,

        NBLIST_TOTNUM_ERROR = 120,

        UNKNOWN_ELEC_ERROR = 201, /*!< Error not yet documented during electrostatic calculations */
        UNKNOWN_VDW_ERROR = 202,  /*!< Error not yet documented during van der Walls calculations */
	
	POLAR_DIAG_ERROR = 251, /*!< There is a zero diagonal element in the polarisability tensor: singular matrix */
	POLAR_VALU_ERROR = 252, /*!< There is an illegal element in the polarisability tensor */

        CONVERG_VEL_ERROR = 310,
        CONVERG_SHAKE_ERROR = 311,

        LNKCELL_CUTOFF_ERROR = 411,
        LNKCELL_NCELLS_ERROR = 412,

        MPI_ERROR = 500, /*!< error concerning MPI parallel implementation : MPI handles by itself its own errors */

        PLUGINS_DLOPEN_ERROR = 601, /*!< There was an error when opening a plugin file */
        PLUGINS_DLSYM_ERROR = 602,  /*!< There was an error when loading a function from a plugin file */

        UNKNOWN_GENERAL_ERROR = 99999 /*!< General error to use for undocumented error (this should never happen in a perfect world !) */
    };

    /* Pointer to the output file. */

    void my_error(enum ERROR_TYPE errorNumber, char file[], int line, int num_optional_args, ...);
  
#ifdef	__cplusplus
}
#endif

#endif
