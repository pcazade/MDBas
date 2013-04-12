/*
 * Copyright (c) 2013 Pierre-Andre Cazade
 * Copyright (c) 2013 Florent hedin
 * 
 * This file is part of MDBas.
 *
 * MDBas is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MDBas is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with MDBas.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef ERRORS_INCLUDED
#define ERRORS_INCLUDED

#ifdef	__cplusplus
extern "C" {
#endif

enum ERROR_TYPE{
  MEM_ALLOC_ERROR = 1,
  
  TOP_FILE_ERROR = 10,
  
  PSF_FILE_ERROR = 20,
  PSF_BADLINE_ERROR = 22,
  PSF_BOND_SEQ_ERROR = 23,
  PSF_ANGL_SEQ_ERROR = 24,
  PSF_DIHE_SEQ_ERROR = 25,
  PSF_IMPR_SEQ_ERROR = 26,
  PSF_UNDEF_BOND_ERROR = 71,
  PSF_UNDEF_ANGL_ERROR = 72,
  PSF_UNDEF_DIHE_ERROR = 73,
  PSF_UNDEF_IMPR_ERROR = 74,
  
  PAR_FILE_ERROR = 30,
  PAR_DIHE_FOURIER_HARMONIC_ERROR = 50,
  
  CONF_FILE_ERROR = 40,
  CONF_ATNUM_ERROR = 41,
  
  RESCONF_FILE_ERROR = 47,
  
  SIMU_FILE_ERROR = 60,
  SIMU_KEYWORD_ERROR = 61,
  SIMU_PARAM_ERROR = 62,
  SIMU_NOTFOUND_ERROR = 63,
  
  DIHE_NONPARAM_ERROR = 110,
  
  EXCLLIST_LASTATOM_ERROR = 111,
  EXCLLIST_SUMATOM_ERROR = 112,
  
  NBLIST_TOTNUM_ERROR = 120,
  
  UNKNOWN_ELEC_ERROR = 201,
  UNKNOWN_VDW_ERROR = 202,
  
  CONVERG_VEL_ERROR = 310,
  CONVERG_SHAKE_ERROR = 311,
  
  LNKCELL_CUTOFF_ERROR = 411,
  LNKCELL_NCELLS_ERROR = 412,
  
  MPI_ERROR,
          
  PLUGINS_DLOPEN_ERROR,
  PLUGINS_DLSYM_ERROR
};

/** Pointer to the output file. **/
extern FILE *outFile;

void my_error(enum ERROR_TYPE errorNumber, char file[], int line, int num_optional_args, ...);

#ifdef	__cplusplus
}
#endif

#endif