/**
 * \file errors.c
 * \brief Contains functions for managing errors.
 * \author Pierre-Andre Cazade and Florent Hedin
 * \version alpha-branch
 * \date 2012
 */

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>

#include "errors.h"

void my_error(enum ERROR_TYPE errorNumber, char file[], int line, int num_optional_args, ...)
{
  /** otpional arguments work **/
  va_list ap;
  va_start(ap,num_optional_args);
  void** optArgs = NULL;
  
  if (num_optional_args > 0)
  {
    // an array of void* will contains adress of the passed parameters
    // note that there is no way to know from those adress what are the corresponding sizes or types of the associated variables !
    optArgs=(void**)malloc(num_optional_args*sizeof(*optArgs));
    
    int i;
    for(i=0; i<num_optional_args; i++)
    {
      optArgs[i] = NULL;
      optArgs[i] = (void*)va_arg(ap,void*);
    }
  }
  
  /** **/
  
  fprintf(outFile,"MDBas failed due to error number: %d\n",errorNumber);
  fprintf(outFile,"For file %s , on line %d\n",file,line);
  
  switch (errorNumber)
  {
    
  case MEM_ALLOC_ERROR:
    fprintf(outFile,"MDBas failed with memory allocation.\n");
  break;
  
  case TOP_FILE_ERROR:
    fprintf(outFile,"MDBas cannot find or open topology file TOP.\n");
    fprintf(outFile,"Most likely, it is not properly named. Please check.\n");
  break;
  
  case PSF_FILE_ERROR:
    fprintf(outFile,"MDBas cannot find or open structure file PSF.\n");
    fprintf(outFile,"Most likely, it is not properly named. Please check.\n");
  break;
  
  case PSF_BADLINE_ERROR:
    fprintf(outFile,"MDBas encountered a problem while reading atomic properties\n");
    fprintf(outFile,"in the PSF file. There is an unexpected line there. Please\n");
    fprintf(outFile,"consult the manual for further details about PSF file\n");
  break;
  
  case PSF_BOND_SEQ_ERROR:
    fprintf(outFile,"There is problem in bonds sequence in the PSF file. Please\n");
    fprintf(outFile,"consult the manual for further details about PSF file\n");
  break;
  
  case PSF_ANGL_SEQ_ERROR:
    fprintf(outFile,"There is problem in angles sequence in the PSF file. Please\n");
    fprintf(outFile,"consult the manual for further details about PSF file\n");
  break;
  
  case PSF_DIHE_SEQ_ERROR:
    fprintf(outFile,"There is problem in dihedrals sequence in the PSF file. Please\n");
    fprintf(outFile,"consult the manual for further details about PSF file\n");
  break;
  
  case PSF_IMPR_SEQ_ERROR:
    fprintf(outFile,"There is problem in improper angles sequence in the PSF file.\n");
    fprintf(outFile,"Please consult the manual for further details about PSF file\n");
  break;
  
  case PAR_FILE_ERROR:
    fprintf(outFile,"MDBas cannot find or open parameter file PAR.\n");
    fprintf(outFile,"Most likely, it is not properly named. Please check.\n");
  break;
  
  case CONF_FILE_ERROR:
    fprintf(outFile,"MDBas cannot find or open configuration file CONF.\n");
    fprintf(outFile,"Most likely, it is not properly named. Please check.\n");
  break;
  
  case CONF_ATNUM_ERROR:
  {
    fprintf(outFile,"MDBas found a different number of atoms in CONF file and\n");
    fprintf(outFile,"in FORF file. Structure does not match configuration.\n");
    fprintf(outFile,"Check carefully these files.\n");
    
    if(num_optional_args>0)
    {
      int *n1 = (int*)optArgs[0];
      int *n2 = (int*)optArgs[1];
      fprintf(outFile,"NATOM is %d (CONF) but the value currently parsed is %d (FORF)\n",*n1,*n2);
    }
  }
  break;
  
  case RESCONF_FILE_ERROR:
    fprintf(outFile,"MDBas cannot open configuration file RESCONF.\n");
  break;
  
  case PAR_DIHE_FOURIER_HARMONIC_ERROR:
    fprintf(outFile,"A dihedral angle is specified as a Fourier series but\n");
    fprintf(outFile,"with one of the component being an harmonic potential.\n");
    fprintf(outFile,"Check in PAR file.\n");
  break;
  
  case SIMU_FILE_ERROR:
    fprintf(outFile,"MDBas cannot find or open simulation file SIMU.\n");
    fprintf(outFile,"Most likely, it is not properly named. Please check.\n");
  break;
  
  case SIMU_KEYWORD_ERROR:
  {
    fprintf(outFile,"MDBas does not recognise a keyword specified in SIMU.\n");
    fprintf(outFile,"Please check SIMU file and the manual for the list of\n");
    fprintf(outFile,"allowed keywords.\n");
    
    if(num_optional_args>0)
    {
      char *s = (char*)optArgs[0];
      fprintf(outFile,"The following keyword is unknown: '%s'\n",s);
    }
  }
  break;
  
  case SIMU_PARAM_ERROR:
  {
    fprintf(outFile,"MDBas does not recognise a parameter specified in SIMU.\n");
    fprintf(outFile,"Please check SIMU file and the manual for the list of\n");
    fprintf(outFile,"allowed keywords and their associated parameters.\n");
  
    if(num_optional_args>0)
    {
      char *s1 = (char*)optArgs[0];
      char *s2 = (char*)optArgs[1];
      fprintf(outFile,"For the following keyword '%s'; the unknown parameter is '%s'\n",s1,s2);
    }
  }
  break;
  
  case SIMU_NOTFOUND_ERROR:
  {
    fprintf(outFile,"MDBas does not find a required parameter in SIMU.\n");
    fprintf(outFile,"Please check SIMU file and the manual for the list of\n");
    fprintf(outFile,"allowed keywords and their associated parameters.\n");
  
    if(num_optional_args>0)
    {
      char *s1 = (char*)optArgs[0];
      char *s2 = (char*)optArgs[1];
      fprintf(outFile,"For the following keyword '%s'; the parameter '%s' was not found.\n",s1,s2);
    }
  }
  break;
  
  case PSF_UNDEF_BOND_ERROR:
    fprintf(outFile,"There is an undefined bond in the PSF. Most likely,\n");
    fprintf(outFile,"there are missing parameters in the PAR file. Please check\n");
  break;
  
  case PSF_UNDEF_ANGL_ERROR:
    fprintf(outFile,"There is an undefined angle in the PSF. Most likely,\n");
    fprintf(outFile,"there are missing parameters in the PAR file. Please check\n");
  break;
  
  case PSF_UNDEF_DIHE_ERROR:
    fprintf(outFile,"There is an undefined dihedral angle in the PSF. Most likely,\n");
    fprintf(outFile,"there are missing parameters in the PAR file. Please check\n");
    break;
  case PSF_UNDEF_IMPR_ERROR:
    fprintf(outFile,"There is an undefined improper angle in the PSF. Most likely,\n");
    fprintf(outFile,"there are missing parameters in the PAR file. Please check\n");
  break;
  
  case DIHE_NONPARAM_ERROR:
    fprintf(outFile,"MDBas found a too many non-parameterised dihedral angles:\n");
    fprintf(outFile,"4*nDihedrals. nDihedrals comes from the value specified\n");
    fprintf(outFile,"in PSF file. Please check in PAR file. If such a number is\n");
    fprintf(outFile,"normal for your simulation, you have to enter list.c to\n");
    fprintf(outFile,"increase the size of the 1-4 pairs array from 5*nDihedrals to\n");
    fprintf(outFile,"the size you really need. Then recompile MDBas.\n");
  break;
  
  case EXCLLIST_LASTATOM_ERROR:
    fprintf(outFile,"MDBas encountered a problem while setting the excluded atoms\n");
    fprintf(outFile,"list. The last atom has exclusion which should not happen. This\n");
    fprintf(outFile,"a bit annoying for there is no simple explanation for this.\n");
    fprintf(outFile,"Maybe an error in one of the input files which is not detected\n");
    fprintf(outFile,"by MDBas. Sorry for the trouble.\n");
  break;
  
  case EXCLLIST_SUMATOM_ERROR:
    fprintf(outFile,"MDBas encountered a problem while setting the excluded atoms\n");
    fprintf(outFile,"list. The total excluded atoms does not match of the sum of\n");
    fprintf(outFile,"excluded atoms for each atom. This a bit annoying for there is\n");
    fprintf(outFile,"no simple explanation for this. Maybe an error in one of the\n");
    fprintf(outFile,"input files which is not detected by MDBas. Sorry for the trouble.\n");
  break;
  
  case NBLIST_TOTNUM_ERROR:
    fprintf(outFile,"The number of neighbours around an atom is larger than the\n");
    fprintf(outFile,"maximum allocated memory in the verList array (default 2048).\n");
    fprintf(outFile,"This can occur fo very large cutoffs or in very heterogeneous\n");
    fprintf(outFile,"systems. This can be fixed by changing the variable MAXLIST in\n");
    fprintf(outFile,"global.h and recompilation of MDBas. An easy way to achieve this\n");
    fprintf(outFile,"without editing the header file, is to add the compilation option\n");
    fprintf(outFile,"-DMAXLIST=X, X being the new size of the array, for example 3072.\n");
  break;
  
  case UNKNOWN_ELEC_ERROR:
    fprintf(outFile,"Unknown electrostatic potential. This is most likely due to an\n");
    fprintf(outFile,"error in the SIMU file. Please check this file and the manual\n");
    fprintf(outFile,"for the list of keywords and available potentials.\n");
  break;
  
  case UNKNOWN_VDW_ERROR:
    fprintf(outFile,"Unknown van der Waals potential. This is most likely due to an\n");
    fprintf(outFile,"error in the SIMU file. Please check this file and the manual\n");
    fprintf(outFile,"for the list of keywords and available potentials.\n");
  break;
  
  case CONVERG_VEL_ERROR:
    fprintf(outFile,"Velocities quenching convergence failure, most likely due a non suitable\n");
    fprintf(outFile,"initial configuration. If not, you can try increasing the number of cycles or\n");
    fprintf(outFile,"make Shake convergence criterion more tolerant. Please check the manual.\n");
  break;
  
  case CONVERG_SHAKE_ERROR:
    fprintf(outFile,"Shake convergence failure, most likely due a non suitable initial\n");
    fprintf(outFile,"configuration. If not, you can try increasing the number of cycles or\n");
    fprintf(outFile,"make Shake convergence criterion more tolerant. Please check the manual.\n");
  break;
  
  case LNKCELL_CUTOFF_ERROR:
    fprintf(outFile,"Too large ratio of the cutoff for the link cell is asked. Maximum is 5.\n");
  break;
  
  case LNKCELL_NCELLS_ERROR:
    fprintf(outFile,"Too few link cells are created compared to the required ratio of the\n");
    fprintf(outFile,"cutoff. Please check first that you cutoff is smaller than half of the\n");
    fprintf(outFile,"smallest lattice parameter of your simulation. Alternatively, you can\n");
    fprintf(outFile,"change the ratio or force the use of the standard neighbour list algorithm.\n");
  break;
  
  case MPI_ERROR:
  {
    fprintf(outFile,"MPI error message above.\n");
  }
  break;
  
  case PLUGINS_DLOPEN_ERROR:
  {
      char *funcname = (char*)optArgs[0];
      char *plugname = (char*)optArgs[1];
      fprintf(outFile,"Error for the function '%s' : \n",funcname);
      fprintf(outFile,"Unable to open the following plugin file : '%s' \n",plugname);
  }
  break;
  
  case PLUGINS_DLSYM_ERROR:
  {
      char *funcname = (char*)optArgs[0];
      char *plugname = (char*)optArgs[1];
      char *symname  = (char*)optArgs[2];
      fprintf(outFile,"Error for the function %s : \n",funcname);
      fprintf(outFile,"Unable to find the following symbol '%s' in file '%s' \n",symname,plugname);
  }
  break;
  
  default:
    fprintf(outFile,"MDBas failed due to unknown error number: %d\n",errorNumber);
    fprintf(outFile,"Reading the manual will not help you. You are by yourself.\n");
    fprintf(outFile,"Errare humanum est.\n");
  break;
  
    
  }
  
  exit(errorNumber);
    
}
